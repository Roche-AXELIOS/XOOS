#include "adapter-design/adapter-design-bundle.h"

#include <xoos/error/error.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/util/sequence-functions.h>
#include <xoos/vfs/vfs.h>

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <unordered_set>

#include <kseq++/seqio.hpp>
#include <taskflow/algorithm/transform.hpp>

#include "adapter-design/generate-edited-sequence.h"
#include "utility/string-util.h"

namespace vfs = xoos::vfs;

using VfsPtr = vfs::VirtualFilesystemPtr;

using fs::path;
using std::tuple;

namespace xoos::demux {
/**
 * Finds the first matcher with a BarcodeType matching @p type, if none is found an error
 * is thrown. The @p name defines which end 5' or 3' the matcher is for and is used to generate the error.
 */
SeqMatcherPtr FindMatcher(const vec<tuple<BarcodeType, SeqMatcherPtr>>& matchers, const std::string& name,
                          const BarcodeType& type);

SeqMatcherPtr AdapterDesignBundle::GetMatcher5p(const BarcodeType& type) const {
  return FindMatcher(adapter_5p, "5'", type);
}

SeqMatcherPtr AdapterDesignBundle::GetMatcher3p(const BarcodeType& type) const {
  return FindMatcher(adapter_3p, "3'", type);
}

/// Four possible values of a DNA nucleotide for generating edited sequences
const std::string kDnaNucleotideAlphabet = "ACTG";

/// A KSeqpp type to read from the virtual filesystem
using KStreamInVfs = klibpp::KStreamIn<vfs::VirtualFileHandlePtr, decltype(&vfs::Read)>;

/**
 * Represents an edited barcode sequence, contains the id of the original barcode, the edit distance
 * of this sequence from the original, and the edited sequence itself
 */
struct EditedBarcode {
  uint barcode_id;
  u32 edit_distance;
  std::string sequence;
};

/**
 * Contains the parameters, intermediate, and final results from generating a LUT from a barcode pool
 * according to a provided definition.
 */
struct GenerateLutResult {
  BarcodeDefinition definition;

  // The transformed barcodes from which to generate the LUT
  BarcodePool barcodes;
  BarcodePool::const_iterator barcodes_begin;
  BarcodePool::const_iterator barcodes_end;

  // The edited barcodes
  vec<vec<EditedBarcode>> edited_barcodes;
  vec<vec<EditedBarcode>>::iterator edited_barcodes_begin;

  // The resulting LUT defined according to the @p definition
  SeqMatcherPtr matcher;
};

using BarcodePoolSource = std::variant<fs::path, BarcodePool>;

/**
 * Contains the FASTA from which to load the barcode pool, the resulting barcodes, their length
 * and additional information to track the generation of any derived LUT.
 */
struct LoadBarcodePoolResult {
  BarcodePoolSource fa;

  BarcodePool barcodes;
  u32 barcode_len{0};

  // All LUT derived from this barcode pool
  vec<GenerateLutResult> generate_lut_results;
};

/**
 * Load a FASTA into a barcode pool, each record in a FASTA has a sequence name which is
 * used as the name in the barcode pool, and the sequence itself.
 */
void LoadBarcodePool(const vfs::VirtualFileHandlePtr& fh, BarcodePool& pool);

/**
 * Determine the length of the barcode sequences in the pool, all sequences must be the same length or a
 * std::nullopt is produced.
 */
std::optional<u32> BarcodePoolLength(const BarcodePool& pool);

/**
 * Reverse all sequences in a barcode pool producing a new barcode pool with sequences that have the same
 * name but a reversed sequence.
 */
BarcodePool Reverse(const BarcodePool& pool);

/**
 * Reverse complement all sequences in a barcode pool producing a new barcode pool with sequences that have the same
 * name but a reverse complemented sequence.
 */
BarcodePool ReverseComplement(const BarcodePool& pool);

/**
 * Determine all barcode definition which are derived from the same source,
 * this is used to enable only loading the FASTA once. A source may be either a FASTA
 * or the provided optional SID pool.
 */
vec<tuple<BarcodePoolSource, vec<BarcodeDefinition>>> GroupBarcodeDefinitionsBySource(
    const AdapterDesign& blueprint, const std::optional<BarcodePool>& sid_pool);

/**
 * Given the @p results of the adapter design bundle loading processing, find the LUT for a given
 * barcode definition.
 */
SeqMatcherPtr FindMatcher(const vec<LoadBarcodePoolResult>& results, const BarcodeDefinition& definition);

/**
 * Given the barcode pool, the edited sequence of that pool, and the definition of the barcode LUT aggreate
 * the edited barcodes into the final LUT.
 */
SeqMatcherPtr AggregateLut(const BarcodePool& barcodes, u32 barcode_len, const vec<vec<EditedBarcode>>& edited_barcodes,
                           const LutDefinition& lut_definition, const SearchDefinition& search_definition);

/**
 * Construct the DAG that will load the given @p fa before deriving a number of different LUT from
 * @p definitions and storing them in @p load_barcode
 */
void BuildLoadBarcodePoolGraph(tf::Taskflow& tf, const VfsPtr& vfs, const BarcodePoolSource& source,
                               const vec<BarcodeDefinition>& definitions,
                               LoadBarcodePoolResult& load_barcode_pool_result);

/**
 * Construct the DAG that will generate a LUT for a set of given barcodes, the DAG will depend on the
 * @p load_barcodes_task.
 */
void BuildGenerateLutGraph(tf::Taskflow& tf, std::optional<tf::Task>& load_barcodes_task,
                           const LoadBarcodePoolResult& load_barcode_pool_result,
                           GenerateLutResult& generate_lut_result);

/**
 * For each FASTA sequence file we generate the following DAG:
 *      1. Load the barcode pool from the FASTA.
 *      2. Transform the barcode pool as required according to the LUT definition, either Reverse or ReverseComplement.
 *      3. For each sequence in the pool generate a task to produce all sequences up
 *         to the defined maximum edit distance.
 *      4. For each LUT definition aggregate the edited sequences into a final LUT.
 *  This DAG is then submitted to taskflow.
 *
 *  If an @p sid_pool is provided, the barcode definitions of type BarcodeType::kSid will be loaded from the sid_pool
 * instead of the FASTA defined in the @p design.
 */
AdapterDesignBundle LoadAdapterDesignBundle(const VfsPtr& vfs, const AdapterDesign& design,
                                            const std::optional<BarcodePool>& sid_pool, const size_t threads) {
  tf::Taskflow tf;

  // First group all barcode definitions by FASTA file, this enables us to only load the
  // FASTA file once before performing any downstream transformations or LUT generation.
  vec<tuple<BarcodePoolSource, vec<BarcodeDefinition>>> barcode_definitions =
      GroupBarcodeDefinitionsBySource(design, sid_pool);

  vec<LoadBarcodePoolResult> load_barcodes(barcode_definitions.size());
  for (size_t i = 0; i < barcode_definitions.size(); ++i) {
    const auto& [source, definitions] = barcode_definitions.at(i);
    auto& load_barcode = load_barcodes.at(i);
    load_barcode.fa = source;
    BuildLoadBarcodePoolGraph(tf, vfs, source, definitions, load_barcode);
  }

  auto executor = tf::Executor(threads);
  executor.run(tf).get();

  vec<tuple<BarcodeType, SeqMatcherPtr>> adapter_5p;
  for (const auto& item : design.adapter_5p) {
    auto matcher = FindMatcher(load_barcodes, item);
    adapter_5p.emplace_back(item.type, matcher);
  }

  vec<tuple<BarcodeType, SeqMatcherPtr>> adapter_3p;
  for (const auto& item : design.adapter_3p) {
    auto matcher = FindMatcher(load_barcodes, item);
    adapter_3p.emplace_back(item.type, matcher);
  }

  return {adapter_5p, adapter_3p};
}

SeqMatcherPtr FindMatcher(const vec<tuple<BarcodeType, SeqMatcherPtr>>& matchers, const std::string& name,
                          const BarcodeType& type) {
  auto matcher = std::find_if(std::cbegin(matchers), std::cend(matchers),
                              [type](const auto& m) { return std::get<0>(m) == type; });
  if (matcher == std::cend(matchers)) {
    // if no matcher is found debug the error by printing all matchers
    std::string all_matchers;
    for (const auto& [barcode_type, _] : matchers) {
      all_matchers += Format(barcode_type) + " ";
    }
    throw error::Error("Missing {} matcher with type '{}' matchers: {}", name, Format(type), all_matchers);
  }
  return std::get<1>(*matcher);
}

SeqMatcherPtr FindMatcher(const vec<LoadBarcodePoolResult>& results, const BarcodeDefinition& definition) {
  for (const auto& load_barcode : results) {
    for (const auto& generate_lut : load_barcode.generate_lut_results) {
      if (generate_lut.definition == definition) {
        return generate_lut.matcher;
      }
    }
  }

  throw error::Error("Unable to find barcode matcher: '{}'", FormatJson(definition));
}

u32 LoadBarcodePool(const VfsPtr& bundle, const fs::path& fasta, BarcodePool& pool) {
  auto fh = bundle->Open(fasta);
  if (fh == nullptr) {
    throw error::Error("FASTA '{}' does not exist in bundle '{}'", fasta.string(), bundle->GetName());
  }
  LoadBarcodePool(fh, pool);

  auto length = BarcodePoolLength(pool);
  if (!length) {
    throw error::Error("Inconsistent sequence lengths in '{}', all sequences must be the same length", fasta.string());
  }

  return *length;
}

void LoadBarcodePool(const vfs::VirtualFileHandlePtr& fh, BarcodePool& pool) {
  auto iss = KStreamInVfs{fh, &vfs::Read};
  klibpp::KSeq record;
  for (uint id = 0;; id++) {
    if (!(iss >> record)) {
      break;
    }
    pool.emplace_back(id, record.seq, record.name);
  }
}

std::optional<u32> BarcodePoolLength(const BarcodePool& pool) {
  std::unordered_set<u32> lengths;
  std::transform(std::cbegin(pool), std::cend(pool), std::inserter(lengths, std::begin(lengths)),
                 [](const auto& barcode) { return barcode.sequence.length(); });
  if (lengths.empty()) {
    return 0;
  }
  if (lengths.size() != 1) {
    return std::nullopt;
  }
  return *std::cbegin(lengths);
}

BarcodePool Reverse(const BarcodePool& pool) {
  BarcodePool result{};
  result.reserve(pool.size());
  for (const auto& item : pool) {
    result.emplace_back(item.id, Reverse(item.sequence), item.name);
  }
  return result;
}

BarcodePool ReverseComplement(const BarcodePool& pool) {
  BarcodePool result{};
  result.reserve(pool.size());
  for (const auto& item : pool) {
    result.emplace_back(item.id, sequence::ReverseComplementCaseSensitive(item.sequence), item.name);
  }
  return result;
}

BarcodePool Transform(LutTransform transform, const BarcodePool& barcodes) {
  if (transform == LutTransform::kReverseComplement) {
    return ReverseComplement(barcodes);
  }
  if (transform == LutTransform::kReverse) {
    return Reverse(barcodes);
  }
  return barcodes;
}

vec<tuple<BarcodePoolSource, vec<BarcodeDefinition>>> GroupBarcodeDefinitionsBySource(
    const AdapterDesign& blueprint, const std::optional<BarcodePool>& sid_pool) {
  std::map<BarcodePoolSource, vec<BarcodeDefinition>> barcode_definitions;
  for (const auto& barcodes : {blueprint.adapter_5p, blueprint.adapter_3p}) {
    for (const auto& barcode : barcodes) {
      if (sid_pool && barcode.type == BarcodeType::kSid) {
        barcode_definitions[*sid_pool].emplace_back(barcode);
      } else {
        barcode_definitions[barcode.lut.sequences].emplace_back(barcode);
      }
    }
  }

  vec<tuple<BarcodePoolSource, vec<BarcodeDefinition>>> result(barcode_definitions.size());
  std::copy(std::cbegin(barcode_definitions), std::cend(barcode_definitions), std::begin(result));
  return result;
}

AdapterDesignBundle LoadAdapterDesignBundle(const fs::path& bundle, const AdapterDesign& design,
                                            const std::optional<BarcodePool>& sid_pool, const size_t threads) {
  const VfsPtr vfs = vfs::Open(bundle);
  return LoadAdapterDesignBundle(vfs, design, sid_pool, threads);
}

SeqMatcherPtr AggregateLut(const BarcodePool& barcodes, u32 barcode_len, const vec<vec<EditedBarcode>>& edited_barcodes,
                           const LutDefinition& lut_definition, const SearchDefinition& search_definition) {
  // KZ2023 switched to LUT indexed by 64-bit integer with only one barcode. For every
  // sequence length, we need a separate hash table. Support sequence lengths up to
  // 32. As input argument, we need the maximum length of the sequence.
  // We will now convert the "legacy" entries into an int64-based LUT.

  auto max_edit_distance{0u};
  for (const auto& barcodes_entries : edited_barcodes) {
    for (const auto& b : barcodes_entries) {
      max_edit_distance = std::max(max_edit_distance, b.edit_distance);
    }
  }
  IntLut int_lut(barcode_len + max_edit_distance);

  // now convert the string-based LUT to an 64-bit bit integer LUT.
  for (const auto& seq_edited_barcodes : edited_barcodes) {
    for (const auto& edited_barcode : seq_edited_barcodes) {
      // calculate hash value and insert one barcode into output lut
      const auto hash_val{
          SeqLut::HashValue(edited_barcode.sequence, 0, static_cast<u32>(edited_barcode.sequence.length()))};
      // Note: barcode type is currently not used, so we can use the same barcode type for all barcodes

      int_lut.Add(static_cast<u32>(edited_barcode.sequence.length()), hash_val,
                  {edited_barcode.barcode_id, edited_barcode.edit_distance, BarcodeMatch::kUMI5p});
    }
  }

  auto seq_lut = std::make_shared<SeqLut>(int_lut, barcodes);
  return std::make_shared<SeqMatcher>(barcode_len, lut_definition.max_edit_distance, search_definition.max_wiggle_left,
                                      search_definition.max_wiggle_right, seq_lut);
}

void BuildLoadBarcodePoolGraph(tf::Taskflow& tf, const VfsPtr& vfs, const BarcodePoolSource& source,
                               const vec<BarcodeDefinition>& definitions,
                               LoadBarcodePoolResult& load_barcode_pool_result) {
  // Load the barcode pool, and determine the length of the barcode.
  std::optional<tf::Task> load_barcode_pool_task;
  if (std::holds_alternative<fs::path>(source)) {
    load_barcode_pool_task = tf.emplace([&vfs, &fa = std::get<fs::path>(source), &load_barcode_pool_result]() {
      load_barcode_pool_result.barcode_len = LoadBarcodePool(vfs, fa, load_barcode_pool_result.barcodes);
    });
  } else {
    const auto& barcode_pool = std::get<BarcodePool>(source);
    auto length = BarcodePoolLength(barcode_pool);
    if (!length) {
      throw error::Error("Inconsistent sequence lengths in the sample sheet, all sequences must be the same length");
    }
    load_barcode_pool_result.barcode_len = length.value();
    load_barcode_pool_result.barcodes = barcode_pool;
    load_barcode_pool_task = std::nullopt;
  }

  // For each LUT derived from this barcode pool, generate a LUT
  load_barcode_pool_result.generate_lut_results.resize(definitions.size());
  for (size_t j = 0; j < definitions.size(); ++j) {
    auto& generate_lut_result = load_barcode_pool_result.generate_lut_results.at(j);
    generate_lut_result.definition = definitions.at(j);
    BuildGenerateLutGraph(tf, load_barcode_pool_task, load_barcode_pool_result, generate_lut_result);
  }
}

void BuildGenerateLutGraph(tf::Taskflow& tf, std::optional<tf::Task>& load_barcodes_task,
                           const LoadBarcodePoolResult& load_barcode_pool_result,
                           GenerateLutResult& generate_lut_result) {
  // Transform the barcode pool as defined in the LUT definition
  auto transform_barcodes = [&generate_lut_result, &load_barcode_pool_result]() {
    generate_lut_result.barcodes =
        Transform(generate_lut_result.definition.lut.transform, load_barcode_pool_result.barcodes);

    generate_lut_result.barcodes_begin = std::cbegin(generate_lut_result.barcodes);
    generate_lut_result.barcodes_end = std::cend(generate_lut_result.barcodes);

    generate_lut_result.edited_barcodes.resize(generate_lut_result.barcodes.size());
    generate_lut_result.edited_barcodes_begin = std::begin(generate_lut_result.edited_barcodes);
  };
  auto transform_barcodes_task = tf.emplace(transform_barcodes);
  if (load_barcodes_task) {
    transform_barcodes_task.succeed(*load_barcodes_task);
  }

  // Generate all edited sequences within the max edit distance from the transformed barcodes
  auto generated_edited_sequence = [ed = generate_lut_result.definition.lut.max_edit_distance](const Barcode& barcode) {
    vec<EditedBarcode> result;
    auto emplace_edited_barcode = [&result, barcode_id = barcode.id](u32 edit_distance, const std::string& seq) {
      result.emplace_back(barcode_id, edit_distance, seq);
    };
    GenerateEditedSequence(kDnaNucleotideAlphabet, barcode.sequence, ed, emplace_edited_barcode);
    return result;
  };
  auto generate_edited_sequences_task =
      tf.transform(std::ref(generate_lut_result.barcodes_begin), std::ref(generate_lut_result.barcodes_end),
                   std::ref(generate_lut_result.edited_barcodes_begin), generated_edited_sequence);
  generate_edited_sequences_task.succeed(transform_barcodes_task);

  // Aggregate the barcodes together into a LUT
  auto aggregate_lut = [&load_barcode_pool_result, &generate_lut_result]() {
    generate_lut_result.matcher = AggregateLut(generate_lut_result.barcodes, load_barcode_pool_result.barcode_len,
                                               generate_lut_result.edited_barcodes, generate_lut_result.definition.lut,
                                               generate_lut_result.definition.search);
  };
  auto aggregate_lut_task = tf.emplace(aggregate_lut);
  aggregate_lut_task.succeed(generate_edited_sequences_task);
}
}  // namespace xoos::demux
