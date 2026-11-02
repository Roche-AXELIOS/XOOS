#include "adapter-design/adapter-design-bundle.h"

#include <algorithm>
#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <variant>

#include <kseq++/seqio.hpp>

#include <taskflow/algorithm/transform.hpp>

#include <xoos/error/error.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/util/sequence-functions.h>
#include <xoos/vfs/vfs.h>

#include "adapter-design/adapter-design-constants.h"
#include "adapter-design/generate-edited-sequence.h"
#include "sequence/sequence-enums.h"
#include "utility/string-util.h"

namespace xoos::demux {
using fs::path;
using std::tuple;

namespace {
void ValidateFixedLiteralDefinition(const FixedComponentDefinition& definition);

/**
 * @brief Materialize a fixed literal straight from its definition.
 *
 * Fixed components are inline literals, so they are validated and stored without LUT/search
 * settings. Variable SID/UMI components are the only components materialized as SeqMatcher objects.
 */
FixedSequence MaterializeFixed(const FixedComponentDefinition& definition) {
  ValidateFixedLiteralDefinition(definition);
  return FixedSequence{definition.sequence};
}

const AdapterComponent& ComponentAt(const vec<AdapterComponent>& components, const std::size_t index) {
  if (index >= components.size()) {
    throw error::Error("adapter has no component at position {}", index);
  }
  return components[index];
}
}  // namespace

const FixedSequence& FixedAt(const vec<AdapterComponent>& components, const std::size_t index) {
  const AdapterComponent& component = ComponentAt(components, index);
  if (component.matcher != nullptr) {
    throw error::Error("adapter component at position {} is not a fixed literal", index);
  }
  return component.fixed;
}

SeqMatcherPtr MatcherAt(const vec<AdapterComponent>& components, const std::size_t index) {
  const AdapterComponent& component = ComponentAt(components, index);
  if (component.matcher == nullptr) {
    throw error::Error("adapter component at position {} has no matcher (it is a fixed literal)", index);
  }
  return component.matcher;
}

const BarcodePool& FirstMatcherPool(const vec<AdapterComponent>& components) {
  const auto it = std::ranges::find_if(components, [](const AdapterComponent& c) { return c.HasMatcher(); });
  if (it == components.end()) {
    throw error::Error("adapter has no variable component to source a barcode pool from");
  }
  return it->matcher->Pool();
}

SeqMatcherPtr AdapterDesignBundle::MatcherAt(const std::size_t index) const {
  return demux::MatcherAt(components, index);
}

const FixedSequence& AdapterDesignBundle::FixedAt(const std::size_t index) const {
  return demux::FixedAt(components, index);
}

/// @brief A KSeqpp type to read from the virtual filesystem
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
  MatcherComponentDefinition definition;

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
using MatcherDefinitionsBySource = std::map<BarcodePoolSource, vec<MatcherComponentDefinition>>;

namespace {
/**
 * Find the barcode-pool override, if any, targeting the component at @p index.
 */
std::optional<BarcodePool> FindOverridePool(const ComponentBarcodePoolOverrides& sid_pool_overrides,
                                            const std::size_t index) {
  const auto override_it =
      std::ranges::find_if(sid_pool_overrides, [index](const ComponentBarcodePoolOverride& sid_pool_override) {
        return sid_pool_override.component_index == index;
      });
  return override_it == sid_pool_overrides.end() ? std::optional<BarcodePool>{}
                                                 : std::optional<BarcodePool>{override_it->pool};
}

/**
 * Resolve the barcode-pool source for a matcher component: a sample-sheet override pool when one
 * targets this SID component, otherwise the FASTA path declared in the definition. This source is the
 * identity used both to load each pool once and to look a component's matcher back up by position, so
 * grouping and materialization must resolve it identically. Two SID components can share the same
 * definition (e.g. both `Sid("sid.fa")`) yet carry different override pools; disambiguating on the
 * source keeps each component wired to its own matcher.
 */
BarcodePoolSource ResolveMatcherSource(const MatcherComponentDefinition& matcher,
                                       const std::optional<BarcodePool>& override_pool,
                                       const std::size_t index) {
  if (override_pool && matcher.kind == MatcherKind::kSid) {
    return *override_pool;
  }
  if (matcher.sequences.empty()) {
    throw error::Error("Matcher component at position {} has no FASTA path or sample-sheet override", index);
  }
  return matcher.sequences;
}
}  // namespace

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

namespace {
/**
 * Load a FASTA into a barcode pool, each record in a FASTA has a sequence name which is
 * used as the name in the barcode pool, and the sequence itself.
 */
void LoadBarcodePool(const vfs::VirtualFileHandlePtr& fh, BarcodePool& pool);
u32 LoadBarcodePool(const VfsPtr& bundle, const fs::path& fasta, BarcodePool& pool);

/**
 * Determine the length of the barcode sequences in the pool, all sequences must be the same length or a
 * std::nullopt is produced.
 */
std::optional<u32> BarcodePoolLength(const BarcodePool& pool);

/// The untransformed FASTA pool used by SID matcher components in @p design, if its FASTA is available.
std::optional<BarcodePool> TryLoadCanonicalSidPool(const VfsPtr& vfs, const AdapterDesign& design) {
  for (const AdapterComponentDefinition& component : design.components) {
    if (const auto* const matcher = std::get_if<MatcherComponentDefinition>(&component);
        matcher != nullptr && matcher->kind == MatcherKind::kSid) {
      const auto fasta = vfs->Open(matcher->sequences);
      if (fasta == nullptr) {
        return std::nullopt;
      }
      BarcodePool canonical_sids;
      LoadBarcodePool(fasta, canonical_sids);
      if (!BarcodePoolLength(canonical_sids).has_value()) {
        throw error::Error("Inconsistent sequence lengths in '{}', all sequences must be the same length",
                           matcher->sequences.string());
      }
      return canonical_sids;
    }
  }
  throw error::Error("Adapter '{}': no SID matcher available to validate the supplied SID pool", design.name);
}

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
 * Determine all matcher definitions which are derived from the same source,
 * this is used to enable only loading the FASTA once. A source may be either a FASTA
 * or the provided optional SID pool.
 */
vec<tuple<BarcodePoolSource, vec<MatcherComponentDefinition>>> GroupMatcherDefinitionsBySource(
    const AdapterDesign& blueprint, const ComponentBarcodePoolOverrides& sid_pool_overrides);

/**
 * Given the @p results of the adapter design bundle loading processing, find the matcher for the
 * component that loads @p source and matches @p definition. Matching on the source as well as the
 * definition is required because two components can share an identical definition (e.g. both
 * `Sid("sid.fa")`) yet be backed by different override pools.
 */
SeqMatcherPtr FindMatcher(const vec<LoadBarcodePoolResult>& results,
                          const BarcodePoolSource& source,
                          const MatcherComponentDefinition& definition);

struct MaterializeComponentVisitor {
  AdapterComponent& component;
  const vec<LoadBarcodePoolResult>& load_barcodes;
  const std::optional<BarcodePool>& sid_pool_override;
  std::size_t component_index;

  void operator()(const FixedComponentDefinition& definition) const {
    component.fixed = MaterializeFixed(definition);
  }

  void operator()(const MatcherComponentDefinition& definition) const {
    component.bait = definition.bait;
    const auto source = ResolveMatcherSource(definition, sid_pool_override, component_index);
    component.matcher = FindMatcher(load_barcodes, source, definition);
  }
};

struct GroupMatcherDefinitionVisitor {
  MatcherDefinitionsBySource& matcher_definitions;
  const std::optional<BarcodePool>& sid_pool_override;
  std::size_t component_index;

  void operator()(const FixedComponentDefinition& definition) const {
    ValidateFixedLiteralDefinition(definition);
  }

  void operator()(const MatcherComponentDefinition& matcher) const {
    const auto source = ResolveMatcherSource(matcher, sid_pool_override, component_index);
    matcher_definitions[source].emplace_back(matcher);
  }
};

/**
 * Given the barcode pool, the edited sequence of that pool, and the definition of the barcode LUT aggreate
 * the edited barcodes into the final LUT.
 */
SeqMatcherPtr AggregateLut(const BarcodePool& barcodes,
                           u32 barcode_len,
                           const vec<vec<EditedBarcode>>& edited_barcodes,
                           const MatcherComponentDefinition& definition);

/**
 * Construct the DAG that will load the given @p fa before deriving a number of different LUT from
 * @p definitions and storing them in @p load_barcode
 */
void BuildLoadBarcodePoolGraph(tf::Taskflow& tf,
                               const VfsPtr& vfs,
                               const BarcodePoolSource& source,
                               const vec<MatcherComponentDefinition>& definitions,
                               LoadBarcodePoolResult& load_barcode_pool_result);

/**
 * Construct the DAG that will generate a LUT for a set of given barcodes, the DAG will depend on the
 * @p load_barcodes_task.
 */
void BuildGenerateLutGraph(tf::Taskflow& tf,
                           std::optional<tf::Task>& load_barcodes_task,
                           const LoadBarcodePoolResult& load_barcode_pool_result,
                           GenerateLutResult& generate_lut_result);
}  // namespace

/**
 * For each FASTA sequence file we generate the following DAG:
 *      1. Load the barcode pool from the FASTA.
 *      2. Transform the barcode pool as required according to the LUT definition, either Reverse or ReverseComplement.
 *      3. For each sequence in the pool generate a task to produce all sequences up
 *         to the defined maximum edit distance.
 *      4. For each LUT definition aggregate the edited sequences into a final LUT.
 *  This DAG is then submitted to taskflow.
 *
 *  If an @p sid_pool is provided, SID matcher definitions will be loaded from the sid_pool instead
 *  of the FASTA defined in the @p design.
 */
AdapterDesignBundle LoadAdapterDesignBundle(const VfsPtr& vfs,
                                            const AdapterDesign& design,
                                            const std::optional<BarcodePool>& sid_pool,
                                            const size_t threads) {
  auto sid_pool_overrides = ComponentBarcodePoolOverrides{};
  if (sid_pool) {
    for (std::size_t i = 0; i < design.components.size(); ++i) {
      if (const auto* const matcher = std::get_if<MatcherComponentDefinition>(&design.components.at(i));
          matcher != nullptr && matcher->kind == MatcherKind::kSid) {
        sid_pool_overrides.emplace_back(i, *sid_pool);
      }
    }
  }
  return LoadAdapterDesignBundle(vfs, design, sid_pool_overrides, threads);
}

AdapterDesignBundle LoadAdapterDesignBundle(const VfsPtr& vfs,
                                            const AdapterDesign& design,
                                            const ComponentBarcodePoolOverrides& sid_pool_overrides,
                                            const size_t threads) {
  // Affixes (prefix / suffix) arrive as literals already derived by the DSL assembler; the
  // loader consumes them directly, so no pre-pass over the design is needed here.
  const auto& resolved_design = design;

  tf::Taskflow tf;

  // First group all barcode definitions by FASTA file, this enables us to only load the
  // FASTA file once before performing any downstream transformations or LUT generation.
  vec<tuple<BarcodePoolSource, vec<MatcherComponentDefinition>>> matcher_definitions =
      GroupMatcherDefinitionsBySource(resolved_design, sid_pool_overrides);

  vec<LoadBarcodePoolResult> load_barcodes(matcher_definitions.size());
  for (size_t i = 0; i < matcher_definitions.size(); ++i) {
    const auto& [source, definitions] = matcher_definitions.at(i);
    auto& load_barcode = load_barcodes.at(i);
    load_barcode.fa = source;
    BuildLoadBarcodePoolGraph(tf, vfs, source, definitions, load_barcode);
  }

  auto executor = tf::Executor(threads);
  executor.run(tf).get();

  // Build the single read-order component vector. Variable components (kSid / kUmi) get their LUT
  // matcher; fixed literals get a materialized sequence and no matcher. The trimmers address these
  // by position.
  vec<AdapterComponent> components;
  components.reserve(resolved_design.components.size());
  for (std::size_t i = 0; i < resolved_design.components.size(); ++i) {
    AdapterComponent component;
    const auto sid_pool_override = FindOverridePool(sid_pool_overrides, i);
    std::visit(MaterializeComponentVisitor{component, load_barcodes, sid_pool_override, i},
               resolved_design.components.at(i));
    components.emplace_back(std::move(component));
  }

  return {std::move(components), resolved_design};
}

namespace {
SeqMatcherPtr FindMatcher(const vec<LoadBarcodePoolResult>& results,
                          const BarcodePoolSource& source,
                          const MatcherComponentDefinition& definition) {
  for (const auto& load_barcode : results) {
    if (load_barcode.fa != source) {
      continue;
    }
    for (const auto& generate_lut : load_barcode.generate_lut_results) {
      if (generate_lut.definition == definition) {
        return generate_lut.matcher;
      }
    }
  }

  throw error::Error("Unable to find barcode matcher: '{}'", FormatJson(AdapterComponentDefinition{definition}));
}

u32 LoadBarcodePool(const VfsPtr& bundle, const fs::path& fasta, BarcodePool& pool) {
  const auto fh = bundle->Open(fasta);
  if (fh == nullptr) {
    throw error::Error("FASTA '{}' does not exist in bundle '{}'", fasta.string(), bundle->GetName());
  }
  LoadBarcodePool(fh, pool);

  const auto length = BarcodePoolLength(pool);
  if (!length.has_value()) {
    throw error::Error("Inconsistent sequence lengths in '{}', all sequences must be the same length", fasta.string());
  }

  return *length;
}

void LoadBarcodePool(const vfs::VirtualFileHandlePtr& fh, BarcodePool& pool) {
  auto iss = KStreamInVfs{fh, &vfs::Read};
  klibpp::KSeq record;
  uint id = 0;
  while ((iss >> record).operator bool()) {
    pool.emplace_back(id, record.seq, record.name);
    ++id;
  }
}

std::optional<u32> BarcodePoolLength(const BarcodePool& pool) {
  std::unordered_set<u32> lengths;
  std::ranges::transform(
      pool, std::inserter(lengths, std::begin(lengths)), [](const auto& barcode) { return barcode.sequence.length(); });
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
    result.emplace_back(item.id, demux::Reverse(item.sequence), item.name);
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

BarcodePool Transform(const SeqTransform transform, const BarcodePool& barcodes) {
  if (transform == SeqTransform::kReverseComplement) {
    return ReverseComplement(barcodes);
  }
  if (transform == SeqTransform::kReverse) {
    return Reverse(barcodes);
  }
  return barcodes;
}

/**
 * @brief Validate that a DNA sequence contains only uppercase bases (A, C, G, T).
 *
 * @param seq         The sequence to validate.
 * @param field_name  Name of the field or source, used in the error message on failure.
 * @throws error::Error if any character in @p seq is not in kDnaAlphabet.
 */
void ValidateDnaSequence(const std::string_view seq, const std::string_view field_name) {
  for (auto c : seq) {
    if (!sequence::IsACGT(c)) {
      throw error::Error("Adapter sequence field '{}' contains invalid character '{}', expected only [{}]",
                         field_name,
                         c,
                         sequence::kDnaAlphabet);
    }
  }
}

void ValidateFixedLiteralDefinition(const FixedComponentDefinition& definition) {
  ValidateDnaSequence(definition.sequence, kSequenceKey);
}

vec<tuple<BarcodePoolSource, vec<MatcherComponentDefinition>>> GroupMatcherDefinitionsBySource(
    const AdapterDesign& blueprint, const ComponentBarcodePoolOverrides& sid_pool_overrides) {
  MatcherDefinitionsBySource matcher_definitions;
  for (std::size_t i = 0; i < blueprint.components.size(); ++i) {
    const auto sid_pool_override = FindOverridePool(sid_pool_overrides, i);
    std::visit(GroupMatcherDefinitionVisitor{matcher_definitions, sid_pool_override, i}, blueprint.components.at(i));
  }

  vec<tuple<BarcodePoolSource, vec<MatcherComponentDefinition>>> result(matcher_definitions.size());
  std::ranges::copy(matcher_definitions, std::begin(result));
  return result;
}
}  // namespace

AdapterDesignBundle LoadAdapterDesignBundle(const fs::path& bundle,
                                            const AdapterDesign& design,
                                            const std::optional<BarcodePool>& sid_pool,
                                            const size_t threads) {
  const VfsPtr vfs = vfs::Open(bundle);
  return LoadAdapterDesignBundle(vfs, design, sid_pool, threads);
}

AdapterDesignBundle LoadAdapterDesignBundle(const fs::path& bundle,
                                            const AdapterDesign& design,
                                            const ComponentBarcodePoolOverrides& sid_pool_overrides,
                                            const size_t threads) {
  const VfsPtr vfs = vfs::Open(bundle);
  return LoadAdapterDesignBundle(vfs, design, sid_pool_overrides, threads);
}

std::optional<BarcodePool> TryLoadCanonicalSidPool(const fs::path& bundle, const AdapterDesign& design) {
  return TryLoadCanonicalSidPool(vfs::Open(bundle), design);
}

namespace {
SeqMatcherPtr AggregateLut(const BarcodePool& barcodes,
                           u32 barcode_len,
                           const vec<vec<EditedBarcode>>& edited_barcodes,
                           const MatcherComponentDefinition& definition) {
  // switched to LUT indexed by 64-bit integer with only one barcode. For every
  // sequence length, we need a separate hash table. Support sequence lengths up to
  // 32. As input argument, we need the maximum length of the sequence.
  // We will now convert the "legacy" entries into an int64-based LUT.

  auto max_edit_distance{0U};
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
      int_lut.Add(static_cast<u32>(edited_barcode.sequence.length()),
                  hash_val,
                  {static_cast<u16>(edited_barcode.barcode_id), static_cast<u16>(edited_barcode.edit_distance)});
    }
  }

  auto seq_lut = std::make_shared<SeqLut>(int_lut, barcodes);
  return std::make_shared<SeqMatcher>(barcode_len,
                                      definition.max_edit_distance,
                                      definition.search.max_wiggle_left,
                                      definition.search.max_wiggle_right,
                                      seq_lut);
}

void BuildLoadBarcodePoolGraph(tf::Taskflow& tf,
                               const VfsPtr& vfs,
                               const BarcodePoolSource& source,
                               const vec<MatcherComponentDefinition>& definitions,
                               LoadBarcodePoolResult& load_barcode_pool_result) {
  // Load the barcode pool, and determine the length of the barcode.
  std::optional<tf::Task> load_barcode_pool_task;
  if (std::holds_alternative<fs::path>(source)) {
    load_barcode_pool_task = tf.emplace([&vfs, &fa = std::get<fs::path>(source), &load_barcode_pool_result]() {
      load_barcode_pool_result.barcode_len = LoadBarcodePool(vfs, fa, load_barcode_pool_result.barcodes);
    });
  } else {
    const auto& barcode_pool = std::get<BarcodePool>(source);
    const auto length = BarcodePoolLength(barcode_pool);
    if (!length.has_value()) {
      throw error::Error("Inconsistent sequence lengths in barcode pool, all sequences must be the same length");
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
}  // namespace

/**
 * @brief Prepend a prefix and/or append a suffix to each barcode in the pool.
 *
 * Both are applied before any transform (e.g. reverse_complement) so that the transform covers
 * the full sequence.
 *
 * @param pool   Barcode pool whose entries will be copied with affixes applied.
 * @param prefix DNA sequence to prepend to each barcode.
 * @param suffix DNA sequence to append to each barcode.
 * @return New BarcodePool with the affixed sequences.
 */
static BarcodePool ApplyAffixes(const BarcodePool& pool, const std::string_view prefix, const std::string_view suffix) {
  ValidateDnaSequence(prefix, kPrefixKey);
  ValidateDnaSequence(suffix, kSuffixKey);

  BarcodePool result;
  result.reserve(pool.size());
  for (const auto& barcode : pool) {
    std::string affixed_seq;
    affixed_seq.reserve(prefix.size() + barcode.sequence.size() + suffix.size());
    affixed_seq += prefix;
    affixed_seq += barcode.sequence;
    affixed_seq += suffix;
    result.emplace_back(barcode.id, std::move(affixed_seq), barcode.name);
  }
  return result;
}

namespace {
void BuildGenerateLutGraph(tf::Taskflow& tf,
                           std::optional<tf::Task>& load_barcodes_task,
                           const LoadBarcodePoolResult& load_barcode_pool_result,
                           GenerateLutResult& generate_lut_result) {
  // Apply prefix/suffix (if any) then transform the barcode pool as defined by the matcher definition.
  auto transform_barcodes = [&generate_lut_result, &load_barcode_pool_result]() {
    const auto& prefix = generate_lut_result.definition.prefix;
    const auto& suffix = generate_lut_result.definition.suffix;
    const auto& pool = load_barcode_pool_result.barcodes;
    const auto affixed = (prefix.empty() && suffix.empty()) ? pool : ApplyAffixes(pool, prefix, suffix);
    generate_lut_result.barcodes = Transform(generate_lut_result.definition.transform, affixed);

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
  auto generated_edited_sequence = [ed = generate_lut_result.definition.max_edit_distance](const Barcode& barcode) {
    vec<EditedBarcode> result;
    auto emplace_edited_barcode = [&result, barcode_id = barcode.id](u32 edit_distance, const std::string& seq) {
      result.emplace_back(barcode_id, edit_distance, seq);
    };
    GenerateEditedSequence(sequence::kDnaAlphabet, barcode.sequence, ed, emplace_edited_barcode);
    return result;
  };
  auto generate_edited_sequences_task = tf.transform(std::ref(generate_lut_result.barcodes_begin),
                                                     std::ref(generate_lut_result.barcodes_end),
                                                     std::ref(generate_lut_result.edited_barcodes_begin),
                                                     generated_edited_sequence);
  generate_edited_sequences_task.succeed(transform_barcodes_task);

  // Aggregate the barcodes together into a LUT. Use the effective barcode length
  // which includes any prefix and suffix appended above.
  auto aggregate_lut = [&load_barcode_pool_result, &generate_lut_result]() {
    const auto affix_len =
        static_cast<u32>(generate_lut_result.definition.prefix.size() + generate_lut_result.definition.suffix.size());
    const auto effective_barcode_len = load_barcode_pool_result.barcode_len + affix_len;

    // IntLut indexes _hash_tables by length, which has kHashTableSlotsByLength (32) slots.
    // With max_edit_distance, the IntLut constructor receives effective_barcode_len + max_edit_distance,
    // so the effective barcode length alone must be strictly less than 32 to leave room for edits.
    // SeqLut::HashValue also limits to 32 bases (fits in a uint64_t).
    constexpr u32 kMaxLutSequenceLength = 32;
    if (effective_barcode_len >= kMaxLutSequenceLength) {
      throw error::Error("LUT effective barcode length {} (core {} + affixes {}) exceeds maximum of {} for '{}'",
                         effective_barcode_len,
                         load_barcode_pool_result.barcode_len,
                         affix_len,
                         kMaxLutSequenceLength - 1,
                         generate_lut_result.definition.sequences.string());
    }
    generate_lut_result.matcher = AggregateLut(generate_lut_result.barcodes,
                                               effective_barcode_len,
                                               generate_lut_result.edited_barcodes,
                                               generate_lut_result.definition);
  };
  auto aggregate_lut_task = tf.emplace(aggregate_lut);
  aggregate_lut_task.succeed(generate_edited_sequences_task);
}
}  // namespace

}  // namespace xoos::demux
