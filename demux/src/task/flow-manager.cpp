#include "task/flow-manager.h"

#include <cmath>
#include <optional>
#include <utility>

#include <fmt/format.h>

#include <xoos/error/error.h>
#include <xoos/io/metadata-util.h>
#include <xoos/log/logging.h>
#include <xoos/types/str-container.h>

#include "adapter-design/adapter-component-indices.h"
#include "adapters/simplex/dual/demux-and-trim-simplex-dual-paired.h"
#include "adapters/simplex/paired/demux-and-trim-simplex-paired.h"
#include "csv.hpp"
#include "io/sample-sheet/sample-sheet.h"
#include "lut-bundle/lut-bundle.h"
#include "metrics/duplex-metrics.h"
#include "metrics/metrics-constants.h"
#include "metrics/simplex-metrics.h"
#include "task/flow-context.h"
#include "utility/stop-watch.h"

namespace xoos::demux {

// NOTE: A previous ExecutorShadow struct replicated tf::Executor's internal
// memory layout to access private thread handles for CPU affinity. It was
// removed because it is ABI-fragile — any Taskflow update that changes the
// internal layout would cause this cast to read garbage memory. If thread
// affinity is needed, it must be set from within the worker thread using
// platform APIs on the current thread. tf::Executor::this_worker_id() can be
// used only to identify the calling worker and map it to a desired CPU set;
// otherwise, request an upstream API addition.

FlowManager::FlowManager(const DemuxAndTrimParam& param, const AdapterDesign& design)
    : _adapter_type(design.type),
      _is_duplex([this]() {
        switch (_adapter_type) {
          using enum AdapterType;
          case kDuplex:
          case kDuplexUMI:
          case kDuplexStem:
            return true;
          default:
            return false;
        }
      }()),
      _param(param),
      _sid_pool(LoadSidPool(param,
                            design,
                            LutBundleRefs{.ysu = _lut_bundle_ysu,
                                          .ys = _lut_bundle_ys,
                                          .simplex = _design_bundle_simplex,
                                          .simplex_10x = _lut_bundle_simplex_10x,
                                          .duplex = _lut_bundle_duplex,
                                          .duplex_umi = _lut_bundle_duplex_umi,
                                          .duplex_stem = _lut_bundle_duplex_stem},
                            _dual_sid_matching_policy)) {
  if (_is_duplex && param.read_length_mode == ReadLengthMode::kAllSplit) {
    throw error::Error("--read-length-mode all-split is only supported for simplex adapter types");
  }

  // initialize Metrics configuration
  metrics_constraints::max_sid_id_index = static_cast<u32>(_sid_pool.size() - 1);
  metrics_constraints::max_logged_read_length = param.length_distribution_report_max;

  InitDemuxObject(design);
  _start_time = std::chrono::high_resolution_clock::now();

  // Chunk input files so that each concurrently-processed file gets worker_threads_per_input threads.
  const size_t chunk_size = std::max(1UL, param.threads / param.worker_threads_per_input);
  Logging::Info("Demultiplexing with {} threads, {} worker threads per input, {} writer threads per sample.",
                param.threads,
                param.worker_threads_per_input,
                param.writing_threads_per_sample);
  _executor = std::make_unique<tf::Executor>(param.threads);

  for (size_t i = 0; i < param.inputs.size(); i += chunk_size) {
    for (size_t j = i; j < param.inputs.size() && j < (i + chunk_size); ++j) {
      // Creating a new flow context will add a graph to the executor.
      _flow_contexts.emplace_back(std::make_unique<FlowContext>(*this, param.inputs[j], j));
      _total_file_size += std::filesystem::file_size(param.inputs[j]);
      ++_nr_files;
    }
    // While we're processing the input, TF will keep generating new tasks for every output file. Eventually,
    // those tasks will be done and the TF will shut down. We need to wait for that to happen.
    _executor->wait_for_all();
    for (const auto& flow_context : _flow_contexts) {
      _nr_sequences += flow_context->NumSequences();
      _nr_bases += flow_context->NumBases();
      _nr_bytes += flow_context->NumBytes();
    }
    // If we arrive here, all tasks are done. We can now safely remove the graphs from the executor.
    _flow_contexts.clear();
    if (Task::HasException()) {
      std::rethrow_exception(Task::GetException());
    }
    if (param.min_concord_dp_bases.has_value()) {
      auto min_concord_dp_bases = DuplexMetrics::MinConcordDupBases();

      // Check if we need to early stop
      if (param.min_concord_dp_bases.value() < min_concord_dp_bases) {
        Logging::Info("Early stopping due to minimum concordant duplex bases reached at {} bases",
                      min_concord_dp_bases);
        break;
      }
    }
  }
  _executor.reset();
  MaterializeEmptySamples();
  RenameOutputFiles();
  _end_time = std::chrono::high_resolution_clock::now();
}

void FlowManager::InitDemuxObject(const AdapterDesign& design) {
  switch (design.type) {
    using enum AdapterType;
    case kYsu: {
      _demux_and_trim_ysu = std::make_unique<DemuxAndTrimYsu>(*_lut_bundle_ysu);
      break;
    }
    case kYs: {
      _demux_and_trim_ys = std::make_unique<DemuxAndTrimYs>(*_lut_bundle_ys);
      break;
    }
    case kSimplex:
      if (_dual_sid_matching_policy.Enabled()) {
        _demux_and_trim_simplex = std::make_unique<DemuxAndTrimSimplexDualPaired>(
            _param, _design_bundle_simplex->components, _dual_sid_matching_policy);
      } else {
        _demux_and_trim_simplex =
            std::make_unique<DemuxAndTrimSimplexPaired>(_param, _design_bundle_simplex->components);
      }
      break;
    case kSimplex10x: {
      _demux_and_trim_simplex_10x = std::make_unique<DemuxAndTrimSimplex10x>(_param, *_lut_bundle_simplex_10x);
      const auto& b = *_lut_bundle_simplex_10x;
      _demux_and_trim_simplex_10x->InitFastPath(b.fixed_1.Seq(), b.sid_2_matcher.Pool(), b.fixed_3.Seq());
      break;
    }
    case kDuplex:
      _hairpin_finder = std::make_unique<HairpinFinder>(*_lut_bundle_duplex);
      _demux_and_trim_duplex = std::make_unique<DemuxAndTrimDuplex>(*_lut_bundle_duplex);
      break;
    case kDuplexUMI:
      _hairpin_finder = std::make_unique<HairpinFinder>(*_lut_bundle_duplex_umi);
      _demux_and_trim_duplex_umi = std::make_unique<DemuxAndTrimDuplexUmi>(*_lut_bundle_duplex_umi);
      break;
    case kDuplexStem:
      _hairpin_finder = std::make_unique<HairpinFinder>(*_lut_bundle_duplex_stem);
      _demux_and_trim_duplex_stem = std::make_unique<DemuxAndTrimDuplexStem>(*_lut_bundle_duplex_stem);
      break;
    default:
      throw error::Error("Unsupported adapter type: {}", static_cast<s32>(design.type));
  }
}

void FlowManager::RegisterOutputPrefix(const std::string& output_prefix, const std::string_view input_filepath) const {
  std::scoped_lock lock(_provenance_mutex);
  _output_prefix_to_input[output_prefix] = input_filepath;
}

void FlowManager::RegisterOutputFile(const fs::path& output_path, const std::string& input_name_prefix) const {
  std::scoped_lock lock(_provenance_mutex);
  _output_files[output_path.parent_path()].emplace_back(output_path, input_name_prefix);
}

namespace {

// Returns the FASTQ extension including compression suffix (e.g., ".fastq.gz", ".fastq.zst", ".fastq").
std::string GetFastqExtension(const fs::path& path) {
  const auto filename = path.filename().string();
  for (const auto* ext : {".fastq.zst", ".fastq.gz", ".fastq"}) {
    if (filename.ends_with(ext)) {
      return ext;
    }
  }
  return {};
}

/**
 * @brief Derives the sample name from the output directory structure.
 *
 * If the files are written under a mode-specific subdirectory such as "<sample>/full" or
 * "<sample>/partial", the parent directory name is the sample name. Otherwise, the current
 * directory name is used directly.
 *
 * @param[in] files      Sorted list of (output_path, input_name_prefix) pairs.
 * @param[in] sample_dir Directory containing the output files for this sample.
 * @return The derived sample name.
 */
std::string DeriveSampleName(const vec<OutputFileEntry>& files, const fs::path& sample_dir) {
  static_cast<void>(files);

  const auto dirname = sample_dir.filename().string();
  if ((dirname == kFullReadSubdir || dirname == kPartialReadSubdir) && sample_dir.has_parent_path()) {
    // current directory name is a mode-specific subdirectory, use parent directory name as sample name
    return sample_dir.parent_path().filename().string();
  }

  // return folder name as normal
  return sample_dir.filename().string();
}

}  // namespace

void FlowManager::MaterializeEmptySamples() const {
  // Only guarantee empty output when the user supplied an explicit sample sheet. Without one, the SID pool comes from
  // the adapter-design bundle's default sheet and materializing every entry would create spurious directories for
  // samples the user never asked about.
  if (!_param.sample_sheet.has_value()) {
    return;
  }

  // With no input files there is nothing to demultiplex and no input to attribute provenance to, so there is no empty
  // output to guarantee. (This also avoids dereferencing inputs.front() below.)
  if (_param.inputs.empty()) {
    return;
  }

  // A sample is considered written if any output file was registered under its directory. _output_files accumulates
  // across all inputs and chunks, so this is the global view. In split mode the primary output lives under
  // "<sample>/full"; otherwise it is directly "<sample>". Build the set of sample names that already have output.
  const bool split_mode = _param.read_length_mode == ReadLengthMode::kAllSplit;
  StrSet written_samples;
  for (const auto& [sample_dir, files] : _output_files) {
    const auto dirname = sample_dir.filename().string();
    if (split_mode && (dirname == kFullReadSubdir || dirname == kPartialReadSubdir) && sample_dir.has_parent_path()) {
      written_samples.insert(sample_dir.parent_path().filename().string());
    } else {
      written_samples.insert(dirname);
    }
  }

  // Use the first input's stem/index only to build a unique temporary prefix; the empty output is not derived from any
  // input, so its provenance input_filepath is recorded as "NA".
  const auto& first_input = _param.inputs.front();
  const auto input_index_zero = std::string{"0"};

  for (const auto& sid : _sid_pool) {
    if (written_samples.contains(sid.name)) {
      continue;
    }

    // Materialize a single valid-but-empty primary (full-read) output for this globally-empty sample. In split mode the
    // file goes under "<sample>/full" to match how a written sample is laid out; partial output is intentionally left
    // absent for an empty sample.
    const auto sample_dir = split_mode ? fmt::format("{}/{}", sid.name, kFullReadSubdir) : sid.name;
    const auto prefix = fmt::format("{}.{}_{}", sid.name, first_input.stem().string(), input_index_zero);
    RegisterOutputPrefix(prefix, kNA);

    // max_number_of_workers = 1 so the emitted filename is "<prefix>.fastq[.gz|.zst]" (the single-writer layout), which
    // RenameOutputFiles() then rewrites to "<sample>-000001.<ext>".
    Sink sink(*this,
              _param.compression_type,
              _param.out_dir,
              prefix,
              sample_dir,
              sid.id + 1,
              _param.compression_level,
              /*max_number_of_workers=*/1);
    sink.EnsureFileCreated();
    sink.FlushAndClose();
  }
}

void FlowManager::RenameOutputFiles() {
  constexpr size_t kDefaultPartIndexWidth = 6;

  // Upper bound on the number of output files per sample: each input file can produce up to
  // writing_threads_per_sample files.
  const auto max_file_count = _param.inputs.size() * std::max(_param.writing_threads_per_sample, 1UL);

  // Width is at least kDefaultPartIndexWidth, increased if needed.
  const auto part_index_width =
      std::max(kDefaultPartIndexWidth, static_cast<size_t>(std::log10(std::max(max_file_count, 1UL))) + 1);

  // Rename output files in each sample directory to the Everest naming convention
  // (<sample_name>-<part_index>.<extension>) and write a provenance TSV mapping each
  // renamed file back to its source input file.
  for (auto& [sample_dir, files] : _output_files) {
    const auto sample_name = DeriveSampleName(files, sample_dir);

    // Sort lexicographically to group files from the same input file together.
    std::ranges::sort(files);

    // Build provenance records and rename each file.
    // (new_filename, input_filepath)
    vec<std::pair<std::string, std::string>> provenance;

    // Assign contiguous 1-based part indices to each output file and rename it to the Everest format.
    for (size_t part_index = 1; part_index <= files.size(); ++part_index) {
      const auto& [old_path, input_name_prefix] = files[part_index - 1];
      const auto extension = GetFastqExtension(old_path);
      const auto new_filename = fmt::format("{}-{:0{}}{}", sample_name, part_index, part_index_width, extension);
      const auto new_path = sample_dir / new_filename;

      // Look up the input file path from the registered prefix.
      const auto it = _output_prefix_to_input.find(input_name_prefix);
      if (it == _output_prefix_to_input.end()) {
        throw error::Error("No input file mapping found for output prefix '{}'", input_name_prefix);
      }
      provenance.emplace_back(new_filename, it->second);

      if (fs::exists(new_path)) {
        if (_param.overwrite) {
          fs::remove(new_path);
        } else {
          throw error::Error(
              "Cannot rename '{}' to '{}': target already exists. "
              "Specify --overwrite to allow overwriting existing files.",
              old_path.string(),
              new_path.string());
        }
      }
      fs::rename(old_path, new_path);
    }

    // Write provenance TSV in the sample directory.
    if (!provenance.empty()) {
      WriteProvenanceTsv(sample_dir, provenance);
    }
  }
}

void FlowManager::WriteProvenanceTsv(const fs::path& sample_dir,
                                     const vec<std::pair<std::string, std::string>>& provenance) const {
  const auto tsv_path = sample_dir / kProvenanceFilename;
  if (fs::exists(tsv_path) && !_param.overwrite) {
    throw error::Error(
        "Cannot write '{}': file already exists. "
        "Specify --overwrite to allow overwriting existing files.",
        tsv_path.string());
  }
  std::ofstream tsv(tsv_path);
  io::WriteTsvMetadata(tsv, _param.command_line_info);
  auto writer = csv::make_tsv_writer(tsv);
  writer << vec<std::string>{"output_filename", "input_filepath"};
  for (const auto& [output, input] : provenance) {
    writer << vec{output, input};
  }
}

namespace {

const MatcherComponentDefinition& RequireMatcherDefinition(const AdapterDesign& design, const std::size_t index) {
  if (index >= design.components.size()) {
    throw error::Error("Adapter design '{}' has no component at position {}", design.name, index);
  }
  const auto* const matcher = std::get_if<MatcherComponentDefinition>(&design.components.at(index));
  if (matcher == nullptr) {
    throw error::Error("Adapter design '{}' component at position {} is not a matcher", design.name, index);
  }
  return *matcher;
}

bool IsDualSidSimplexDesign(const AdapterDesign& design) {
  if (design.type != AdapterType::kSimplex) {
    return false;
  }
  using Index = adapter_component_index::Simplex;
  const auto& sid_5p = RequireMatcherDefinition(design, Index::kSidFwd);
  const auto& sid_3p = RequireMatcherDefinition(design, Index::kSidRev);
  if (sid_5p.kind != MatcherKind::kSid || sid_3p.kind != MatcherKind::kSid) {
    return false;
  }
  return sid_5p.sequences.empty() || sid_3p.sequences.empty() || sid_5p.sequences != sid_3p.sequences;
}

void RequireSampleSheetForDualSidSimplex(const DemuxAndTrimParam& param, const AdapterDesign& design) {
  if (param.sample_sheet.has_value() || !IsDualSidSimplexDesign(design)) {
    return;
  }
  throw error::Error(
      "Adapter design '{}' requires --sample-sheet because dual-SID assignment needs explicit 5p-3p SID pairings",
      design.name);
}

ComponentBarcodePoolOverrides MakeSimplexDualSidOverrides(const SampleSheetSidPools& sample_sheet_pools) {
  using Index = adapter_component_index::Simplex;
  return {
      ComponentBarcodePoolOverride{Index::kSidFwd, sample_sheet_pools.sid_5p_pool},
      ComponentBarcodePoolOverride{Index::kSidRev, sample_sheet_pools.sid_3p_pool},
  };
}

std::optional<BarcodePool> MakeSingleSidOverride(const std::optional<SampleSheetSidPools>& sample_sheet_pools) {
  if (!sample_sheet_pools) {
    return std::nullopt;
  }
  return sample_sheet_pools->sid_5p_pool;
}

void ValidateSampleSheetCompatibleWithDesign(const std::optional<SampleSheetSidPools>& sample_sheet_pools,
                                             const AdapterDesign& design) {
  if (!sample_sheet_pools || !sample_sheet_pools->has_dual_sids || design.type == AdapterType::kSimplex) {
    return;
  }
  throw error::Error("Dual-SID sample sheets are currently supported only for simplex adapter designs");
}

void LogDualSidPartialMatches(const SidPool& sid_pool, const DualSidMatchingPolicy& matching_policy) {
  for (std::size_t i = 0; i < sid_pool.size() && i < matching_policy.partial_matches.size(); ++i) {
    const auto& partial_match = matching_policy.partial_matches.at(i);
    if (!partial_match.allow_5p_partial_match && !partial_match.allow_3p_partial_match) {
      continue;
    }

    std::string allowed_sides;
    if (partial_match.allow_5p_partial_match) {
      allowed_sides += "5p";
    }
    if (partial_match.allow_3p_partial_match) {
      if (!allowed_sides.empty()) {
        allowed_sides += ", ";
      }
      allowed_sides += "3p";
    }
    Logging::Info("{} ({}): allows {} SID partial match", sid_pool.at(i).name, sid_pool.at(i).sequence, allowed_sides);
  }
}

}  // namespace

SidPool FlowManager::LoadSimplexSidPool(const DemuxAndTrimParam& param,
                                        const AdapterDesign& design,
                                        const std::optional<SampleSheetSidPools>& sample_sheet_pools,
                                        const std::optional<BarcodePool>& single_sid_pool,
                                        std::unique_ptr<AdapterDesignBundle>& design_bundle_simplex,
                                        DualSidMatchingPolicy& dual_sid_matching_policy) {
  const bool has_dual_sid = IsDualSidSimplexDesign(design) || (sample_sheet_pools && sample_sheet_pools->has_dual_sids);
  if (has_dual_sid) {
    RequireSampleSheetForDualSidSimplex(param, design);
  }
  SidPool sid_pool;
  if (sample_sheet_pools && has_dual_sid) {
    design_bundle_simplex = std::make_unique<AdapterDesignBundle>(LoadAdapterDesignBundle(
        param.adapter_design_bundle, design, MakeSimplexDualSidOverrides(*sample_sheet_pools), param.threads));
    sid_pool = sample_sheet_pools->sid_pool;
    dual_sid_matching_policy = sample_sheet_pools->matching_policy;
  } else {
    design_bundle_simplex = std::make_unique<AdapterDesignBundle>(
        LoadAdapterDesignBundle(param.adapter_design_bundle, design, single_sid_pool, param.threads));
    // The assembled matcher pool may be transformed for its read-side orientation. Preserve the
    // normalized sample-sheet sequence for user-visible sample metrics and output metadata.
    sid_pool = single_sid_pool.value_or(FirstMatcherPool(design_bundle_simplex->components));
  }
  if (dual_sid_matching_policy.Enabled()) {
    LogDualSidPartialMatches(sid_pool, dual_sid_matching_policy);
  }
  return sid_pool;
}

// Allow SID pool to be specified in the constructor so it can exist a const object.
SidPool FlowManager::LoadSidPool(const DemuxAndTrimParam& param,
                                 const AdapterDesign& design,
                                 LutBundleRefs lut_bundles,
                                 DualSidMatchingPolicy& dual_sid_matching_policy) {
  const auto load_bundle_sw = StopWatch{};
  const auto sample_sheet_pools =
      param.sample_sheet ? std::make_optional(LoadSampleSheetSidPools(*param.sample_sheet, param.dual_sid_mode))
                         : std::nullopt;
  ValidateSampleSheetCompatibleWithDesign(sample_sheet_pools, design);
  auto single_sid_pool = MakeSingleSidOverride(sample_sheet_pools);
  const auto canonical_sid_pool = TryLoadCanonicalSidPool(param.adapter_design_bundle, design);
  if (!single_sid_pool.has_value() && !canonical_sid_pool.has_value() && !IsDualSidSimplexDesign(design)) {
    throw error::Error("Adapter '{}': canonical SID FASTA is unavailable and no sample sheet was supplied",
                       design.name);
  }
  if (single_sid_pool.has_value() && (!sample_sheet_pools.has_value() || !sample_sheet_pools->has_dual_sids)) {
    if (canonical_sid_pool.has_value()) {
      NormalizeSidPoolToCanonical(*single_sid_pool, *canonical_sid_pool, design);
    } else {
      Logging::Warn(
          "Adapter '{}': canonical SID FASTA is unavailable; preserving sample-sheet SIDs as supplied "
          "before applying the adapter DSL.",
          design.name);
    }
  }
  SidPool sid_pool;
  switch (design.type) {
    using enum AdapterType;
    case kYsu:
      lut_bundles.ysu = std::make_unique<LutBundleYsu>(
          LoadLutBundle<LutBundleYsu>(param.adapter_design_bundle, design, single_sid_pool, param.threads));
      // The LUT stores SID1 in the orientation required by the matcher. Metrics and output metadata
      // must instead retain the canonical sample-sheet sequence supplied to users.
      sid_pool = single_sid_pool.value_or(canonical_sid_pool.value_or(lut_bundles.ysu->Sid1Pool()));
      break;
    case kYs:
      lut_bundles.ys = std::make_unique<LutBundleYs>(
          LoadLutBundle<LutBundleYs>(param.adapter_design_bundle, design, single_sid_pool, param.threads));
      // SID2 includes the matcher transform and, for legacy designs, structural bases. Keep the
      // canonical sample sequence in user-visible metrics and output metadata.
      sid_pool = single_sid_pool.value_or(canonical_sid_pool.value_or(lut_bundles.ys->Sid2Pool()));
      break;
    case kSimplex:
      sid_pool = LoadSimplexSidPool(
          param, design, sample_sheet_pools, single_sid_pool, lut_bundles.simplex, dual_sid_matching_policy);
      break;
    case kSimplex10x:
      lut_bundles.simplex_10x = std::make_unique<LutBundleSimplex10x>(
          LoadLutBundle<LutBundleSimplex10x>(param.adapter_design_bundle, design, single_sid_pool, param.threads));
      sid_pool = lut_bundles.simplex_10x->sid_2_matcher.Pool();
      break;
    case kDuplex:
      lut_bundles.duplex = std::make_unique<LutBundleDuplex>(
          LoadLutBundle<LutBundleDuplex>(param.adapter_design_bundle, design, single_sid_pool, param.threads));
      sid_pool = lut_bundles.duplex->Sid2Pool();
      break;
    case kDuplexStem:
      lut_bundles.duplex_stem = std::make_unique<LutBundleDuplexStem>(
          LoadLutBundle<LutBundleDuplexStem>(param.adapter_design_bundle, design, single_sid_pool, param.threads));
      sid_pool = lut_bundles.duplex_stem->Sid2Pool();
      break;
    case kDuplexUMI:
      lut_bundles.duplex_umi = std::make_unique<LutBundleDuplexUmi>(
          LoadLutBundle<LutBundleDuplexUmi>(param.adapter_design_bundle, design, single_sid_pool, param.threads));
      sid_pool = lut_bundles.duplex_umi->Sid2Pool();
      break;
    default:
      throw error::Error("Unsupported adapter type: {}", static_cast<s32>(design.type));
  }
  const auto load_bundle_seconds = load_bundle_sw.ElapsedTime<std::chrono::milliseconds>();
  Logging::Info("Loaded adapter design bundle in {} ms", load_bundle_seconds);

  if (sid_pool.empty()) {
    throw error::Error("SID pool is empty after loading adapter design bundle; expected at least one SID.");
  }
  // copy elision should be happening here now
  return sid_pool;
}

FlowManager::~FlowManager() {
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(_end_time - _start_time).count();
  if (_is_duplex) {
    // Aggregate metrics for duplex
    const auto global_results = DuplexMetrics::SumTotal();
    global_results.WriteMetrics(_param, _sid_pool, _adapter_type);
    // if we have strand detection turn on log info
    if (_param.strand_detector.has_value()) {
      ReportStrandMetrics(global_results);
    }
  } else {
    const auto simplex_metrics_total = SimplexMetrics::SumTotal();
    simplex_metrics_total.WriteMetrics(_param, _sid_pool);
  }
  // Calculate Gbases per minute.
  const auto gbases = static_cast<double>(_nr_bases) * 1e-9;
  const auto gbases_per_minute = 60000.0 * gbases / static_cast<double>(duration);
  const auto mega_sequences = 1e-6 * static_cast<double>(_nr_sequences);
  const auto mega_sequences_second = 1000.0 * mega_sequences / static_cast<double>(duration);

  Logging::InfoNoExcept(
      "\nThroughput: processed {} sequences, {} bases, {} bytes in {} ms (using {} workers)\nNet throughput: {:.5} "
      "Gbases/minute, {:.5} MSequences/s",
      _nr_sequences,
      _nr_bases,
      _nr_bytes,
      duration,
      _param.threads,
      gbases_per_minute,
      mega_sequences_second);
}

const HairpinFinder& FlowManager::HairpinFinderObj() const {
  if (!_hairpin_finder) {
    throw error::Error("FlowManager::HairpinFinderObj() called for non-duplex adapter type. Saw adapter type: {}",
                       std::to_underlying(_adapter_type));
  }
  return *_hairpin_finder;
}

const DemuxAndTrimDuplex& FlowManager::DemuxObjectDuplex() const {
  switch (_adapter_type) {
    using enum AdapterType;
    case kDuplex:
      return *_demux_and_trim_duplex;
    case kDuplexStem:
      return *_demux_and_trim_duplex_stem;
    case kDuplexUMI:
      return *_demux_and_trim_duplex_umi;
    default:
      throw error::Error("FlowManager::DemuxObjectDuplex() called for non-duplex adapter type. Saw adapter type: {}",
                         std::to_underlying(_adapter_type));
  }
}

const DemuxAndTrimYsu& FlowManager::DemuxObjectYsu() const {
  if (!_demux_and_trim_ysu) {
    throw error::Error("FlowManager::DemuxObjectYsu() called but YSU demux object is not initialized");
  }
  return *_demux_and_trim_ysu;
}

const DemuxAndTrimYs& FlowManager::DemuxObjectYs() const {
  if (!_demux_and_trim_ys) {
    throw error::Error("FlowManager::DemuxObjectYs() called but YS demux object is not initialized");
  }
  return *_demux_and_trim_ys;
}

const DemuxAndTrimSimplex& FlowManager::DemuxObjectSimplex() const {
  if (!_demux_and_trim_simplex) {
    throw error::Error("FlowManager::DemuxObjectSimplex() called but simplex demux object is not initialized");
  }
  return *_demux_and_trim_simplex;
}

const DemuxAndTrimSimplex10x& FlowManager::DemuxObjectSimplex10x() const {
  if (!_demux_and_trim_simplex_10x) {
    throw error::Error("FlowManager::DemuxObjectSimplex10x() called but simplex-10x demux object is not initialized");
  }
  return *_demux_and_trim_simplex_10x;
}

}  // namespace xoos::demux
