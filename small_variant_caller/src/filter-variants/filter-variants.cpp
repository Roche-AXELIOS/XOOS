#include "filter-variants.h"

#include <atomic>
#include <filesystem>
#include <map>
#include <memory>
#include <ranges>
#include <thread>
#include <utility>

#include <fmt/ranges.h>
#include <htslib/bgzf.h>
#include <htslib/vcf.h>

#include <taskflow/taskflow.hpp>

#include <xoos/enum/enum-util.h>
#include <xoos/error/error.h>
#include <xoos/io/bed-region.h>
#include <xoos/io/fasta-reader.h>
#include <xoos/io/metadata-util.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/io/vcf/vcf-writer.h>
#include <xoos/log/logging.h>
#include <xoos/sex_predict/sex.h>
#include <xoos/sex_predict/vcf-sex-predictor.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>
#include <xoos/util/container-functions.h>
#include <xoos/util/parse-int.h>
#include <xoos/util/string-functions.h>
#include <xoos/util/tmp-dir.h>

#include "compute-bam-features/compute-bam-features.h"
#include "compute-bam-features/progress-meter.h"
#include "compute-vcf-features/compute-vcf-features.h"
#include "compute-vcf-features/vcf-header-util.h"
#include "core/cli-option-names.h"
#include "core/filtering.h"
#include "core/model-metadata.h"
#include "core/vcf-fields.h"
#include "file_type_inference/detect-sample-type.h"
#include "mitochondria-processing.h"
#include "shap-value-tsv.h"
#include "util/bgzf-utils.h"
#include "util/cli-util.h"
#include "util/file-util.h"
#include "util/log-util.h"
#include "util/model-resolver.h"
#include "util/parallel-compute-utils.h"
#include "vcf-metrics.h"

namespace xoos::svc {

using BedRegion = io::BedRegion;

void ReplacePopafInfoLine(const io::VcfHeaderPtr& header) {
  header->RemoveInfoLine(kFieldPopaf);
  header->AddInfoLine(
      {.id = kFieldPopaf, .description = kPopafDesc, .number = io::kNumberOne, .type = io::FieldType::kFloat});
}

const StrUnorderedMap<u32>& SelectSexPredictionMedianDp(const ChromMedianDepth& chrom_median_dp,
                                                        const Workflow workflow) {
  return workflow == Workflow::kTumorNormalWgs ? chrom_median_dp.tumor : chrom_median_dp.normal;
}

namespace {
/**
 * Loads a block list from file
 * @param block_list The blocklist in chr_pos_ref_alt format
 * @return The parsed left-padded blocklist
 */
StrUnorderedSet LoadBlockList(const fs::path& block_list) {
  StrUnorderedSet result;
  std::ifstream block_list_stream(block_list);
  std::string line;
  while (getline(block_list_stream, line)) {
    // Need to left pad the blocklist so it matches the behavior of GetVariantCorrelationKey
    vec<std::string> split_line;
    std::string seg;
    std::stringstream input_line(line);
    while (std::getline(input_line, seg, '_')) {
      split_line.push_back(seg);
    }
    constexpr size_t kBlocklistFieldCount = 4;
    if (split_line.size() == kBlocklistFieldCount) {
      const std::string key =
          GetVariantCorrelationKey(split_line[0], util::ParseU64(split_line[1]) - 1, split_line[2], split_line[3]);
      result.insert(key);
    }
  }
  return result;
}

/**
 * @brief Helper function to extract sex chromosome name and PARs from a BED file.
 * @param bed_path Path of BED file containing PARs
 * @param default_chrom Default chromosome name if not specified in BED file
 * @return Chromosome name and vector of PARs.
 */
std::pair<std::string, vec<Interval>> ExtractSexChromPAR(const std::optional<fs::path>& bed_path,
                                                         const std::string& default_chrom) {
  const auto chromosome_label = default_chrom.starts_with("chr") ? default_chrom.substr(3) : default_chrom;
  std::string name;
  vec<Interval> par;
  if (bed_path.has_value()) {
    auto chr_to_intervals = GetChromIntervalMap(bed_path).value();
    if (chr_to_intervals.empty()) {
      WarnAsErrorIfSet("No intervals found in {}", bed_path.value());
    } else if (chr_to_intervals.size() == 1) {
      name = chr_to_intervals.begin()->first;
      par = chr_to_intervals.begin()->second;
    } else {
      vec<std::string> ref_names;
      ref_names.reserve(chr_to_intervals.size());
      for (const auto& [chrom, intervals] : chr_to_intervals) {
        ref_names.emplace_back(chrom);
      }
      WarnAsErrorIfSet("PARs BED file has multiple reference names: {}", string::Join(ref_names, ","));
    }
  }
  if (par.empty()) {
    WarnAsErrorIfSet("No PARs available for chromosome {}", chromosome_label);
  }
  if (name.empty()) {
    name = default_chrom;
    WarnAsErrorIfSet("Using default name for chromosome {}: {}", chromosome_label, name);
  }
  return std::make_pair(name, par);
}

/**
 * @brief Decide whether tumor WGS sample-type detection must be deferred to VCF pre-scan.
 * @details The CLI parse phase cannot infer sample type without reading the VCF. For tumor WGS workflows,
 * BWA, auto-detect, and no user-provided model, this helper enables collection of sample-type metrics in
 * the pre-scan so the VCF is not parsed in a separate detector pass.
 * @param param Parsed filter_variants parameters.
 * @return True when pre-scan should collect sample-type metrics and resolve the model afterward.
 */
bool RequiresDeferredSampleTypeDetection(const FilterVariantsParam& param) {
  const bool supports_detection =
      param.workflow == Workflow::kTumorNormalWgs || param.workflow == Workflow::kTumorOnlyWgs;
  return supports_detection && param.sample_type == SampleType::kAutoDetect && param.aligner == AlignerType::kBwa &&
         !param.model_provided;
}

/**
 * @brief Determine the VCF sample index used for sample-type detection.
 * @details Tumor-only VCFs always use sample index 0. Tumor-normal VCFs use `tumor_sample` metadata when
 * available and otherwise fall back to sample index 0.
 * @param workflow Active filtering workflow.
 * @param header VCF header to inspect.
 * @param vcf_file VCF path used in warning messages.
 * @return Zero-based VCF sample index to evaluate.
 */
s32 DetermineTumorSampleIndex(const Workflow workflow, const io::VcfHeaderPtr& header, const fs::path& vcf_file) {
  if (workflow == Workflow::kTumorOnlyWgs) {
    return 0;
  }

  const auto sample_indexes = GetTumorNormalSampleIndexes(header);
  if (sample_indexes.has_value()) {
    return sample_indexes->tumor_sample_idx;
  }

  Logging::Warn("VCF header '{}' does not define tumor_sample; using sample index 0 for sample-type detection",
                vcf_file);
  return 0;
}

/**
 * @brief Compute the denominator length for sample-type detection in filter_variants.
 * @details BED targets take precedence because only target-overlapping records are evaluated. Without BED
 * targets, the function sums VCF header contig lengths and falls back to the default whole-genome length
 * when header lengths are absent or sum to zero.
 * @param header VCF header containing contig length metadata.
 * @param target_regions Optional BED intervals used for both filtering and denominator length.
 * @return Total length to use in the sample-type phred denominator.
 */
u64 ComputeSampleTypeDetectionLength(const io::VcfHeaderPtr& header,
                                     const std::optional<ChromIntervalsMap>& target_regions) {
  if (target_regions.has_value()) {
    // BED-restricted detection must use BED length as the denominator because records outside the BED do
    // not contribute observations.
    return SumChromIntervalLengths(*target_regions);
  }

  const u64 header_length = SumVcfHeaderContigLengths(header);
  return header_length > 0 ? header_length : kDefaultWholeGenomeLengthForSampleTypeDetection;
}

/**
 * @brief Apply model-specific thresholds after deferred sample-type resolution.
 * @details When the model was resolved from sample type and aligner, threshold defaults may also need to
 * change. Explicit user-provided threshold options are preserved; only unset values are replaced.
 * @param param filter_variants parameters updated in place.
 */
void ApplyResolvedModelThresholds(FilterVariantsParam& param) {
  if (param.model_provided || param.aligner == AlignerType::kCustom) {
    return;
  }

  const auto thresholds = ResolveModelThresholds(param.sample_type, param.aligner, param.config.model_thresholds);
  if (!thresholds.has_value()) {
    return;
  }
  if (!param.snv_min_ml_score_provided && thresholds->snv_min_ml_score.has_value()) {
    param.snv_min_ml_score = *thresholds->snv_min_ml_score;
  }
  if (!param.indel_min_ml_score_provided && thresholds->indel_min_ml_score.has_value()) {
    param.indel_min_ml_score = *thresholds->indel_min_ml_score;
  }
  if (param.workflow == Workflow::kTumorOnlyWgs && !param.min_af_provided && thresholds->min_af.has_value()) {
    param.min_af = *thresholds->min_af;
  }
  if (param.workflow == Workflow::kTumorNormalWgs && !param.min_tumor_af_provided &&
      thresholds->min_tumor_af.has_value()) {
    param.min_tumor_af = *thresholds->min_tumor_af;
  }
}

/**
 * @brief Synchronize rendered tumor WGS command-line metadata after deferred resolution.
 * @details Updates sample type, aligner-derived model, and threshold options in the command line stored in
 * FilterVariantsParam. This keeps output VCF metadata aligned with the runtime configuration.
 * @param param filter_variants parameters whose command_line field is updated in place.
 */
void SyncDeferredModelCommandLine(FilterVariantsParam& param) {
  if (!param.command_line.has_value()) {
    return;
  }

  auto& command_line = param.command_line->command_line;
  ReplaceOrAppendCommandLineOption(command_line, cli_opt_name::kAligner, enum_util::FormatEnumName(param.aligner));
  ReplaceOrAppendCommandLineOption(
      command_line, cli_opt_name::kSampleType, enum_util::FormatEnumName(param.sample_type));
  ReplaceOrAppendCommandLineOption(command_line, cli_opt_name::kModel, param.model.string());
  ReplaceOrAppendCommandLineOption(
      command_line, cli_opt_name::kSnvMinMlScore, fmt::format("{}", param.snv_min_ml_score));
  ReplaceOrAppendCommandLineOption(
      command_line, cli_opt_name::kIndelMinMlScore, fmt::format("{}", param.indel_min_ml_score));
  if (param.workflow == Workflow::kTumorOnlyWgs) {
    ReplaceOrAppendCommandLineOption(command_line, cli_opt_name::kMinAf, fmt::format("{}", param.min_af));
  } else {
    ReplaceOrAppendCommandLineOption(command_line, cli_opt_name::kMinTumorAf, fmt::format("{}", param.min_tumor_af));
  }
  Logging::Info("Resolved command line: {}", command_line);
}

/**
 * @brief Resolve sample type, model, thresholds, and command-line metadata after VCF pre-scan.
 * @details Consumes sample-type metrics collected during ParallelVcfPreScan(). If detection was not
 * enabled, the function is a no-op. When enabled, it classifies the sample type, resolves the corresponding
 * workflow-specific model, reapplies model thresholds while respecting explicit CLI values, and updates output
 * metadata.
 * @param param filter_variants parameters updated in place.
 * @param pre_scan Pre-scan result that must contain sample-type metrics when detection is enabled.
 * @param sample_type_detection_enabled Whether deferred detection was requested for this invocation.
 * @throws xoos::error::Error if detection is enabled but pre-scan metrics are missing.
 */
void ResolveDeferredSampleType(FilterVariantsParam& param,
                               const VcfPreScanResult& pre_scan,
                               const bool sample_type_detection_enabled) {
  if (!sample_type_detection_enabled) {
    return;
  }
  if (!pre_scan.sample_type_metrics.has_value()) {
    throw error::Error("Internal error: sample-type metrics were not collected during VCF pre-scan");
  }

  const auto& metrics = *pre_scan.sample_type_metrics;
  param.sample_type = ClassifySampleType(metrics.phred, param.sample_type_max_ffpe_phred);
  Logging::Info(
      "Inferred sample type '{}' from VCF phred {:.4f}, max FFPE phred {}, median DP {}, error count {}, "
      "total length {}",
      enum_util::FormatEnumName(param.sample_type),
      metrics.phred,
      param.sample_type_max_ffpe_phred,
      metrics.median_dp,
      metrics.error_count,
      metrics.total_length);
  if (param.workflow == Workflow::kTumorOnlyWgs) {
    param.model = ResolveTumorOnlyWgsModel(param.sample_type, param.aligner);
  } else {
    param.model = ResolveTumorNormalModel(param.sample_type, param.aligner);
  }
  // Threshold resolution happens here rather than in CLI parsing because the model is unknown until the
  // sample type is inferred from pre-scan metrics.
  ApplyResolvedModelThresholds(param);
  SyncDeferredModelCommandLine(param);
}

/**
 * @brief Extract sex chromosome information and predicted sex using pre-computed chromosomal median DP values.
 * @param param CLI parameters for `filter_variants`
 * @param chr_med_dp Pre-computed chromosome median depth values from the VCF pre-scan
 * @return Tuple: predicted sex, chrX name, chrY name, chrX PAR intervals, chrY PAR intervals, chromosomal median DPs
 */
std::tuple<sex_predict::Sex, std::string, std::string, vec<Interval>, vec<Interval>, ChromMedianDepth> GetSexInfo(
    const FilterVariantsParam& param, const ChromMedianDepth& chr_med_dp) {
  using enum sex_predict::Sex;
  // Extract the name and PARs for chromosomes X and Y
  auto [chr_x_name, chr_x_par] = ExtractSexChromPAR(param.par_bed_x, kDefaultChrXName);
  auto [chr_y_name, chr_y_par] = ExtractSexChromPAR(param.par_bed_y, kDefaultChrYName);
  auto sex = kUnknown;

  // Tumor-normal workflows predict biological sex from the tumor sample; all other workflows use the
  // normal or single-sample DP values. Prediction depends only on the per-chromosome median DP, the
  // sex-determination autosome name, and the chrX name; PAR BED files only drive downstream calling.
  const auto& sex_prediction_median_dp = SelectSexPredictionMedianDp(chr_med_dp, param.workflow);
  if (sex_prediction_median_dp.empty()) {
    // Not fatal on its own: median DP is only required when feature normalization is enabled, and
    // `ValidateChromMedianDepthForNormalization` reports that case with an actionable message.
    Logging::Warn("No median DP values extracted; biological sex cannot be predicted");
  } else if (!chr_x_name.empty()) {
    const auto get_median_dp = [&sex_prediction_median_dp](const std::string& chrom) {
      const auto itr = sex_prediction_median_dp.find(chrom);
      return itr == sex_prediction_median_dp.end() ? 0U : itr->second;
    };
    Logging::Info("Chromosome median DP: {}={}, {}={}",
                  param.sd_chr_name,
                  get_median_dp(param.sd_chr_name),
                  chr_x_name,
                  get_median_dp(chr_x_name));
    sex = sex_predict::PredictSex(sex_prediction_median_dp, param.sd_chr_name, chr_x_name);
    Logging::Info("Sex: {}", GetDescriptionForSex(sex));
    if (sex == kUnknown) {
      WarnAsErrorIfSet("Cannot determine the biological sex of input sample");
    }
  }
  return std::make_tuple(sex, chr_x_name, chr_y_name, chr_x_par, chr_y_par, chr_med_dp);
}
}  // namespace

void FilterVariantsClass::SetGlobalContext() {
  using enum Workflow;
  // The GlobalContext struct contains various parameters required for workflow-specific filtering, and they do not
  // change between regions. The exact parameters may differ for each workflow, but certain parameters are only
  // applicable to specific workflows.

  _global_ctx.model_config = _model_config;

  // BAM feature extraction parameters
  _global_ctx.bam_feat_params =
      ComputeBamFeaturesParams{.feature_cols = _model_config.feature_cols,
                               .min_bq = _param.min_bq,
                               .min_mapq = _param.min_mapq,
                               .min_allowed_distance_from_end = _param.min_allowed_distance_from_end,
                               .min_family_size = _min_family_size,
                               .max_read_variant_count = _param.max_read_variant_count,
                               .max_read_variant_count_normalized = _param.max_read_variant_count_normalized,
                               .min_homopolymer_length = _param.min_homopolymer_length,
                               .sequencing_protocol = _param.sequencing_protocol,
                               .filter_homopolymer = _param.filter_homopolymer,
                               .tumor_sample_name = _param.tumor_sample_name,
                               .tumor_rg_ids = GetReadGroupIdsForSample(
                                   _alignment_reader_cache.Open(_param.bam_files, "r"), _param.tumor_sample_name),
                               .decode_yc = _param.decode_yc,
                               .missing_yc_tag_warning_emitted = std::make_shared<std::atomic_bool>(false),
                               .min_base_type = _param.min_base_type,
                               .duplex_lowbq = _param.duplex_lowbq};
  _global_ctx.skip_variants_vcf = _param.skip_variants_vcf;

  // reference genome and region parameters
  _global_ctx.genome = _param.genome;
  _global_ctx.ref_seqs = &_ref_seqs;
  _global_ctx.bed_regions =
      _param.bed_file.has_value() ? GetChromIntervalMap(_param.bed_file).value() : ChromIntervalsMap{};
  _global_ctx.interest_regions =
      _param.interest_regions.has_value() ? _param.interest_regions.value() : ChromIntervalsMap{};

  _global_ctx.hdr = _hdr;
  std::tie(_global_ctx.vcf_info_metadata, _global_ctx.vcf_fmt_metadata) = _hdr->GetFieldMetadata();

  // Sex chromosome and normalization parameters — uses pre-computed median DP from the VCF pre-scan
  const auto& [sex, chr_x_name, chr_y_name, chr_x_par, chr_y_par, chr_med_dp] =
      GetSexInfo(_param, _pre_scan.chrom_median_dp);
  if (_param.normalize_features == FeatureNormalization::kMedianDp) {
    // Check whether the calculated median depth values are valid for normalization.
    ValidateChromMedianDepthForNormalization(chr_med_dp, FeatureNamesHaveSampleContext(_model_config.workflow));
  }
  // Store chromosome median DP when feature normalization is enabled or when the normalized DP ratio filter is active,
  // since the latter needs per-chromosome median DP values to compute the ratio.
  const bool needs_chrom_median_dp =
      _param.normalize_features == FeatureNormalization::kMedianDp || _param.min_dp_ratio > 0;
  _global_ctx.normalize_targets = needs_chrom_median_dp ? chr_med_dp : ChromMedianDepth{};

  // Pre-computed variant density from the VCF pre-scan (avoids boundary effects from per-region computation)
  _global_ctx.pre_scan_variants = &_pre_scan.chrom_variants;
  _global_ctx.pre_scan_density = &_pre_scan.chrom_variant_density;
  _global_ctx.normalize_scoring_features = (_param.normalize_features == FeatureNormalization::kMedianDp);
  _global_ctx.sex = sex;
  _global_ctx.chr_x_name = chr_x_name;
  _global_ctx.chr_y_name = chr_y_name;
  _global_ctx.chr_x_par = chr_x_par;
  _global_ctx.chr_y_par = chr_y_par;

  // Workflow specific parameters
  SetWorkflowSpecificContext();
}

void FilterVariantsClass::SetWorkflowSpecificContext() {
  using enum Workflow;
  switch (_model_config.workflow) {
    case kGermlineWgs: {
      break;
    }
    case kTumorOnlyTe: {
      StrUnorderedSet block_list;
      if (_param.block_list) {
        block_list = LoadBlockList(*_param.block_list);
      }
      StrUnorderedSet hotspots;
      if (_param.hotspot_list) {
        hotspots = LoadHotspotVariants(*_param.hotspot_list);
      }
      _global_ctx.phased = _param.phased;
      _global_ctx.force_calls = GetChromIntervalMap(_param.forcecall_list);
      _global_ctx.hotspots = hotspots;
      _global_ctx.block_list = block_list;
      _global_ctx.min_af = _param.min_af;
      _global_ctx.weighted_counts_threshold = _param.weighted_counts_threshold;
      _global_ctx.hotspot_weighted_counts_threshold = _param.hotspot_weighted_counts_threshold;
      _global_ctx.ml_threshold = _param.ml_threshold;
      _global_ctx.hotspot_ml_threshold = _param.hotspot_ml_threshold;
      _global_ctx.min_phased_allele_freq = _param.min_phased_allele_freq;
      _global_ctx.max_phased_allele_freq = _param.max_phased_allele_freq;
      _global_ctx.min_alt_counts = _param.min_alt_counts;
      break;
    }
    case kGermlineTagging: {
      const auto tn_sample_idx = GetTumorNormalSampleIndexes(_hdr);
      if (tn_sample_idx.has_value()) {
        _global_ctx.vcf_tumor_index = tn_sample_idx->tumor_sample_idx;
        _global_ctx.vcf_normal_index = tn_sample_idx->normal_sample_idx;
      }
      _global_ctx.snv_min_ml_score = _param.snv_min_ml_score;
      _global_ctx.indel_min_ml_score = _param.indel_min_ml_score;
      _global_ctx.is_germline_tagging = true;
      break;
    }
    case kTumorNormalWgs: {
      _global_ctx.snv_min_ml_score = _param.snv_min_ml_score;
      _global_ctx.indel_min_ml_score = _param.indel_min_ml_score;
      _global_ctx.min_tumor_support = _param.min_tumor_support;
      _global_ctx.max_normal_support = _param.max_normal_support;
      _global_ctx.min_tumor_af = _param.min_tumor_af;
      // Relaxed catalog thresholds are only meaningful when a catalog VCF is supplied.
      if (_param.catalog_vcf.has_value()) {
        _global_ctx.catalog_min_tumor_af = _param.catalog_min_tumor_af;
        _global_ctx.catalog_snv_min_ml_score = _param.catalog_snv_min_ml_score;
        _global_ctx.catalog_indel_min_ml_score = _param.catalog_indel_min_ml_score;
      }
      _global_ctx.min_dp_ratio = _param.min_dp_ratio;
      _global_ctx.max_indel_size = _param.max_indel_size;
      const auto tn_sample_idx = GetTumorNormalSampleIndexes(_hdr);
      if (tn_sample_idx.has_value()) {
        _global_ctx.vcf_tumor_index = tn_sample_idx->tumor_sample_idx;
        _global_ctx.vcf_normal_index = tn_sample_idx->normal_sample_idx;
      }
      break;
    }
    case kTumorOnlyWgs: {
      _global_ctx.snv_min_ml_score = _param.snv_min_ml_score;
      _global_ctx.indel_min_ml_score = _param.indel_min_ml_score;
      _global_ctx.min_support = _param.min_support;
      _global_ctx.min_af = _param.min_af;
      _global_ctx.max_indel_size = _param.max_indel_size;
      break;
    }
    case kGermlineMitochondria: {
      _global_ctx.min_af = _param.min_af;
      _global_ctx.min_dp = _param.min_dp;
      _global_ctx.min_hom_alt_af = _param.min_hom_alt_af;
      _global_ctx.blocklist_regions =
          _param.blocklist_bed.has_value() ? GetChromIntervalMap(_param.blocklist_bed).value() : ChromIntervalsMap{};
      break;
    }
    default: {
      break;
    }
  }
}

namespace {
/**
 * @brief Helper function to create score calculators based on workflow and model configuration.
 * @param param CLI parameters for `filter_variants`
 * @param model_config Model configuration for scoring
 * @return Vector of ScoreCalculator instances
 */
vec<ScoreCalculator> CreateScoreCalculators(const FilterVariantsParam& param, const SVCConfig& model_config) {
  using enum Workflow;
  vec<ScoreCalculator> calculators;
  // The germline-mitochondria workflow filters on allele frequency and depth only; it has no ML model.
  if (param.workflow == kGermlineMitochondria) {
    return calculators;
  }
  if (param.workflow == kGermlineWgs) {
    calculators.emplace_back(param.snv_model,
                             model_config.snv_scoring_cols.size(),
                             model_config.snv_model_lgbm_prediction_params,
                             param.snv_shap_value_tsv.has_value());
    calculators.emplace_back(param.indel_model,
                             model_config.indel_scoring_cols.size(),
                             model_config.indel_model_lgbm_prediction_params,
                             param.indel_shap_value_tsv.has_value());
  } else {
    calculators.emplace_back(param.model,
                             model_config.scoring_cols.size(),
                             model_config.model_lgbm_prediction_params,
                             param.shap_value_tsv.has_value());
  }
  return calculators;
}

/**
 * @brief Run a single filter task for one genomic region.
 *
 * Dispatches to the appropriate workflow-specific filter method, collects the resulting
 * VCF records and SHAP rows, and pushes them to the writer queue. If the queue has been
 * shut down (e.g., due to a writer failure), the task returns immediately.
 *
 * @param worker Per-thread filter worker that owns the readers and calculators.
 * @param workflow Active workflow type (germline, tumor-only, tumor-normal, etc.).
 * @param region Genomic region to filter.
 * @param region_index Index of this region in the partitioned region list, used to
 *                     associate the result with the correct output file.
 * @param queue Thread-safe queue to push the completed result to.
 * @param progress Shared progress tracker for logging.
 * @throws error::Error if filtering fails for the region.
 */
void RunFilterTask(FilterRegionClass& worker,
                   const Workflow workflow,
                   const TargetRegion& region,
                   const size_t region_index,
                   RegionOutputQueue& queue,
                   Progress& progress) {
  using enum Workflow;
  if (queue.IsShutdown()) {
    return;
  }
  try {
    progress.UpdateAndLog(log::LogLevel::kInfo);
    RegionResult result;
    switch (workflow) {
      case kGermlineWgs: {
        worker.FilterGermlineRegion(region, result);
        break;
      }
      case kGermlineTagging: {
        worker.FilterGermlineTaggingRegion(region, result);
        break;
      }
      case kTumorOnlyTe: {
        worker.FilterTumorOnlyTeRegion(region, result);
        break;
      }
      case kTumorNormalWgs: {
        worker.FilterTumorNormalRegion(region, result);
        break;
      }
      case kTumorOnlyWgs: {
        worker.FilterTumorOnlyWgsRegion(region, result);
        break;
      }
      case kGermlineMitochondria: {
        worker.FilterMitochondriaRegion(region, result);
        break;
      }
      default: {
        break;
      }
    }
    queue.Push({.region_index = region_index, .result = std::move(result)});
  } catch (std::exception& e) {
    throw error::Error(
        "Error filtering variants in region '{}:{}-{}': {}", region.chrom, region.start, region.end, e.what());
  }
}

/**
 * @brief Writer thread body: consume completed regions and write each to a headerless BGZF file.
 *
 * Pops RegionOutput items from the queue, writes VCF records as raw BGZF data (no VCF header)
 * using vcf_format + bgzf_write, and appends any SHAP value rows to the corresponding TSV writers.
 * On failure, captures the exception and shuts down the queue to signal producer tasks to stop.
 *
 * @param num_regions Total number of regions expected. The function returns after consuming
 *                    this many items or when the queue is shut down.
 * @param queue Thread-safe queue to pop completed region results from.
 * @param region_files Ordered list of per-region output file paths, indexed by region_index.
 * @param hdr VCF header used to format records via VcfWriter::WriteBgzfRecord.
 * @param shap_writers Optional SHAP value TSV writers (combined, SNV-only, InDel-only).
 * @param exception_out Output parameter set to the caught exception if the writer fails.
 */
void WriteRegionResults(const size_t num_regions,
                        RegionOutputQueue& queue,
                        const vec<fs::path>& region_files,
                        const io::VcfHeaderPtr& hdr,
                        const ShapWriters& shap_writers,
                        std::exception_ptr& exception_out) {
  try {
    size_t regions_written = 0;
    while (regions_written < num_regions) {
      auto item = queue.Pop();
      if (!item.has_value()) {
        break;
      }
      const auto& output = item.value();
      const auto& file_path = region_files.at(output.region_index);

      // Write this region's records as raw BGZF data (no VCF header).
      // vcf_format + bgzf_write bypasses bcf_write, which requires bcf_hdr_write first.
      if (!output.result.out_records.empty()) {
        const BgzfPtr bgzf(bgzf_open(file_path.c_str(), "w"));
        if (bgzf == nullptr) {
          throw error::Error("Failed to open BGZF file: {}", file_path.string());
        }
        for (const auto& record : output.result.out_records) {
          io::VcfWriter::WriteBgzfRecord(bgzf.get(), hdr, record);
        }
      }

      if (shap_writers.combined && !output.result.shap_value_rows.empty()) {
        shap_writers.combined->AppendRows(output.result.shap_value_rows);
      }
      if (shap_writers.snv && !output.result.snv_shap_value_rows.empty()) {
        shap_writers.snv->AppendRows(output.result.snv_shap_value_rows);
      }
      if (shap_writers.indel && !output.result.indel_shap_value_rows.empty()) {
        shap_writers.indel->AppendRows(output.result.indel_shap_value_rows);
      }
      ++regions_written;
    }
  } catch (const std::runtime_error&) {
    exception_out = std::current_exception();
    queue.Shutdown();
  }
}
}  // namespace

void FilterVariantsClass::WriteVcfOutput(const fs::path& header_file, const vec<fs::path>& region_files) const {
  // Concatenate header + per-region BGZF data files at the block level (no record re-reading).
  CreateParentDirectoryIfNotExists(_param.vcf_output);
  ConcatenateBgzfFiles(_param.vcf_output, header_file, region_files);
  if (bcf_index_build(_param.vcf_output.c_str(), 0) != 0) {
    Logging::Warn("Failed to build index for {}", _param.vcf_output.string());
  }

  // A non-empty partition means the input contained variants in the requested regions and filtering ran.
  // Distinguish that case from the earlier zero-region path: an empty final VCF now means every variant
  // was filtered out, which is a valid result but should be visible to the caller.
  if (!_partitioned_regions.empty()) {
    const io::VcfReader output_reader(_param.vcf_output);
    if (output_reader.GetNextRecord() == nullptr) {
      WarnAsErrorIfSet("No variants found in VCF after applying filters. Writing empty outputs.");
    }
  }
}

void FilterVariantsClass::WriteMetrics() const {
  const auto& metrics_output = _param.metrics_output;
  // No metrics path set (e.g. an internal caller constructing the param struct directly): skip.
  if (metrics_output.empty()) {
    return;
  }

  VcfMetrics metrics;
  metrics.num_input_variants = CountInputVariants(_param.vcf_file);

  // Accumulate output-side metrics from the final output VCF (single pass).
  const io::VcfReader reader(_param.vcf_output);
  while (const auto record = reader.GetNextRecord()) {
    metrics.AccumulateOutputRecord(*record);
  }

  // Combine all sample names (VCF column order) into a single comma-separated ##Sample value.
  const auto sample_indexes = reader.GetHeader()->GetSampleIndexes();
  vec<std::string> sample_names(sample_indexes.size());
  for (const auto& [name, index] : sample_indexes) {
    sample_names[ToUnsigned(index)] = name;
  }
  const std::string sample_name = string::Join(sample_names, ",");

  const io::CommandLineInfo command_line = _param.command_line.value_or(io::CommandLineInfo{});
  // Predicted biological sex comes from the global context, populated by SetGlobalContext(). On the
  // empty-output path (no variants), metrics are written before the context is set, so the sex is
  // its default of Unknown.
  WriteVcfMetricsFile(metrics_output, metrics, command_line, sample_name, _global_ctx.sex, _param.sample_type);
  Logging::Info("Wrote QC metrics to {}", metrics_output.string());
}

void FilterVariantsClass::ParallelFiltering() {
  // Create a unique temp directory for per-region BGZF data files.
  // TmpDir uses mkdtemp for atomic, collision-free creation and installs signal handlers for cleanup.
  const TmpDir tmp_dir_guard(_param.vcf_output.parent_path(), ".svc_tmp");
  const auto& tmp_dir = tmp_dir_guard.Path();

  // Write the VCF header to a temp file. This triggers bcf_hdr_sync which rebuilds the
  // header's internal dictionary after AddFilterLine/AddFormatLine calls. Must happen
  // before SetGlobalContext() which reads the dictionary via GetFieldMetadata().
  const auto header_file = tmp_dir / "header.vcf.gz";
  {
    const io::VcfWriter header_writer(header_file, _hdr, io::NoIndex{});
    header_writer.WriteHeader();
  }

  // Graceful empty output: there is nothing to filter, either because the input VCF holds no
  // records or because no variant survived the target regions. Write a valid header-only VCF and
  // return successfully instead of running the filtering machinery over zero regions. Only a
  // corrupt or unreadable input is an error; that is rejected earlier by the input validators.
  if (_partitioned_regions.empty()) {
    WarnAsErrorIfSet("No variants to filter in VCF file '{}'. Writing empty outputs.", _param.vcf_file);
    WriteVcfOutput(header_file, {});
    WriteMetrics();
    return;
  }

  // Initialize the global context — shared read-only state across all regions and tasks.
  SetGlobalContext();
  if (auto num_warnings = CheckVcfFeatureResources(_model_config.vcf_feature_cols,
                                                   _global_ctx.interest_regions,
                                                   _hdr,
                                                   _param.pop_af_vcf.has_value(),
                                                   _param.catalog_vcf.has_value());
      num_warnings > 0) {
    WarnAsErrorIfSet(
        "There were {} warnings while checking VCF feature resources; VCF feature extraction may not be accurate!",
        num_warnings);
  }

  const auto num_regions = _partitioned_regions.size();
  Progress progress(num_regions);

  // Build the ordered list of per-region temp file paths.
  vec<fs::path> region_files;
  region_files.reserve(num_regions);
  for (const auto idx : std::views::iota(size_t{0}, num_regions)) {
    region_files.emplace_back(tmp_dir / fmt::format("region_{}.vcf.gz", idx));
  }

  RegionOutputQueue queue;

  // Set up per-thread worker contexts.
  _workers.reserve(_param.threads);
  for (u32 i = 0; i < _param.threads; ++i) {
    try {
      const AnnotationVcfPaths annotation_vcfs{
          .popaf = _param.pop_af_vcf, .pon = _param.pon_vcf, .catalog = _param.catalog_vcf};
      auto worker_ctx =
          std::make_unique<WorkerContext>(_param.vcf_file, annotation_vcfs, _param.bam_files, _alignment_reader_cache);
      worker_ctx->calculators = CreateScoreCalculators(_param, _model_config);
      _workers.emplace_back(_global_ctx, worker_ctx, _hdr);
    } catch (std::exception& e) {
      throw error::Error("Error creating WorkerContext: {}", e.what());
    }
  }

  // Build Taskflow graph: one filter task per region.
  tf::Taskflow taskflow;
  tf::Executor executor{_param.threads};
  const auto workflow = _model_config.workflow;

  for (size_t idx = 0; idx < num_regions; ++idx) {
    taskflow.emplace([&workers = _workers,
                      workflow,
                      &region = _partitioned_regions[idx],
                      region_idx = idx,
                      &queue,
                      &executor,
                      &progress]() {
      auto& worker = workers.at(static_cast<size_t>(executor.this_worker_id()));
      RunFilterTask(worker, workflow, region, region_idx, queue, progress);
    });
  }

  // Dedicated writer thread: pops completed regions and writes each to a headerless BGZF data file.
  const ShapWriters shap_writers{
      .combined = _shap_value_writer, .snv = _snv_shap_value_writer, .indel = _indel_shap_value_writer};
  std::exception_ptr writer_exception;
  std::jthread writer_thread([&num_regions, &queue, &region_files, &hdr = _hdr, &shap_writers, &writer_exception]() {
    WriteRegionResults(num_regions, queue, region_files, hdr, shap_writers, writer_exception);
  });

  Logging::Info("Filtering and writing '{}' regions...", progress.total_region_count);

  // Run all filter tasks. On failure, shut down the queue so the writer thread exits.
  std::exception_ptr filter_exception;
  try {
    executor.run(taskflow).get();
  } catch (const std::runtime_error&) {
    filter_exception = std::current_exception();
    queue.Shutdown();
  }

  writer_thread.join();

  if (!filter_exception && !writer_exception) {
    WriteVcfOutput(header_file, region_files);
    WriteMetrics();
  }

  // Propagate exceptions — filter exceptions take priority.
  // tmp_dir is cleaned up by TmpDir destructor on scope exit.
  if (filter_exception) {
    std::rethrow_exception(filter_exception);
  }
  if (writer_exception) {
    std::rethrow_exception(writer_exception);
  }

  if (num_regions != progress.total_region_count) {
    throw error::Error("Number of regions {} and tasks executed {} are not the same!",
                       num_regions,
                       std::to_string(progress.total_region_count));
  }
}

void FilterVariantsClass::SetReferenceSequences() {
  StrSet chrom_names{};
  for (const auto& r : _partitioned_regions) {
    // Only store the chromosomes that have variants in the VCF file that need to be processed
    chrom_names.insert(r.chrom);
  }
  io::FastaReader fasta_reader(_param.genome);
  for (const auto& [chrom, length] : _hdr->GetContigLengths()) {
    if (chrom_names.contains(chrom)) {
      _ref_seqs[chrom] = fasta_reader.GetSequence(chrom, 0, static_cast<s32>(length));
    }
  }
}

void FilterVariantsClass::SetShapValueTsvWriters() {
  // Set up SHAP values output writers
  if (_param.shap_value_tsv.has_value()) {
    _shap_value_writer = std::make_shared<LockedTsvWriter>(_param.shap_value_tsv.value());
    WriteShapValueTsvHeader(*_shap_value_writer, _model_config.scoring_names, _param.command_line);
  }
  if (_param.snv_shap_value_tsv.has_value()) {
    _snv_shap_value_writer = std::make_shared<LockedTsvWriter>(_param.snv_shap_value_tsv.value());
    WriteShapValueTsvHeader(*_snv_shap_value_writer, _model_config.snv_scoring_names, _param.command_line);
  }
  if (_param.indel_shap_value_tsv.has_value()) {
    _indel_shap_value_writer = std::make_shared<LockedTsvWriter>(_param.indel_shap_value_tsv.value());
    WriteShapValueTsvHeader(*_indel_shap_value_writer, _model_config.indel_scoring_names, _param.command_line);
  }
}

void FilterVariantsClass::AddDuplexFormatLines() const {
  using enum io::FieldType;
  if (IsDuplexProtocol(_param.sequencing_protocol)) {
    _hdr->AddFormatLine({.id = kDuplexConcordantCountsId,
                         .description = kDuplexConcordantCountsDesc,
                         .number = io::kNumberR,
                         .type = kInteger});
    _hdr->AddFormatLine({.id = kDuplexSimplexCountsId,
                         .description = kDuplexSimplexCountsDesc,
                         .number = io::kNumberR,
                         .type = kInteger});
    _hdr->AddFormatLine({.id = kDuplexDiscordantCountsId,
                         .description = kDuplexDiscordantCountsDesc,
                         .number = io::kNumberR,
                         .type = kInteger});
    _hdr->AddFormatLine(
        {.id = kDuplexLowbqCountsId, .description = kDuplexLowbqCountsDesc, .number = io::kNumberR, .type = kInteger});
  }
}

void FilterVariantsClass::FilterGermline() {
  // Verify that the model files and scoring columns are compatible. The specified scoring columns must match the
  // features used in the model both in name and order. Any mismatch means the models are not compatible with the
  // requested workflow based on the provided config and will not produce reliable results.
  VerifyModelCompatibility(_param.snv_model, _model_config.snv_scoring_cols);
  VerifyModelCompatibility(_param.indel_model, _model_config.indel_scoring_cols);

  Logging::Info("Loading VCF file {}", _param.vcf_file);
  const io::VcfReader reader(_param.vcf_file);
  if (1 != reader.GetHeader()->GetNumSamples()) {
    // Germline filtering is designed for a single sample only.
    throw error::Error("VCF file must contain exactly one sample");
  }

  CreateParentDirectoryIfNotExists(_param.vcf_output);

  Logging::Info("Writing output to VCF file {}", _param.vcf_output.string());
  // Add germline specific header information
  _hdr = reader.GetHeader()->Clone();
  ReplacePopafInfoLine(_hdr);
  if (_param.command_line) {
    _hdr->AddCustomMetaDataLine(io::GetCommandInfo(*_param.command_line));
  }
  _hdr->AddFilterLine({.id = kFilteringPassId, .description = kFilteringGermlinePassDesc});
  _hdr->AddFilterLine({.id = kFilteringFailId, .description = kFilteringGermlineFailDesc});
  _hdr->AddFilterLine({.id = kFilteringMissingFeatureId, .description = kFilteringMissingFeatureDesc});
  _hdr->AddFilterLine({.id = kFilteringFalsePositiveId, .description = kFilteringFalsePositiveDesc});
  _hdr->AddFilterLine({.id = kFilteringMultialleleFormatId, .description = kFilteringMultialleleFormatDesc});
  _hdr->AddFilterLine({.id = kFilteringMultiallelePartnerId, .description = kFilteringMultiallelePartnerDesc});
  _hdr->AddFilterLine({.id = kFilteringMultialleleConflictId, .description = kFilteringMultialleleConflictDesc});
  _hdr->AddFilterLine({.id = kFilteringNonAcgtRefAltId, .description = kFilteringNonAcgtRefAltDesc});
  _hdr->AddFormatLine(
      {.id = kGermlineMLId, .description = kGermlineMLDesc, .number = io::kNumberOne, .type = io::FieldType::kInteger});
  _hdr->AddFormatLine({.id = kMachineLearningId,
                       .description = kMachineLearningDesc,
                       .number = io::kNumberOne,
                       .type = io::FieldType::kFloat});
  AddDuplexFormatLines();
  _hdr->AddFormatLine({.id = kGermlineGnomadAFId,
                       .description = kGermlineGnomadAFDesc,
                       .number = io::kNumberEachAllele,
                       .type = io::FieldType::kFloat});
  _hdr->AddFormatLine({.id = kGermlineRefAvgMapqId,
                       .description = kGermlineRefAvgMapqDesc,
                       .number = io::kNumberEachAllele,
                       .type = io::FieldType::kFloat});
  _hdr->AddFormatLine({.id = kGermlineAltAvgMapqId,
                       .description = kGermlineAltAvgMapqDesc,
                       .number = io::kNumberEachAllele,
                       .type = io::FieldType::kFloat});
  _hdr->AddFormatLine({.id = kGermlineRefAvgDistId,
                       .description = kGermlineRefAvgDistDesc,
                       .number = io::kNumberEachAllele,
                       .type = io::FieldType::kFloat});
  _hdr->AddFormatLine({.id = kGermlineAltAvgDistId,
                       .description = kGermlineAltAvgDistDesc,
                       .number = io::kNumberEachAllele,
                       .type = io::FieldType::kFloat});
  _hdr->AddFormatLine({.id = kGermlineDensity100BPId,
                       .description = kGermlineDensity100BPDesc,
                       .number = io::kNumberOne,
                       .type = io::FieldType::kInteger});

  // Single parallel VCF pre-scan: collects DP values, variant positions, and variant density
  _pre_scan = ParallelVcfPreScan(_param.vcf_file, _param.threads, std::nullopt);

  // Partition using DP-weighted regions for balanced parallel filtering
  _partitioned_regions = PartitionVcfRegionsByDp(_pre_scan, _param.threads);

  // Load reference sequences based on which chromosomes have variants in the VCF that need to be processed
  SetReferenceSequences();

  SetShapValueTsvWriters();

  // Filter the input VCF in parallel and write the output
  ParallelFiltering();

  Logging::Info("Wrote output VCF file {}", _param.vcf_output);
}

void FilterVariantsClass::FilterGermlineTagging() {
  VerifyModelCompatibility(_param.model, _model_config.scoring_cols);

  CreateParentDirectoryIfNotExists(_param.vcf_output);

  Logging::Info("Loading VCF file {}", _param.vcf_file);
  const io::VcfReader reader(_param.vcf_file);
  Logging::Info("Writing output to VCF file {}", _param.vcf_output.string());
  // Add germline specific header information
  _hdr = reader.GetHeader()->Clone();
  ReplacePopafInfoLine(_hdr);
  if (_param.command_line) {
    _hdr->AddCustomMetaDataLine(io::GetCommandInfo(*_param.command_line));
  }
  _hdr->AddFormatLine({.id = kMachineLearningId,
                       .description = kMachineLearningDesc,
                       .number = io::kNumberOne,
                       .type = io::FieldType::kFloat});
  _hdr->AddInfoLine({.id = kGermlineTaggingInfoGermlineId,
                     .description = kGermlineTaggingInfoGermlineDesc,
                     .number = io::kNumberZero,
                     .type = io::FieldType::kFlag});
  _hdr->AddInfoLine({.id = kGermlineTaggingInfoSomaticId,
                     .description = kGermlineTaggingInfoSomaticDesc,
                     .number = io::kNumberZero,
                     .type = io::FieldType::kFlag});
  if (!_hdr->HasInfoField(kFilteringPassId)) {
    _hdr->AddFilterLine({.id = kFilteringPassId, .description = kFilteringPassDesc});
  }
  if (!_hdr->HasInfoField(kFilteringFailId)) {
    _hdr->AddFilterLine({.id = kFilteringFailId, .description = kFilteringFailDesc});
  }

  // Single parallel VCF pre-scan: collects DP values, variant positions, and variant density
  _pre_scan = ParallelVcfPreScan(_param.vcf_file, _param.threads, std::nullopt);

  // Partition using DP-weighted regions for balanced parallel filtering
  _partitioned_regions = PartitionVcfRegionsByDp(_pre_scan, _param.threads);

  // Load reference sequences based on which chromosomes have variants in the VCF that need to be processed
  SetReferenceSequences();

  SetShapValueTsvWriters();

  // Filter the input VCF in parallel and write the output
  ParallelFiltering();

  Logging::Info("Wrote output VCF file {}", _param.vcf_output);
}

void FilterVariantsClass::FilterTumorNormal() {
  using enum io::FieldType;

  CreateParentDirectoryIfNotExists(_param.vcf_output);

  Logging::Info("Loading VCF file {}", _param.vcf_file);
  const io::VcfReader reader(_param.vcf_file);
  Logging::Info("Writing output to VCF file {}", _param.vcf_output.string());
  // Add germline specific header information
  _hdr = reader.GetHeader()->Clone();

  const auto target_regions = GetChromIntervalMap(_param.bed_file);
  const bool detect_sample_type = RequiresDeferredSampleTypeDetection(_param);
  std::optional<SampleTypePreScanOptions> sample_type_pre_scan_options;
  if (detect_sample_type) {
    // Reuse the filter pre-scan VCF pass for sample-type metrics to avoid a separate detector pass over
    // the same indexed VCF.
    sample_type_pre_scan_options = SampleTypePreScanOptions{
        .sample_idx = DetermineTumorSampleIndex(_param.workflow, reader.GetHeader(), _param.vcf_file),
        .max_af = _param.sample_type_max_af,
        .total_length = ComputeSampleTypeDetectionLength(reader.GetHeader(), target_regions),
        .bed_regions = target_regions,
    };
  }

  // Pre-scan the full VCF so chromosome median DP and variant density include off-target chromosomes.
  // Sample-type metrics remain restricted by SampleTypePreScanOptions::bed_regions below.
  _pre_scan = ParallelVcfPreScan(_param.vcf_file,
                                 {
                                     .threads = _param.threads,
                                     .partition_bed_regions = std::nullopt,
                                     .sample_type = std::move(sample_type_pre_scan_options),
                                 });
  ResolveDeferredSampleType(_param, _pre_scan, detect_sample_type);

  VerifyModelCompatibility(_param.model, _model_config.scoring_cols);

  if (_param.command_line) {
    _hdr->AddCustomMetaDataLine(io::GetCommandInfo(*_param.command_line));
  }
  _hdr->AddFilterLine({.id = kFilteringPassId, .description = kFilteringGermlinePassDesc});
  _hdr->AddFilterLine({.id = kFilteringFailId, .description = kFilteringGermlineFailDesc});
  _hdr->AddFilterLine({.id = kFilteringFailMinTumorAfId, .description = kFilteringFailMinTumorAfDesc});
  _hdr->AddFilterLine({.id = kFilteringMLScoreId, .description = kFilteringMLScoreDesc});
  _hdr->AddFilterLine({.id = kFilteringFailMinTumorSupportId, .description = kFilteringFailMinTumorSupportDesc});
  _hdr->AddFilterLine({.id = kFilteringFailMaxNormalSupportId, .description = kFilteringFailMaxNormalSupportDesc});
  _hdr->AddFilterLine({.id = kFilteringFailMinDpRatioId, .description = kFilteringMinDpRatioDesc});
  _hdr->AddFilterLine({.id = kFilteringFailMaxIndelSizeId, .description = kFilteringFailMaxIndelSizeDesc});
  _hdr->AddFilterLine({.id = kFilteringMissingFeatureId, .description = kFilteringMissingFeatureDesc});
  _hdr->AddFilterLine({.id = kFilteringNonAcgtRefAltId, .description = kFilteringNonAcgtRefAltDesc});
  _hdr->AddInfoLine({.id = kRefBQId, .description = kRefBQDesc, .number = io::kNumberOne, .type = kFloat});
  _hdr->AddInfoLine({.id = kRefMQId, .description = kRefMQDesc, .number = io::kNumberOne, .type = kFloat});
  _hdr->AddInfoLine({.id = kAltBQId, .description = kAltBQDesc, .number = io::kNumberOne, .type = kFloat});
  _hdr->AddInfoLine({.id = kAltMQId, .description = kAltMQDesc, .number = io::kNumberOne, .type = kFloat});
  _hdr->AddInfoLine({.id = kSubtypeId, .description = kSubtypeDesc, .number = io::kNumberOne, .type = kFloat});
  _hdr->AddInfoLine({.id = kContextId, .description = kContextDesc, .number = io::kNumberOne, .type = kFloat});
  ReplacePopafInfoLine(_hdr);
  _hdr->AddFormatLine(
      {.id = kMachineLearningId, .description = kMachineLearningDesc, .number = io::kNumberOne, .type = kFloat});
  _hdr->AddFormatLine({.id = kBaseqQualId, .description = kBaseQualDesc, .number = io::kNumberOne, .type = kFloat});
  _hdr->AddFormatLine({.id = kMapQualId, .description = kMapQualDesc, .number = io::kNumberOne, .type = kFloat});
  _hdr->AddFormatLine({.id = kDistanceId, .description = kDistanceDesc, .number = io::kNumberOne, .type = kFloat});
  AddDuplexFormatLines();

  // Partition using DP-weighted regions for balanced parallel filtering, restricted to target regions.
  _partitioned_regions = PartitionVcfRegionsByDp(_pre_scan, _param.threads, target_regions);

  // Load reference sequences based on which chromosomes have variants in the VCF that need to be processed
  SetReferenceSequences();

  SetShapValueTsvWriters();

  // Filter the input VCF in parallel and write the output
  ParallelFiltering();

  Logging::Info("Wrote output VCF file {}", _param.vcf_output);
}

void FilterVariantsClass::FilterTumorOnlyWgs() {
  using enum io::FieldType;

  // Panel of Normals is required for tumor-only-wgs (enforced by the CLI): PON features are part of the scoring schema.
  const bool detect_sample_type = RequiresDeferredSampleTypeDetection(_param);
  if (!detect_sample_type) {
    VerifyModelCompatibility(_param.model, _model_config.scoring_cols);
  }

  CreateParentDirectoryIfNotExists(_param.vcf_output);

  Logging::Info("Loading VCF file {}", _param.vcf_file);
  const io::VcfReader reader(_param.vcf_file);
  if (1 != reader.GetHeader()->GetNumSamples()) {
    throw error::Error("VCF file must contain exactly one sample");
  }
  Logging::Info("Writing output to VCF file {}", _param.vcf_output.string());
  _hdr = reader.GetHeader()->Clone();

  const auto target_regions = GetChromIntervalMap(_param.bed_file);
  std::optional<SampleTypePreScanOptions> sample_type_pre_scan_options;
  if (detect_sample_type) {
    sample_type_pre_scan_options = SampleTypePreScanOptions{
        .sample_idx = DetermineTumorSampleIndex(_param.workflow, reader.GetHeader(), _param.vcf_file),
        .max_af = _param.sample_type_max_af,
        .total_length = ComputeSampleTypeDetectionLength(reader.GetHeader(), target_regions),
        .bed_regions = target_regions,
    };
  }

  // Reuse the filtering pre-scan to collect sample-type metrics without another pass over the VCF.
  _pre_scan = ParallelVcfPreScan(_param.vcf_file,
                                 {
                                     .threads = _param.threads,
                                     .partition_bed_regions = std::nullopt,
                                     .sample_type = std::move(sample_type_pre_scan_options),
                                 });
  ResolveDeferredSampleType(_param, _pre_scan, detect_sample_type);
  if (detect_sample_type) {
    VerifyModelCompatibility(_param.model, _model_config.scoring_cols);
  }

  if (_param.command_line) {
    _hdr->AddCustomMetaDataLine(io::GetCommandInfo(*_param.command_line));
  }
  _hdr->AddFilterLine({.id = kFilteringPassId, .description = kFilteringGermlinePassDesc});
  _hdr->AddFilterLine({.id = kFilteringFailId, .description = kFilteringGermlineFailDesc});
  _hdr->AddFilterLine({.id = kFilteringMLScoreId, .description = kFilteringMLScoreDesc});
  _hdr->AddFilterLine({.id = kFilteringFailMinSupportId, .description = kFilteringFailMinSupportDesc});
  _hdr->AddFilterLine({.id = kFilteringFailMinAfId, .description = kFilteringFailMinAfDesc});
  _hdr->AddFilterLine({.id = kFilteringFailMaxIndelSizeId, .description = kFilteringFailMaxIndelSizeDesc});
  _hdr->AddFilterLine({.id = kFilteringMissingFeatureId, .description = kFilteringMissingFeatureDesc});
  _hdr->AddFilterLine({.id = kFilteringNonAcgtRefAltId, .description = kFilteringNonAcgtRefAltDesc});
  _hdr->AddInfoLine({.id = kRefBQId, .description = kRefBQDesc, .number = io::kNumberOne, .type = kFloat});
  _hdr->AddInfoLine({.id = kRefMQId, .description = kRefMQDesc, .number = io::kNumberOne, .type = kFloat});
  _hdr->AddInfoLine({.id = kAltBQId, .description = kAltBQDesc, .number = io::kNumberOne, .type = kFloat});
  _hdr->AddInfoLine({.id = kAltMQId, .description = kAltMQDesc, .number = io::kNumberOne, .type = kFloat});
  _hdr->AddInfoLine({.id = kSubtypeId, .description = kSubtypeDesc, .number = io::kNumberOne, .type = kInteger});
  _hdr->AddInfoLine({.id = kContextId, .description = kContextDesc, .number = io::kNumberOne, .type = kString});
  ReplacePopafInfoLine(_hdr);
  _hdr->AddFormatLine(
      {.id = kMachineLearningId, .description = kMachineLearningDesc, .number = io::kNumberOne, .type = kFloat});
  _hdr->AddFormatLine({.id = kBaseqQualId, .description = kBaseQualDesc, .number = io::kNumberOne, .type = kFloat});
  _hdr->AddFormatLine({.id = kMapQualId, .description = kMapQualDesc, .number = io::kNumberOne, .type = kFloat});
  _hdr->AddFormatLine({.id = kDistanceId, .description = kDistanceDesc, .number = io::kNumberOne, .type = kFloat});
  AddDuplexFormatLines();

  // Partition using DP-weighted regions for balanced parallel filtering, restricted to target regions.
  _partitioned_regions = PartitionVcfRegionsByDp(_pre_scan, _param.threads, target_regions);

  SetReferenceSequences();

  SetShapValueTsvWriters();

  ParallelFiltering();

  Logging::Info("Wrote output VCF file {}", _param.vcf_output);
}

void FilterVariantsClass::FilterTumorOnlyTe() {
  VerifyModelCompatibility(_param.model, _model_config.scoring_cols);

  // load the VCF File
  Logging::Info("Loading VCF file {}", _param.vcf_file);
  const io::VcfReader reader(_param.vcf_file);
  if (1 != reader.GetHeader()->GetNumSamples()) {
    throw error::Error("VCF file must contain exactly one sample");
  }

  CreateParentDirectoryIfNotExists(_param.vcf_output);

  Logging::Info("Writing output to VCF file {}", _param.vcf_output.string());

  _hdr = reader.GetHeader()->Clone();
  if (_param.command_line) {
    _hdr->AddCustomMetaDataLine(io::GetCommandInfo(*_param.command_line));
  }
  // Append descriptions for the FILTER fields
  _hdr->AddFilterLine({.id = kFilteringPassId, .description = kFilteringPassDesc});
  _hdr->AddFilterLine({.id = kFilteringMapQualityId, .description = kFilteringMapQualityDesc});
  _hdr->AddFilterLine({.id = kFilteringBlocklistedId, .description = kFilteringBlocklistedDesc});
  _hdr->AddFilterLine({.id = kFilteringMLScoreId, .description = kFilteringMLScoreDesc});
  _hdr->AddFilterLine({.id = kFilteringCountsId, .description = kFilteringCountsDesc});
  _hdr->AddFilterLine({.id = kFilteringBaseQualityId, .description = kFilteringBaseQualityDesc});
  _hdr->AddFilterLine({.id = kFilteringForcedId, .description = kFilteringForcedDesc});
  _hdr->AddFilterLine({.id = kFilteringAFId, .description = kFilteringAFDesc});
  _hdr->AddFilterLine({.id = kFilteringMinAltCountsId, .description = kFilteringMinAltCountsDesc});

  // Append descriptions for the FORMAT fields
  _hdr->AddFormatLine({.id = kWeightedCountsId,
                       .description = kWeightedCountsDesc,
                       .number = io::kNumberEachAllele,
                       .type = io::FieldType::kFloat});
  _hdr->AddFormatLine({.id = kMachineLearningId,
                       .description = kMachineLearningDesc,
                       .number = io::kNumberOne,
                       .type = io::FieldType::kFloat});
  _hdr->AddFormatLine({.id = kNonDuplexCountsId,
                       .description = kNonDuplexCountsDesc,
                       .number = io::kNumberOne,
                       .type = io::FieldType::kInteger});
  _hdr->AddFormatLine({.id = kDuplexCountsId,
                       .description = kDuplexCountsDesc,
                       .number = io::kNumberOne,
                       .type = io::FieldType::kInteger});
  _hdr->AddFormatLine({.id = kStrandBiasMetricId,
                       .description = kStrandBiasMetricDesc,
                       .number = io::kNumberOne,
                       .type = io::FieldType::kFloat});
  _hdr->AddFormatLine({.id = kPlusOnlyCountsId,
                       .description = kPlusOnlyCountsDesc,
                       .number = io::kNumberOne,
                       .type = io::FieldType::kInteger});
  _hdr->AddFormatLine({.id = kMinusOnlyCountsId,
                       .description = kMinusOnlyCountsDesc,
                       .number = io::kNumberOne,
                       .type = io::FieldType::kInteger});
  _hdr->AddFormatLine({.id = kSequenceContextId,
                       .description = kSequenceContextDesc,
                       .number = io::kNumberOne,
                       .type = io::FieldType::kString});
  _hdr->AddFormatLine({.id = kPhysicalPhasedId,
                       .description = kPhysicalPhasedDesc,
                       .number = io::kNumberOne,
                       .type = io::FieldType::kString});
  _hdr->AddFormatLine(
      {.id = kHotspotId, .description = kHotspotDesc, .number = io::kNumberOne, .type = io::FieldType::kInteger});
  AddDuplexFormatLines();

  // Load BED regions if provided
  const auto bed_regions = GetChromIntervalMap(_param.bed_file);

  // Pre-scan the full VCF (without BED filtering) so that median DP and variant density are
  // computed from all variants. BED filtering is applied only during region partitioning.
  _pre_scan = ParallelVcfPreScan(_param.vcf_file, _param.threads, std::nullopt);

  // Partition using DP-weighted regions, restricted to BED intervals when provided
  _partitioned_regions = PartitionVcfRegionsByDp(_pre_scan, _param.threads, bed_regions);

  // Load reference sequences based on which chromosomes have variants in the VCF that need to be processed
  SetReferenceSequences();

  SetShapValueTsvWriters();

  // Filter E2E by region and write the output
  ParallelFiltering();

  Logging::Info("Wrote output VCF file {}", _param.vcf_output);
}

void FilterVariantsClass::CheckTumorNormalSampleName() const {
  if (!_param.tumor_sample_name.has_value() || _param.tumor_sample_name->empty()) {
    throw error::Error("Tumor sample name for read groups is not specified");
  }
  Logging::Info("Tumor sample name for read groups: {}", _param.tumor_sample_name.value());
}

void FilterVariantsClass::VerifyParameters() const {
  using enum Workflow;
  switch (_model_config.workflow) {
    case kGermlineTagging:
    case kTumorNormalWgs: {
      CheckTumorNormalSampleName();
      break;
    }
    case kTumorOnlyTe:
    case kTumorOnlyWgs:
    case kGermlineWgs:
    case kGermlineMitochondria:
      break;
    default:
      throw error::Error("Unsupported workflow");
  }

  if (_param.bam_files.empty()) {
    throw error::Error("No BAM file provided");
  }
}

void FilterVariantsClass::FilterMitochondria() {
  Logging::Info("Starting {} workflow", enum_util::FormatEnumName(Workflow::kGermlineMitochondria));

  CreateParentDirectoryIfNotExists(_param.vcf_output);

  Logging::Info("Loading VCF file {}", _param.vcf_file);
  const io::VcfReader reader(_param.vcf_file);
  if (reader.GetHeader()->GetNumSamples() != 1) {
    throw error::Error("VCF file must contain exactly one sample");
  }
  Logging::Info("Writing output to VCF file {}", _param.vcf_output.string());

  _hdr = reader.GetHeader()->Clone();
  // The germline-mitochondria output emits only the fields SVC produces. Drop every original
  // (input-inherited) INFO, FORMAT, and FILTER definition before re-adding the SVC FILTER/FORMAT lines
  // below, so the output header carries no upstream INFO/FORMAT/FILTER lines (contigs, reference, and
  // other meta are retained). Output records are built fresh in FilterMitochondriaRegion, so no
  // inherited per-record INFO/FORMAT values carry through either.
  _hdr->RemoveAllInfoLines();
  _hdr->RemoveAllFormatLines();
  _hdr->RemoveAllFilterLines();
  if (_param.command_line) {
    _hdr->AddCustomMetaDataLine(io::GetCommandInfo(*_param.command_line));
  }
  // Append descriptions for the FILTER fields applied by AF/depth filtering (PASS + mito low-AF/low-DP).
  // These are the only FILTER lines present after stripping the input definitions above.
  AddMitochondriaFilterHeaderLines(_hdr, _param.min_af, _param.min_dp);
  if (_param.blocklist_bed.has_value()) {
    AddMitochondriaBlocklistFilterHeaderLine(_hdr);
  }
  // The AF/depth filtering writes AD/DP/AF (and duplex count) FORMAT fields, matching the germline workflow's
  // data fields. These are the only FORMAT lines present after stripping the input definitions above.
  AddMitochondriaDataFormatLines(_hdr, IsDuplexProtocol(_param.sequencing_protocol));
  // Rebuild the header's internal dictionary after the remove/add edits above, before any record is
  // created against it (CreateFromHeader in FilterMitochondriaRegion) or its metadata is read. The
  // header is also synced later by WriteHeader(), but the removals above re-index the dictionary, so
  // sync explicitly here rather than relying on that downstream write.
  _hdr->Sync();

  // Single parallel VCF pre-scan: collects DP values, variant positions, and variant density
  _pre_scan = ParallelVcfPreScan(_param.vcf_file, _param.threads, std::nullopt);

  // An input VCF with no variants is handled gracefully by `ParallelFiltering`, which writes a
  // header-only output VCF.
  if (_pre_scan.chrom_variants.size() > 1) {
    Logging::Warn("Multiple chromosomes found in VCF file. Chromosomes found: {}",
                  fmt::join(std::views::keys(_pre_scan.chrom_variants), ", "));
  }

  // Partition using DP-weighted regions for balanced parallel filtering
  _partitioned_regions = PartitionVcfRegionsByDp(_pre_scan, _param.threads);

  // Load reference sequences based on which chromosomes have variants in the VCF that need to be processed
  SetReferenceSequences();

  SetShapValueTsvWriters();

  // Filter the input VCF in parallel and write the output
  ParallelFiltering();

  Logging::Info("Wrote output VCF file {}", _param.vcf_output);
}

/**
 * The main entry point for the filter_variants module. This function performs checks to ensure the correct number of
 * models are passed based on the specified workflow before calling a workflow specific filtering function.
 * @param param A FilterVariantsParam struct that contains required input/output filepaths and filtering settings
 */
void FilterVariants(const FilterVariantsParam& param) {
  using enum Workflow;
  // This function acts as the entry point for the filter-variants tool. Based on the workflow the respective filtering
  // process is run

  // Create Filtering object
  FilterVariantsClass filtering(param);

  // Verify that the correct number of models have been passed based on the workflow, and that BAM input has been passed
  filtering.VerifyParameters();

  switch (param.workflow) {
    case kGermlineWgs: {
      filtering.FilterGermline();
      break;
    }
    case kGermlineTagging: {
      filtering.FilterGermlineTagging();
      break;
    }
    case kTumorOnlyTe: {
      filtering.FilterTumorOnlyTe();
      break;
    }
    case kTumorNormalWgs: {
      filtering.FilterTumorNormal();
      break;
    }
    case kTumorOnlyWgs: {
      filtering.FilterTumorOnlyWgs();
      break;
    }
    case kGermlineMitochondria: {
      filtering.FilterMitochondria();
      break;
    }
    default:
      throw error::Error("Specified Workflow has not yet been implemented");
  }
  Logging::Info("Filtered all variants");
}

}  // namespace xoos::svc
