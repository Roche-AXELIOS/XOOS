#include "cli/copy-number-caller-cli.h"

#include <algorithm>
#include <array>
#include <optional>
#include <string>

#include <fmt/format.h>

#include <xoos/enum/enum-util.h>
#include <xoos/io/library-type.h>
#include <xoos/log/logging.h>
#include <xoos/util/file-functions.h>

#include "cli/copy-number-caller-cli-option-names.h"
#include "cli/copy-number-caller-subcommands.h"
#include "cli/define-copy-number-caller-options.h"
#include "coverage/calculate-coverage.h"
#include "gc-correct/gc-correct.h"
#include "germline-normal-wgs/germline-normal-wgs.h"
#include "intervals/augment-intervals.h"
#include "likelihood/likelihood.h"
#include "merge-segments/merge-segments-main.h"
#include "one-sample-logr/one-sample-logr.h"
#include "predict-sex/predict-sex.h"
#include "purity-ploidy-search/purity-ploidy-search.h"
#include "seg-to-vcf/seg-to-vcf.h"
#include "segmentation/pscbs.h"
#include "somatic-tumor-normal-wgs/somatic-tumor-normal-wgs.h"
#include "tumor-only/tumor-only-bafs.h"
#include "tumor-only/tumor-only-wgs.h"
#include "two-sample-logr/two-sample-logr.h"

namespace xoos::cnc {

namespace opt = cli_opt_name;

constexpr std::array<const char*, 12> kHiddenSubcommands = {
    kAugmentIntervalsSubcommand,
    kCalculateCoverageSubcommand,
    kGCCorrectSubcommand,
    kTwoSampleLogRSubcommand,
    kOneSampleLogRSubcommand,
    kSegmentationSubcommand,
    kPurityPloidySearchSubcommand,
    kPredictGermlineCNVSubcommand,
    kMergeSegmentsSubcommand,
    kSegToVcfSubcommand,
    kPredictSexSubcommand,
    kTumorOnlyBafsSubcommand,
};

namespace {
constexpr auto* kDefaultHg38GermlineBlocklistBed = "/resources/hg38_decoy_intervals.bed";
constexpr auto* kDefaultHg38TumorNormalBlocklistBed = "/resources/hg38_blocklist.bed";
constexpr auto* kDefaultHg38TumorOnlyBlocklistBed = "/resources/hg38_tumor_only_wgs_blocklist.bed";
constexpr auto* kDefaultHg38SnpBlocklistBed = "/resources/hg38_snp_blocklist.bed";

bool IsHiddenSubcommand(const std::string& name) {
  return std::ranges::find(kHiddenSubcommands, name) != std::ranges::end(kHiddenSubcommands);
}

void ApplyHiddenGroup(CLI::App* const subcommand) {
  if (subcommand != nullptr && IsHiddenSubcommand(subcommand->get_name())) {
    subcommand->group("");
  }
}

void ConfigureMappabilityBigWig(const CLI::App* const app, CopyNumberCallerOptions& options, const bool is_sbx_s) {
  if (app->count(opt::kMappabilityBigwig) > 0) {
    Logging::Info("Using user-specified mappability BigWig '{}'.", options.mappability_bigwig_fname.string());
    return;
  }
  options.mappability_bigwig_fname = is_sbx_s ? kHg38SbxSMappabilityBigWig : kHg38MappabilityBigWig;
  Logging::Info("Using default mappability BigWig '{}'.", options.mappability_bigwig_fname.string());
}

fs::path GetResourcePath(const fs::path& resource) {
  return resource;
}

fs::path GetResourcePath(const std::optional<fs::path>& resource) {
  return *resource;
}

template <typename Resource>
void ConfigureResource(const CLI::App* const app,
                       const char* const option_name,
                       Resource& resource,
                       const fs::path& bundled_resource) {
  const bool uses_bundled_resource = app->count(option_name) == 0;
  if (uses_bundled_resource) {
    resource = bundled_resource;
  }
  const auto resource_path = GetResourcePath(resource);
  if (uses_bundled_resource) {
    Logging::Info("Using bundled resource '{}'.", resource_path.string());
  } else {
    Logging::Info("Using user-specified resource '{}'.", resource_path.string());
  }
}

void AppendResolvedPathOptionToCommandLine(const CLI::App* const app,
                                           CopyNumberCallerOptions& options,
                                           const char* const option_name,
                                           const fs::path& resource) {
  const auto* const option = app->get_option_no_throw(option_name);
  if (option != nullptr && option->empty() && !resource.empty()) {
    options.command_line_info.command_line += fmt::format(" {} {}", option_name, resource.string());
  }
}

void AppendResolvedPathOptionToCommandLine(const CLI::App* const app,
                                           CopyNumberCallerOptions& options,
                                           const char* const option_name,
                                           const std::optional<fs::path>& resource) {
  if (resource.has_value()) {
    AppendResolvedPathOptionToCommandLine(app, options, option_name, *resource);
  }
}

void SetEffectiveIntervalSizeDefault(const CLI::App* const app, const CopyNumberCallerOptions& options) {
  const auto* const interval_size_option = app->get_option_no_throw(opt::kIntervalSize);
  if (interval_size_option == nullptr || !interval_size_option->empty()) {
    return;
  }

  // CLI11 exposes the callback app as const, but the app itself is mutable here.
  // This only updates CLI11's metadata default for rendering.
  const_cast<CLI::Option*>(interval_size_option)  // NOLINT(cppcoreguidelines-pro-type-const-cast) NOSONAR
      ->default_str(std::to_string(options.augment_intervals_options.whole_genome_interval_size));
}

}  // namespace

bool IsSbxSLibraryType(std::optional<io::LibraryType>& library_type, const fs::path& bam_path) {
  if (library_type.has_value()) {
    Logging::Info("Using user-specified library type '{}'.", enum_util::FormatEnumName(*library_type));
  } else {
    library_type = io::DetectLibraryType(bam_path);
    if (!library_type.has_value()) {
      Logging::Warn(
          "Could not determine library type from '{}'. Library-dependent defaults will use the non-SBX-S "
          "fallback values.",
          bam_path.string());
    }
  }
  return library_type == io::LibraryType::kSbxS;
}

void ConfigureSomaticLikelihoodOptions(CopyNumberCallerOptions& options) {
  options.likelihood_options.mode = LikelihoodMode::kSomatic;
  options.likelihood_options.mapq_cutoff_for_calls = kLikelihoodDefaultSomaticMAPQCutoffForCalls;
  using enum SegmentType;
  options.likelihood_options.input_segments_type = kBaf;
  if (options.bafs_fname.has_value()) {
    options.likelihood_options.output_segments_type =
        options.likelihood_options.tumor_only ? kSomaticTumorOnlyWithBafLikelihood : kSomaticWithBafLikelihood;
  } else {
    options.likelihood_options.output_segments_type = kSomaticNoBafLikelihood;
  }
}

void ConfigureGermlineLikelihoodOptions(CopyNumberCallerOptions& options) {
  options.sample_metadata_options.purity = kGermlinePurity;
  options.sample_metadata_options.ploidy = kDiploidPloidy;
  options.likelihood_options.mode = LikelihoodMode::kGermline;
  options.likelihood_options.input_segments_type = SegmentType::kLogROnly;
  options.likelihood_options.output_segments_type = SegmentType::kGermlineLikelihood;
}

void ConfigureAugmentIntervalsOptions(
    const CLI::App* const app,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    CopyNumberCallerOptions& options) {  // NOSONAR cpp:S995 cpp:M23_131
                                         // for later population during parse; opts must remain mutable
  if (app->count(opt::kIntervalSize) == 0) {
    options.augment_intervals_options.whole_genome_interval_size = kAugmentIntervalsDefaultWholeGenomeIntervalSize;
  }
  if (app->count(opt::kMinIntervalMappability) == 0) {
    options.augment_intervals_options.whole_genome_min_mappability = kAugmentIntervalsDefaultWholeGenomeMinMappability;
  }
}

void ConfigureSomaticTumorNormalWGSOptions(
    const CLI::App* const app,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    CopyNumberCallerOptions& options) {  // NOSONAR cpp:S995 cpp:M23_131
                                         // for later population during parse; opts must remain mutable
  options.augment_intervals_options.generate_whole_genome_intervals = true;
  options.mode = CopyNumberCallerModes::kSomaticTumorNormalWGS;
  options.segmentation_options.segmentation_mode = SegmentationMode::kSomatic;
  // An explicit mappability override does not require library detection, so avoid reading the BAM in that case.
  const bool is_sbx_s =
      app->count(opt::kMappabilityBigwig) == 0 && IsSbxSLibraryType(options.library_type, options.normal_bam_fname);
  ConfigureMappabilityBigWig(app, options, is_sbx_s);
  ConfigureResource(app, opt::kSeedSegments, options.seed_segments_fname, kHg38SeedSegments);
  // Restore somatic-specific defaults for bindings shared with the germline subcommand.
  // CLI11 writes default_val to the binding at registration time (last writer wins), so
  // whichever subcommand was registered last clobbers the other's defaults.
  if (app->count(opt::kSegmentUndoing) == 0) {
    options.segmentation_options.undo_logr_mode = SegmentUndoMode::kIterative;
  }
  if (app->count(opt::kIntervalSize) == 0) {
    options.augment_intervals_options.whole_genome_interval_size = kAugmentIntervalsDefaultWholeGenomeIntervalSize;
  }
  if (app->count(opt::kMinIntervalMappability) == 0) {
    options.augment_intervals_options.whole_genome_min_mappability = kAugmentIntervalsDefaultWholeGenomeMinMappability;
  }
  if (app->count(opt::kMinVcfNormalDepth) == 0) {
    options.baf_filter_options.normal_sample_min_depth = kVcfParsingOptionsDefaultSomaticNormalSampleMinDepth;
  }
  if (app->count(opt::kBafLikelihoodWeightForPurityPloidySearchRerun) == 0) {
    options.purity_ploidy_search_options.baf_likelihood_weight = kBafLikelihoodWeight;
  }
  ConfigureResource(app, opt::kBlocklistBed, options.blocklist_bed_fname, kDefaultHg38TumorNormalBlocklistBed);
}

void ConfigureSomaticTumorOnlyWGSOptions(
    const CLI::App* const app,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    CopyNumberCallerOptions& options) {  // NOSONAR cpp:S995 cpp:M23_131
                                         // for later population during parse; opts must remain mutable
  options.augment_intervals_options.generate_whole_genome_intervals = true;
  options.mode = CopyNumberCallerModes::kSomaticTumorOnlyWGS;
  options.segmentation_options.segmentation_mode = SegmentationMode::kSomatic;
  options.likelihood_options.mode = LikelihoodMode::kSomatic;
  // Tumor-only mode: purity/ploidy search and likelihood run without a matched normal. This enables
  // the tumor-only priors in PurityPloidySearchWithWGDQC and triggers ClassifyGermlineVsSomaticSegments
  // and SegmentType::kSomaticTumorOnlyWithBafLikelihood in CalculateLikelihoods.
  options.purity_ploidy_search_options.tumor_only = true;
  options.likelihood_options.tumor_only = true;
  // An explicit mappability override does not require library detection, so avoid reading the BAM in that case.
  const bool is_sbx_s =
      app->count(opt::kMappabilityBigwig) == 0 && IsSbxSLibraryType(options.library_type, options.tumor_bam_fname);
  ConfigureMappabilityBigWig(app, options, is_sbx_s);
  ConfigureResource(app, opt::kSeedSegments, options.seed_segments_fname, kHg38SeedSegments);
  // Restore somatic-specific defaults for bindings shared with the germline subcommand.
  // CLI11 writes default_val to the binding at registration time (last writer wins), so
  // whichever subcommand was registered last clobbers the other's defaults.
  if (app->count(opt::kSegmentUndoing) == 0) {
    options.segmentation_options.undo_logr_mode = SegmentUndoMode::kIterative;
  }
  if (app->count(opt::kIntervalSize) == 0) {
    options.augment_intervals_options.whole_genome_interval_size = kAugmentIntervalsDefaultWholeGenomeIntervalSize;
  }
  if (app->count(opt::kMinIntervalMappability) == 0) {
    options.augment_intervals_options.whole_genome_min_mappability = kAugmentIntervalsDefaultWholeGenomeMinMappability;
  }
  if (app->count(opt::kBafLikelihoodWeightForPurityPloidySearchRerun) == 0) {
    options.purity_ploidy_search_options.baf_likelihood_weight = kTumorOnlyWgsBafLikelihoodWeight;
  }
  ConfigureResource(app, opt::kBlocklistBed, options.blocklist_bed_fname, kDefaultHg38TumorOnlyBlocklistBed);
  ConfigureResource(app, opt::kSnpBlocklistBed, options.snp_blocklist_bed_fname, kDefaultHg38SnpBlocklistBed);
}

void ConfigurePurityPloidySearchOptions(const CLI::App* const app, CopyNumberCallerOptions& options) {
  if (app->count(opt::kBafLikelihoodWeightForPurityPloidySearchRerun) == 0) {
    options.purity_ploidy_search_options.baf_likelihood_weight = kBafLikelihoodWeight;
  }
}

void ConfigureGermlineNormalWGSOptions(
    const CLI::App* const app,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    CopyNumberCallerOptions& options) {  // NOSONAR cpp:S995 cpp:M23_131
                                         // for later population during parse; opts must remain mutable
  options.augment_intervals_options.generate_whole_genome_intervals = true;
  options.mode = CopyNumberCallerModes::kGermlineNormalWGS;
  options.sample_metadata_options.ploidy = kDiploidPloidy;
  options.sample_metadata_options.purity = kGermlinePurity;
  options.likelihood_options.mode = LikelihoodMode::kGermline;
  options.segmentation_options.cbs_method = CbsMaxTMethod::kFast;
  options.segmentation_options.truncate_outliers = true;
  options.segmentation_options.segmentation_mode = SegmentationMode::kGermline;
  options.baf_filter_options.normal_sample_min_baf = kVcfParsingOptionsDefaultGermlineNormalSampleMinBAF;
  options.baf_filter_options.normal_sample_max_baf = kVcfParsingOptionsDefaultGermlineNormalSampleMaxBAF;
  // Resolve the library type independently of individual library-specific defaults. An explicit interval size must not
  // prevent SBX-S filter defaults from being selected.
  const bool is_sbx_s = IsSbxSLibraryType(options.library_type, options.normal_bam_fname);
  ConfigureMappabilityBigWig(app, options, is_sbx_s);
  ConfigureResource(app, opt::kSeedSegments, options.seed_segments_fname, kHg38SeedSegments);
  options.likelihood_options.non_integer_tcn_filter_options =
      is_sbx_s ? kLikelihoodSbxSNonIntegerTcnFilterOptions : kLikelihoodDefaultNonIntegerTcnFilterOptions;
  if (is_sbx_s) {
    const auto& filter_options = options.likelihood_options.non_integer_tcn_filter_options;
    if (filter_options.dup_min_expected_tcn.has_value()) {
      Logging::Info(
          "SBX-S (simplex) library; using [{:.1f}, {:.1f}] as the expected TCN range for heterozygous deletions and "
          "{:.2f} as the minimum expected TCN for duplications.",
          filter_options.het_del_min_expected_tcn,
          filter_options.het_del_max_expected_tcn,
          filter_options.dup_min_expected_tcn.value());
    }
  }
  // Restore germline-specific defaults for bindings shared with the somatic subcommand.
  // CLI11 writes default_val to the binding at registration time (last writer wins), so
  // whichever subcommand was registered last clobbers the other's defaults.
  if (app->count(opt::kIntervalSize) == 0) {
    const size_t interval_size =
        is_sbx_s ? kAugmentIntervalsGermlineWGSSbxSIntervalSize : kAugmentIntervalsGermlineWGSDefaultIntervalSize;
    if (is_sbx_s) {
      Logging::Info("SBX-S (simplex) library; using {}-bp whole-genome interval size for germline WGS.", interval_size);
    }
    options.augment_intervals_options.whole_genome_interval_size = interval_size;
  }
  if (app->count(opt::kMinIntervalMappability) == 0) {
    options.augment_intervals_options.whole_genome_min_mappability =
        kAugmentIntervalsGermlineWGSDefaultWholeGenomeMinMappability;
  }
  if (app->count(opt::kMinVcfNormalDepth) == 0) {
    options.baf_filter_options.normal_sample_min_depth = kVcfParsingOptionsDefaultGermlineNormalSampleMinDepth;
  }
  ConfigureResource(app, opt::kBlocklistBed, options.blocklist_bed_fname, kDefaultHg38GermlineBlocklistBed);
  options.likelihood_options.mapq_cutoff_for_calls_is_user_set = app->count(opt::kMinMapqForCalls) > 0;
}

void ConfigureSegmentationOptions(
    const CLI::App* const app,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    CopyNumberCallerOptions& options) {  // NOSONAR cpp:S995 cpp:M23_131
                                         // for later population during parse; opts must remain mutable
  // Restore the segmentation subcommand's iterative default for --segment-undoing.
  // undo_logr_mode is a shared binding: CLI11 writes each subcommand's default_val to it at
  // registration time (last writer wins), and germline-wgs (registered last) writes kDisabled.
  // Re-assert kIterative here unless the user set the flag explicitly.
  if (app->count(opt::kSegmentUndoing) == 0) {
    options.segmentation_options.undo_logr_mode = SegmentUndoMode::kIterative;
  }
}

void CopyNumberCallerCliPreCallback(
    const CLI::App* const app,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    CopyNumberCallerOptions& options) {  // NOSONAR cpp:S995 cpp:M23_131
                                         // for later population during parse; opts must remain mutable
  file::CreateWritableDirectory(fs::absolute(options.output_dir));
  // update sample ID when sample ID is not exposed as an option to the user
  const auto* const sample_id_opt = app->get_option_no_throw(opt::kSampleId);
  if (sample_id_opt == nullptr) {
    if (options.sample_metadata_options.normal_sample_name.has_value() &&
        options.sample_metadata_options.tumor_sample_name.has_value()) {
      // in the case of tumor-normal paired analysis, use the format normal_tumor as the sample ID
      const auto& normal_name = options.sample_metadata_options.normal_sample_name.value();
      const auto& tumor_name = options.sample_metadata_options.tumor_sample_name.value();
      if (!normal_name.empty() && !tumor_name.empty()) {
        options.sample_metadata_options.sample_id = fmt::format("{}_{}", normal_name, tumor_name);
      }
    } else if (options.sample_metadata_options.normal_sample_name.has_value()) {
      // in the case of normal-only analysis, use the normal sample name as the sample ID
      const auto& normal_name = options.sample_metadata_options.normal_sample_name.value();
      if (!normal_name.empty()) {
        options.sample_metadata_options.sample_id = normal_name;
      }
    } else if (options.sample_metadata_options.tumor_sample_name.has_value()) {
      // in the case of tumor-only analysis, use the tumor sample name as the sample ID
      const auto& tumor_name = options.sample_metadata_options.tumor_sample_name.value();
      if (!tumor_name.empty()) {
        options.sample_metadata_options.sample_id = tumor_name;
      }
    } else {
      options.sample_metadata_options.sample_id = kSampleMetadataOptionsDefaultSampleId;
    }
  }
  const auto* const reference_opt = app->get_option_no_throw(opt::kReference);
  if (reference_opt != nullptr && reference_opt->count() > 0) {
    options.reference_genome_fai_fname = options.reference_genome_fname;
    options.reference_genome_fai_fname += ".fai";
  }
  SetEffectiveIntervalSizeDefault(app, options);
  options.command_line_info = cli::GetCommandLineInfo(app);
  AppendResolvedPathOptionToCommandLine(app, options, opt::kMappabilityBigwig, options.mappability_bigwig_fname);
  AppendResolvedPathOptionToCommandLine(app, options, opt::kSeedSegments, options.seed_segments_fname);
  AppendResolvedPathOptionToCommandLine(app, options, opt::kBlocklistBed, options.blocklist_bed_fname);
  AppendResolvedPathOptionToCommandLine(app, options, opt::kSnpBlocklistBed, options.snp_blocklist_bed_fname);
}

/// Register a hidden subcommand with the standard pre-parse callback.
namespace {
template <typename DefineOptsFn, typename MainFn>
void AddHiddenSubcommand(CLI::App* const app,
                         const char* name,
                         DefineOptsFn define_opts,
                         CopyNumberCallerOptionsPtr& options,
                         MainFn main_fn,
                         const char* description) {
  auto* sub = cli::AddSubcommand<CopyNumberCallerOptions>(
      app,
      name,
      define_opts,
      options,
      main_fn,
      description,
      [](const CLI::App* const app, const CopyNumberCallerOptionsPtr& opts) {
        CopyNumberCallerCliPreCallback(app, *opts);
      });
  ApplyHiddenGroup(sub);
}

void AddHiddenSubcommands(CLI::App* const app, CopyNumberCallerOptionsPtr& options) {
  // augment-intervals has a custom pre-parse callback to restore its defaults for options whose
  // bindings are shared with (and clobbered at registration time by) the WGS subcommands.
  auto* const augment_intervals_sub_app = cli::AddSubcommand<CopyNumberCallerOptions>(
      app,
      kAugmentIntervalsSubcommand,
      DefineAugmentIntervalsSubcommandOptions,
      options,
      AugmentIntervalsMain,
      "Generate or augment intervals.",
      [](const CLI::App* const app, const CopyNumberCallerOptionsPtr& opts) {
        ConfigureAugmentIntervalsOptions(app, *opts);
        CopyNumberCallerCliPreCallback(app, *opts);
      });
  ApplyHiddenGroup(augment_intervals_sub_app);
  AddHiddenSubcommand(app,
                      kCalculateCoverageSubcommand,
                      DefineCalculateCoverageSubcommandOptions,
                      options,
                      CalculateCoverageMain,
                      "Calculate coverage and mapping quality metrics.");
  AddHiddenSubcommand(app,
                      kGCCorrectSubcommand,
                      DefineGCCorrectSubcommandOptions,
                      options,
                      GCCorrectMain,
                      "Correct coverage for GC bias.");
  AddHiddenSubcommand(app,
                      kTwoSampleLogRSubcommand,
                      DefineTwoSampleLogRSubcommandOptions,
                      options,
                      TwoSampleLogRMain,
                      "Compute tumor/normal log ratios.");
  AddHiddenSubcommand(app,
                      kOneSampleLogRSubcommand,
                      DefineOneSampleLogRSubcommandOptions,
                      options,
                      OneSampleLogRMain,
                      "Compute log ratios from a single sample.");
  // segmentation has a custom pre-parse callback to restore its iterative --segment-undoing default,
  // whose binding is shared with (and clobbered at registration time by) the WGS subcommands.
  auto* const segmentation_sub_app = cli::AddSubcommand<CopyNumberCallerOptions>(
      app,
      kSegmentationSubcommand,
      DefineSegmentationSubcommandOptions,
      options,
      ParentSpecificBinarySegmentationMain,
      "Segment log ratio data.",
      [](const CLI::App* const app, const CopyNumberCallerOptionsPtr& opts) {
        ConfigureSegmentationOptions(app, *opts);
        CopyNumberCallerCliPreCallback(app, *opts);
      });
  ApplyHiddenGroup(segmentation_sub_app);
  auto* const purity_ploidy_search_sub_app = cli::AddSubcommand<CopyNumberCallerOptions>(
      app,
      kPurityPloidySearchSubcommand,
      DefinePurityPloidySearchSubcommandOptions,
      options,
      PurityPloidySearchMain,
      "Search for purity/ploidy solutions.",
      [](const CLI::App* const app, const CopyNumberCallerOptionsPtr& opts) {
        ConfigurePurityPloidySearchOptions(app, *opts);
        CopyNumberCallerCliPreCallback(app, *opts);
      });
  ApplyHiddenGroup(purity_ploidy_search_sub_app);

  // predict-germline has a custom pre-parse callback
  auto* const predict_germline_sub_app = cli::AddSubcommand<CopyNumberCallerOptions>(
      app,
      kPredictGermlineCNVSubcommand,
      DefinePredictGermlineCNVSubcommandOptions,
      options,
      CalculateLikelihoodsMain,
      "Predict germline CNVs from segments.",
      [](const CLI::App* const app, const CopyNumberCallerOptionsPtr& opts) {
        ConfigureGermlineLikelihoodOptions(*opts);
        CopyNumberCallerCliPreCallback(app, *opts);
      });
  ApplyHiddenGroup(predict_germline_sub_app);

  AddHiddenSubcommand(app,
                      kMergeSegmentsSubcommand,
                      DefineMergeSegmentsSubcommandOptions,
                      options,
                      MergeSegmentsMain,
                      "Merge adjacent segments.");
  AddHiddenSubcommand(app,
                      kSegToVcfSubcommand,
                      DefineSegToVcfSubcommandOptions,
                      options,
                      WriteSegmentsToVcfMain,
                      "Convert SEG files to VCF.");
  AddHiddenSubcommand(app,
                      kPredictSexSubcommand,
                      DefinePredictSexSubcommandOptions,
                      options,
                      PredictSexMain,
                      "Predict sex from coverage data.");
  AddHiddenSubcommand(app,
                      kTumorOnlyBafsSubcommand,
                      DefineTumorOnlyBafsSubcommandOptions,
                      options,
                      TumorOnlyBafsMain,
                      "Extract BAF values from a tumor-only VCF using a population SNP VCF.");
}
}  // namespace

void AddSubcommands(CLI::App* const app, CopyNumberCallerOptionsPtr& options) {
  AddHiddenSubcommands(app, options);

  cli::AddSubcommand<CopyNumberCallerOptions>(app,
                                              kPredictSomaticCNASubcommand,
                                              DefinePredictSomaticCNASubcommandOptions,
                                              options,
                                              CalculateLikelihoodsMain,
                                              "Predict somatic CNAs from segments.",
                                              [](const CLI::App* const app, const CopyNumberCallerOptionsPtr& opts) {
                                                ConfigureSomaticLikelihoodOptions(*opts);
                                                CopyNumberCallerCliPreCallback(app, *opts);
                                              });

  cli::AddSubcommand<CopyNumberCallerOptions>(app,
                                              kSomaticTumorNormalWGSSubcommand,
                                              DefineSomaticTumorNormalWGSOptions,
                                              options,
                                              SomaticTumorNormalWGSMain,
                                              "Run somatic tumor/normal WGS pipeline.",
                                              [](const CLI::App* const app, const CopyNumberCallerOptionsPtr& opts) {
                                                ConfigureSomaticTumorNormalWGSOptions(app, *opts);
                                                CopyNumberCallerCliPreCallback(app, *opts);
                                              });

  cli::AddSubcommand<CopyNumberCallerOptions>(app,
                                              kSomaticTumorOnlyWGSSubcommand,
                                              DefineSomaticTumorOnlyWGSOptions,
                                              options,
                                              SomaticTumorOnlyWGSMain,
                                              "Run somatic tumor-only WGS pipeline.",
                                              [](const CLI::App* const app, const CopyNumberCallerOptionsPtr& opts) {
                                                ConfigureSomaticTumorOnlyWGSOptions(app, *opts);
                                                CopyNumberCallerCliPreCallback(app, *opts);
                                              });

  cli::AddSubcommand<CopyNumberCallerOptions>(app,
                                              kGermlineNormalWGSSubcommand,
                                              DefineGermlineNormalWGSOptions,
                                              options,
                                              GermlineNormalWGSMain,
                                              "Run germline normal WGS pipeline.",
                                              [](const CLI::App* const app, const CopyNumberCallerOptionsPtr& opts) {
                                                ConfigureGermlineNormalWGSOptions(app, *opts);
                                                CopyNumberCallerCliPreCallback(app, *opts);
                                              });
}

s32 CopyNumberCallerCliMain(const s32 argc, char** const argv) {
  const auto app = cli::SetupDefaultCli(PROGRAM_NAME, VERSION);
  app->require_subcommand(1);
  auto options = std::make_shared<CopyNumberCallerOptions>();
  AddSubcommands(app.get(), options);
  return xoos::cli::RunCli(app.get(), argc, argv);
}
}  // namespace xoos::cnc
