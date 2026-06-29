#include "filter-variants-cli.h"

#include <xoos/cli/enum-option-util.h>
#include <xoos/cli/file-option-util.h>
#include <xoos/cli/thread-count-option-util.h>

#include "compute-bam-features/compute-bam-features-cli.h"
#include "compute-vcf-features/compute-vcf-features-cli.h"
#include "core/cli-option-names.h"
#include "core/command-line-info.h"
#include "util/cli-util.h"

namespace xoos::svc {

using cli::AddOptionalEnumOption;
using cli::AddThreadCountOption;
using enum Workflow;

// "custom" workflow not supported in filter_variants
constexpr std::array kSupportedWorkflows = {
    kGermline, kGermlineMultiSample, kTumorOnlyTe, kTumorNormalWgs, kGermlineTagging};

// CLI option group names
constexpr auto* kCliOptGroupBamFeatureExtraction{"BAM feature extraction"};
constexpr auto* kCliOptGroupVcfFeatureExtraction{"VCF feature extraction"};

// Default values for CLI options
constexpr u32 kDefaultMaxVcfRegionSizePerThread = 64'000;
constexpr auto* kDefaultSdChrName = "chr1";
constexpr auto* kDefaultVcfOutput = "output.vcf.gz";
constexpr u32 kAutomaticOutputVcfBufferSize = 0;
constexpr u32 kDefaultOutputVcfBufferSizeMultiplier = 8;

// Default pre-trained model path(s) for each workflow
constexpr auto* kDefaultGermlineSnvModelPath = "/resources/model-germline-sbxd-giraffe-snv.txt.gz";
constexpr auto* kDefaultGermlineIndelModelPath = "/resources/model-germline-sbxd-giraffe-indel.txt.gz";
constexpr auto* kDefaultGermlineMultiSampleSnvModelPath =
    "/resources/model-germline-sbxd-giraffe-multisample-snv.txt.gz";
constexpr auto* kDefaultGermlineMultiSampleIndelModelPath =
    "/resources/model-germline-sbxd-giraffe-multisample-indel.txt.gz";
constexpr auto* kDefaultTumorOnlyTeModelPath = "/resources/model-ffpe-bwa.txt.gz";
constexpr auto* kDefaultTumorNormalWgsModelPath = "/resources/model-somatic-tumor-normal-giraffe.txt.gz";

/**
 * @brief Helper function to define CLI options for the main application.
 * These options are either required or important in this submodule, or they are common to all other submodules in SVC.
 * @param app Main application pointer where CLI options will be defined.
 * @param params Shared pointer to CLI parameters to store parsed option values
 */
static void AddMainOptions(CLI::App* const app, const FilterVariantsParamPtr& params) {
  AddWarnAsErrorOption(app);

  AddThreadCountOption(app, cli_opt_name::kThreads, params->threads);

  app->add_option(cli_opt_name::kConfig, params->config_file, "Path to config JSON file")->check(CLI::ExistingFile);

  auto* const bam_opt = app->add_option(cli_opt_name::kBamInput,
                                        params->bam_files,
                                        "Path(s) to input BAM file(s), produced by GATK HaplotypeCaller/Mutect2")
                            ->required();
  CheckIndexedBamFile(bam_opt);

  auto* const vcf_opt = app->add_option(cli_opt_name::kVcfInput,
                                        params->vcf_file,
                                        "Path to input VCF file, produced by GATK Mutect2/HaplotypeCaller")
                            ->required();
  CheckIndexedVcfFile(vcf_opt);

  auto* const genome_opt =
      app->add_option(cli_opt_name::kGenome, params->genome, "Path to indexed FASTA file for reference genome")
          ->required();
  CheckIndexedFastaFile(genome_opt);

  auto* const regions_opt =
      app->add_option(cli_opt_name::kTargetRegions, params->bed_file, "Path to BED file for 0-based target regions");
  CheckBedFile(regions_opt);

  app->add_option(cli_opt_name::kOutputDir, params->output_dir, "Path to output directory")->default_val(".");

  const auto vcf_output_desc = fmt::format("Path to output VCF file, relative to {}", cli_opt_name::kOutputDir);
  app->add_option(cli_opt_name::kVcfOutput, params->vcf_output, vcf_output_desc)
      ->default_val(kDefaultVcfOutput)
      ->transform([&params](const std::string& path) { return (params->output_dir / fs::path(path)).string(); })
      ->check(CLI::NonexistentPath);
}

/**
 * @brief Helper function to add core CLI options for subcommands. CLI options added here are common across all
 * subcommands.
 * @param sub Subcommand application pointer where CLI options are to be added
 * @param params Shared pointer to CLI parameters to store parsed option values
 */
static void AddCoreOptions(CLI::App* const sub, const FilterVariantsParamPtr& params) {
  sub->add_option(cli_opt_name::kMaxBamRegionSizePerThread,
                  params->max_bam_region_size_per_thread,
                  "Maximum BAM region size per thread during feature extraction (inclusive)")
      ->default_val(kDefaultMaxBamRegionSizePerThread)
      ->check(CLI::PositiveNumber);

  sub->add_option(cli_opt_name::kMaxVcfRegionSizePerThread,
                  params->max_vcf_region_size_per_thread,
                  "Maximum VCF region size per thread during variant filtration (inclusive)")
      ->default_val(kDefaultMaxVcfRegionSizePerThread)
      ->check(CLI::PositiveNumber);

  const auto output_buffer_size_desc =
      fmt::format("Output buffer size for writing VCF file. '{}' sets automatically to {} * {}",
                  kAutomaticOutputVcfBufferSize,
                  kDefaultOutputVcfBufferSizeMultiplier,
                  cli_opt_name::kThreads);
  sub->add_option(cli_opt_name::kOutputVcfBufferSize, params->output_vcf_buffer_size, output_buffer_size_desc)
      ->default_val(kAutomaticOutputVcfBufferSize)
      ->check(CLI::NonNegativeNumber);

  sub->add_option(cli_opt_name::kSdChrName, params->sd_chr_name, "autosome name for sex determination")
      ->default_val(kDefaultSdChrName);

  auto* const par_bed_x_opt = sub->add_option(
      cli_opt_name::kParBedX, params->par_bed_x, "Path to chromosome X pseudoautosommal region BED file");
  CheckBedFile(par_bed_x_opt);

  auto* const par_bed_y_opt = sub->add_option(
      cli_opt_name::kParBedY, params->par_bed_y, "Path to chromosome Y pseudoautosommal region BED file");
  CheckBedFile(par_bed_y_opt);
}

/**
 * @brief Helper function to add tumor-normal-wgs specific CLI options.
 * @param app CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
static void AddTumorNormalWgsSpecificOptions(const CLI::App* const app, const FilterVariantsParamPtr& params) {
  CLI::App* const sub = app->get_subcommand(enum_util::FormatEnumName(kTumorNormalWgs));
  const auto defaults = SVCConfig(kTumorNormalWgs);

  // TODO: add an option for passing in a germline-tagging model and rename this option
  auto* const model_opt = sub->add_option(cli_opt_name::kModel, params->model)
                              ->description("Path to input model file")
                              ->default_val(kDefaultTumorNormalWgsModelPath);
  CheckNonEmptyFile(model_opt);

  sub->add_option(cli_opt_name::kTumorSampleName, params->tumor_sample_name, "tumor sample name for read groups")
      ->required();

  sub->add_option(cli_opt_name::kSnvMinMlScore,
                  params->somatic_tn_snv_ml_threshold,
                  "ML score threshold for Somatic TN for SNVs (inclusive). Variants must score at least this much to "
                  "be included")
      ->default_val(defaults.snv_min_ml_score)
      ->check(kCliRangeFraction);

  sub->add_option(cli_opt_name::kIndelMinMlScore,
                  params->somatic_tn_indel_ml_threshold,
                  "ML score threshold for Somatic TN for InDels (inclusive). Variants must score at least this much to "
                  "be included")
      ->default_val(defaults.indel_min_ml_score)
      ->check(kCliRangeFraction);

  sub->add_option(cli_opt_name::kMinTumorSupport,
                  params->tumor_support_threshold,
                  "Tumor support threshold for Somatic TN (inclusive). Variants must have at least this much support "
                  "to be included")
      ->default_val(defaults.min_tumor_support)
      ->check(CLI::NonNegativeNumber);
}

/**
 * @brief Helper function to add tumor-only-te specific CLI options.
 * @param app CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
static void AddTumorOnlyTeSpecificOptions(const CLI::App* const app, const FilterVariantsParamPtr& params) {
  CLI::App* const sub = app->get_subcommand(enum_util::FormatEnumName(kTumorOnlyTe));
  const auto defaults = SVCConfig(kTumorOnlyTe);

  auto* const model_opt = sub->add_option(cli_opt_name::kModel, params->model)
                              ->description("Path to input model file")
                              ->default_val(kDefaultTumorOnlyTeModelPath);
  CheckNonEmptyFile(model_opt);

  auto* const blocklist_opt = sub->add_option(cli_opt_name::kBlocklist,
                                              params->block_list,
                                              "Text file listing variants to skip. Variants are represented as "
                                              "`chr_pos_ref_alt`, one per line. Variant position is 1-based.");
  CheckNonEmptyFile(blocklist_opt);

  auto* const hotspot_opt = sub->add_option(
      cli_opt_name::kHotspotVcf, params->hotspot_list, "A VCF file containing a list of hotspot variants");
  CheckNonEmptyFile(hotspot_opt);

  auto* const forcecall_opt =
      sub->add_option(cli_opt_name::kForcecallBed,
                      params->forcecall_list,
                      "A 0-based BED file with forced call positions to be included in the output");
  CheckBedFile(forcecall_opt);

  sub->add_option(cli_opt_name::kMinAltCounts,
                  params->min_alt_counts,
                  "Minimum number of alt counts for retaining a variant (inclusive)")
      ->default_val(defaults.min_alt_counts)
      ->check(CLI::NonNegativeNumber);

  sub->add_option(cli_opt_name::kMinAf,
                  params->min_allele_freq_threshold,
                  "Minimum allowed variant allele frequency (inclusive)")
      ->default_val(defaults.min_af)
      ->check(kCliRangeFraction);

  sub->add_option(cli_opt_name::kMinPhasedAf,
                  params->min_phased_allele_freq,
                  "Minimum allowed phased allele frequency for a variant (inclusive)")
      ->default_val(defaults.min_phased_af)
      ->check(kCliRangeFraction);

  sub->add_option(cli_opt_name::kMaxPhasedAf,
                  params->max_phased_allele_freq,
                  "Maximum allowed phased allele frequency for a variant (inclusive)")
      ->default_val(defaults.max_phased_af)
      ->check(kCliRangeFraction);

  sub->add_option(cli_opt_name::kMinWeightedCounts,
                  params->weighted_counts_threshold,
                  "Threshold for weighted counts (inclusive). Variants must score at least this much to be included")
      ->default_val(defaults.min_weighted_counts)
      ->check(CLI::NonNegativeNumber);

  sub->add_option(
         cli_opt_name::kHotspotMinWeightedCounts,
         params->hotspot_weighted_counts_threshold,
         "Threshold for hotspot weighted counts (inclusive). Variants must score at least this much to be included")
      ->default_val(defaults.hotspot_min_weighted_counts)
      ->check(CLI::NonNegativeNumber);

  sub->add_option(cli_opt_name::kMinMlScore,
                  params->ml_threshold,
                  "Threshold for ML score (inclusive). Variants must score at least this much to be included")
      ->default_val(defaults.min_ml_score)
      ->check(kCliRangeFraction);

  sub->add_option(
         cli_opt_name::kHotspotMinMlScore,
         params->hotspot_ml_threshold,
         "Threshold for ML score in hotspots (inclusive). Variants must score at least this much to be included")
      ->default_val(defaults.hotspot_min_ml_score)
      ->check(kCliRangeFraction);

  sub->add_option(cli_opt_name::kPhased, params->phased, "Call phased variants in VCF")->default_val(defaults.phased);
}

/**
 * @brief Helper function to add germline-tagging specific CLI options.
 * @param app CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
static void AddGermlineTaggingSpecificOptions(const CLI::App* const app, const FilterVariantsParamPtr& params) {
  CLI::App* const sub = app->get_subcommand(enum_util::FormatEnumName(kGermlineTagging));

  // no default model for germline-tagging, user must provide one
  // TODO: add a default germline-tagging model once available
  auto* const opt =
      sub->add_option(cli_opt_name::kModel, params->model)->description("Path to input model file")->required();
  CheckNonEmptyFile(opt);
}

/**
 * @brief Helper function to add germline specific CLI options.
 * @param app CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
static void AddGermlineSpecificOptions(const CLI::App* const app, const FilterVariantsParamPtr& params) {
  CLI::App* const sub = app->get_subcommand(enum_util::FormatEnumName(kGermline));

  auto* const snv_model_opt = sub->add_option(cli_opt_name::kSnvModel, params->snv_model)
                                  ->description("Path to input SNV model file")
                                  ->default_val(kDefaultGermlineSnvModelPath);
  CheckNonEmptyFile(snv_model_opt);

  auto* const indel_model_opt = sub->add_option(cli_opt_name::kIndelModel, params->indel_model)
                                    ->description("Path to input indel model file")
                                    ->default_val(kDefaultGermlineIndelModelPath);
  CheckNonEmptyFile(indel_model_opt);
}

/**
 * @brief Helper function to add germline-multi-sample specific CLI options.
 * @param app CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
static void AddGermlineMultiSampleSpecificOptions(const CLI::App* const app, const FilterVariantsParamPtr& params) {
  CLI::App* const sub = app->get_subcommand(enum_util::FormatEnumName(kGermlineMultiSample));

  auto* const snv_model_opt = sub->add_option(cli_opt_name::kSnvModel, params->snv_model)
                                  ->description("Path to input SNV model file")
                                  ->default_val(kDefaultGermlineMultiSampleSnvModelPath);
  CheckNonEmptyFile(snv_model_opt);

  auto* const indel_model_opt = sub->add_option(cli_opt_name::kIndelModel, params->indel_model)
                                    ->description("Path to input indel model file")
                                    ->default_val(kDefaultGermlineMultiSampleIndelModelPath);
  CheckNonEmptyFile(indel_model_opt);
}

void filter_variants::DefineOptions(CLI::App* const app, FilterVariantsParamPtr& params) {
  AddMainOptions(app, params);

  // Add a subcommand for each workflow
  for (const Workflow workflow : kSupportedWorkflows) {
    const std::string name = enum_util::FormatEnumName(workflow);
    const std::string desc = fmt::format("Filter variants for the {} workflow", name);
    CLI::App* const sub = app->add_subcommand(name, desc)->fallthrough();
    // Do not apply force_callback() to subcommand options to avoid overwriting params set by other subcommands
    AddCoreOptions(sub, params);

    // Add options for VCF feature extraction
    for (auto* const opt : compute_vcf_features::AddSharedOptions(sub, params)) {
      opt->group(kCliOptGroupVcfFeatureExtraction);
    }

    const auto defaults = SVCConfig(workflow);

    // Add options for BAM feature extraction
    for (auto* const opt : compute_bam_features::AddSharedOptions(sub, params, defaults)) {
      opt->group(kCliOptGroupBamFeatureExtraction);
    }

    cli::AddEnumOption(sub,
                       cli_opt_name::kNormalizeFeatures,
                       params->normalize_features,
                       "Normalize read-depth related feature values",
                       defaults.normalize_features);
  }
  app->require_subcommand(kMinSubcommands, kMaxSubcommands);

  // Add workflow-specific CLI options for each subcommand
  AddTumorNormalWgsSpecificOptions(app, params);
  AddTumorOnlyTeSpecificOptions(app, params);
  AddGermlineTaggingSpecificOptions(app, params);
  AddGermlineSpecificOptions(app, params);
  AddGermlineMultiSampleSpecificOptions(app, params);

  // hide the `germline-tagging` subcommand by assigning it to an empty string group
  app->get_subcommand(enum_util::FormatEnumName(kGermlineTagging))->group("");
}

/**
 * @brief Helper function to apply workflow config to tumor-normal-wgs unique CLI parameters.
 * @param sub Subcommand CLI application pointer to check which options were specified via CLI
 * @param params CLI parameters shared pointer to apply defaults to
 */
static void ApplyConfigToTumorNormalWgsUniqueOptions(const CLI::App* const sub, const FilterVariantsParamPtr& params) {
  const bool is_not_tn_wgs = (params->config.workflow != kTumorNormalWgs);
  if (is_not_tn_wgs || sub->count(cli_opt_name::kSnvMinMlScore) == 0) {
    params->somatic_tn_snv_ml_threshold = params->config.snv_min_ml_score;
  }
  if (is_not_tn_wgs || sub->count(cli_opt_name::kIndelMinMlScore) == 0) {
    params->somatic_tn_indel_ml_threshold = params->config.indel_min_ml_score;
  }
  if (is_not_tn_wgs || sub->count(cli_opt_name::kMinTumorSupport) == 0) {
    params->tumor_support_threshold = params->config.min_tumor_support;
  }
}

/**
 * @brief Helper function to apply workflow config to tumor-only-te unique CLI parameters.
 * @param sub Subcommand CLI application pointer to check which options were specified via CLI
 * @param params CLI parameters shared pointer to apply defaults to
 */
static void ApplyConfigToTumorOnlyTeUniqueOptions(const CLI::App* const sub, const FilterVariantsParamPtr& params) {
  const bool is_not_to_te = (params->config.workflow != kTumorOnlyTe);
  if (is_not_to_te || sub->count(cli_opt_name::kMinAltCounts) == 0) {
    params->min_alt_counts = params->config.min_alt_counts;
  }
  if (is_not_to_te || sub->count(cli_opt_name::kMinAf) == 0) {
    params->min_allele_freq_threshold = params->config.min_af;
  }
  if (is_not_to_te || sub->count(cli_opt_name::kMinPhasedAf) == 0) {
    params->min_phased_allele_freq = params->config.min_phased_af;
  }
  if (is_not_to_te || sub->count(cli_opt_name::kMaxPhasedAf) == 0) {
    params->max_phased_allele_freq = params->config.max_phased_af;
  }
  if (is_not_to_te || sub->count(cli_opt_name::kMinWeightedCounts) == 0) {
    params->weighted_counts_threshold = params->config.min_weighted_counts;
  }
  if (is_not_to_te || sub->count(cli_opt_name::kHotspotMinWeightedCounts) == 0) {
    params->hotspot_weighted_counts_threshold = params->config.hotspot_min_weighted_counts;
  }
  if (is_not_to_te || sub->count(cli_opt_name::kMinMlScore) == 0) {
    params->ml_threshold = params->config.min_ml_score;
  }
  if (is_not_to_te || sub->count(cli_opt_name::kHotspotMinMlScore) == 0) {
    params->hotspot_ml_threshold = params->config.hotspot_min_ml_score;
  }
  if (is_not_to_te || sub->count(cli_opt_name::kPhased) == 0) {
    params->phased = params->config.phased;
  }
}

static void SetGermlineDefaultModelPaths(const CLI::App* const sub,
                                         const FilterVariantsParamPtr& params,
                                         const char* snv_default_path,
                                         const char* indel_default_path) {
  if (sub->count(cli_opt_name::kSnvModel) == 0) {
    params->snv_model = snv_default_path;
  }
  if (sub->count(cli_opt_name::kIndelModel) == 0) {
    params->indel_model = indel_default_path;
  }
}

static void SetNonGermlineDefaultModelPath(const CLI::App* const sub,
                                           const FilterVariantsParamPtr& params,
                                           const char* default_path) {
  if (sub->count(cli_opt_name::kModel) == 0) {
    params->model = default_path;
  }
}

/**
 * @brief Helper function to set default model paths based on workflow if not provided via CLI.
 * @param sub Subcommand CLI application pointer to check which options were specified via CLI
 * @param params CLI parameters shared pointer to apply defaults to
 */
static void SetDefaultModelPaths(const CLI::App* const sub, const FilterVariantsParamPtr& params) {
  switch (params->workflow) {
    case kGermline:
      SetGermlineDefaultModelPaths(sub, params, kDefaultGermlineSnvModelPath, kDefaultGermlineIndelModelPath);
      break;
    case kGermlineMultiSample:
      SetGermlineDefaultModelPaths(
          sub, params, kDefaultGermlineMultiSampleSnvModelPath, kDefaultGermlineMultiSampleIndelModelPath);
      break;
    case kTumorOnlyTe:
      SetNonGermlineDefaultModelPath(sub, params, kDefaultTumorOnlyTeModelPath);
      break;
    case kTumorNormalWgs:
      SetNonGermlineDefaultModelPath(sub, params, kDefaultTumorNormalWgsModelPath);
      break;
    default:
      SetNonGermlineDefaultModelPath(sub, params, "");
      break;
  }
}

/**
 * @brief Helper function to apply CLI defaults and workflow config presets to CLI parameters if the corresponding CLI
 * options were not specified.
 * @param sub Subcommand CLI application pointer to check which options were specified via CLI
 * @param params CLI parameters shared pointer to apply defaults to
 */
static void ApplyConfig(const CLI::App* const sub, const FilterVariantsParamPtr& params) {
  // set default model paths if not provided via CLI
  SetDefaultModelPaths(sub, params);

  // apply config to compute-bam-features options
  compute_bam_features::ApplyConfig(sub, params);

  if (sub->count(cli_opt_name::kNormalizeFeatures) == 0) {
    params->normalize_features = params->config.normalize_features;
  }

  // apply config to subcommand unique options
  ApplyConfigToTumorNormalWgsUniqueOptions(sub, params);
  ApplyConfigToTumorOnlyTeUniqueOptions(sub, params);
}

/**
 * @brief Preprocess the CLI options and validate the parameters.
 * @param app CLI application pointer
 * @param params CLI parameters to be validated
 * @throws CLI::ValidationError if the parameters are invalid
 */
void filter_variants::PreCallback(const cli::ConstAppPtr app, const FilterVariantsParamPtr& params) {
  params->command_line = GetCommandLineInfo(app);

  if (app->count(cli_opt_name::kVcfOutput) == 0) {
    // update VCF output path to be relative to output directory
    params->vcf_output = fs::path(params->output_dir) / params->vcf_output;
    if (fs::exists(params->vcf_output)) {
      throw CLI::ValidationError(fmt::format("Output VCF file '{}' already exists", params->vcf_output.string()));
    }
  }

  // Check which subcommand was used, set workflow and config accordingly, and apply config defaults as needed
  for (const Workflow workflow : kSupportedWorkflows) {
    const std::string name = enum_util::FormatEnumName(workflow);
    if (app->got_subcommand(name)) {
      params->workflow = workflow;
      params->config = JsonToConfig(params->config_file.value_or(fs::path{}), params->workflow);
      ApplyConfig(app->get_subcommand(name), params);
      break;
    }
  }

  if (params->output_vcf_buffer_size == kAutomaticOutputVcfBufferSize) {
    params->output_vcf_buffer_size = static_cast<u32>(params->threads * kDefaultOutputVcfBufferSizeMultiplier);
  }

  if (params->decode_yc != YcDecodeMethod::kNone && !IsDuplexProtocol(params->sequencing_protocol)) {
    throw CLI::ValidationError("Duplex sequencing protocol must be used when decoding YC tags");
  }
}

}  // namespace xoos::svc
