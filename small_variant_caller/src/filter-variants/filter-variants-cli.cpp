#include "filter-variants-cli.h"

#include <algorithm>
#include <array>
#include <string>
#include <string_view>

#include <fmt/format.h>

#include <xoos/cli/cli.h>
#include <xoos/cli/enum-option-util.h>
#include <xoos/cli/file-option-util.h>
#include <xoos/cli/thread-count-option-util.h>
#include <xoos/cli/validators/bed-validator.h>
#include <xoos/cli/validators/indexed-bam-validator.h>
#include <xoos/cli/validators/indexed-fasta-validator.h>
#include <xoos/cli/validators/indexed-vcf-validator.h>
#include <xoos/cli/validators/nonempty-file-validator.h>
#include <xoos/error/error.h>

#include "../core/workflow.h"
#include "../util/model-resolver.h"
#include "compute-bam-features/compute-bam-features-cli.h"
#include "compute-vcf-features/compute-vcf-features-cli.h"
#include "core/cli-option-names.h"
#include "core/germline-model-type.h"
#include "file_type_inference/sample-type-metrics.h"
#include "util/cli-util.h"
#include "util/model-resolver.h"

namespace xoos::svc {

using cli::AddThreadCountOption;
using enum Workflow;

// "custom" workflow not supported in filter_variants
constexpr std::array kSupportedWorkflows = {
    kGermlineWgs, kTumorOnlyTe, kTumorOnlyWgs, kTumorNormalWgs, kGermlineTagging, kGermlineMitochondria};

// CLI option group names
constexpr auto* kCliOptGroupBamFeatureExtraction{"BAM feature extraction"};
constexpr auto* kCliOptGroupVcfFeatureExtraction{"VCF feature extraction"};

// Default values for CLI options
constexpr u32 kDefaultMaxVcfRegionSizePerThread = 64'000;
constexpr auto* kDefaultSdChrName = "chr1";
constexpr auto* kDefaultVcfOutput = "output.vcf.gz";
constexpr auto* kDefaultMetricsOutput = "metrics/metrics.tsv";
// Default pre-trained model path(s) for each workflow
constexpr auto* kDefaultTumorOnlyTeModelPath = "/resources/model-ffpe-bwa.txt.gz";

namespace {
/**
 * @brief Helper function to add core CLI options for subcommands. CLI options added here are common across all
 * subcommands.
 * @param sub Subcommand application pointer where CLI options are to be added
 * @param params Shared pointer to CLI parameters to store parsed option values
 */
void AddCommonOptions(CLI::App* const sub, const FilterVariantsParamPtr& params) {
  AddWarnAsErrorOption(sub);

  AddThreadCountOption(sub, cli_opt_name::kThreads, params->threads);

  sub->add_option(cli_opt_name::kConfig, params->config_file, "Path to config JSON file")->check(CLI::ExistingFile);

  auto* const bam_opt = sub->add_option(cli_opt_name::kBamInput,
                                        params->bam_files,
                                        "Path(s) to input BAM file(s), produced by GATK HaplotypeCaller/Mutect2")
                            ->required();
  bam_opt->check(cli::IndexedBamFileValidator());

  auto* const vcf_opt = sub->add_option(cli_opt_name::kVcfInput,
                                        params->vcf_file,
                                        "Path to input VCF file, produced by GATK Mutect2/HaplotypeCaller")
                            ->required();
  vcf_opt->check(cli::IndexedVcfFileValidator());

  auto* const genome_opt =
      sub->add_option(cli_opt_name::kGenome, params->genome, "Path to indexed FASTA file for reference genome")
          ->required();
  genome_opt->check(cli::IndexedFastaFileValidator());

  auto* const regions_opt =
      sub->add_option(cli_opt_name::kTargetRegions,
                      params->bed_file,
                      "Path to uncompressed or gzip/BGZF-compressed (.gz) BED file for 0-based target regions");
  regions_opt->check(cli::BedFileValidator());

  sub->add_option(cli_opt_name::kOutputDir, params->output_dir, "Path to output directory")->default_val(".");

  const auto vcf_output_desc = fmt::format("Path to output VCF file, relative to {}", cli_opt_name::kOutputDir);
  sub->add_option(cli_opt_name::kVcfOutput, params->vcf_output, vcf_output_desc)
      ->default_val(kDefaultVcfOutput)
      ->transform([&params](const std::string& path) { return (params->output_dir / fs::path(path)).string(); })
      ->check(CLI::NonexistentPath);

  const auto metrics_output_desc =
      fmt::format("Path to output QC metrics TSV file, relative to {}", cli_opt_name::kOutputDir);
  sub->add_option(cli_opt_name::kOutputMetrics, params->metrics_output, metrics_output_desc)
      ->default_val(kDefaultMetricsOutput)
      ->transform([&params](const std::string& path) { return (params->output_dir / fs::path(path)).string(); });

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
}

/**
 * @brief Helper function to add sex-chromosome CLI options. These cover sex determination and PAR-aware variant calling
 * on the sex chromosomes. They are not applicable to the `germline-mitochondria` workflow, which filters mitochondrial
 * variants on AF/depth.
 * @param sub Subcommand application pointer where CLI options are to be added
 * @param params Shared pointer to CLI parameters to store parsed option values
 */
void AddSexChromosomeOptions(CLI::App* const sub, const FilterVariantsParamPtr& params) {
  sub->add_option(cli_opt_name::kSdChrName, params->sd_chr_name, "autosome name for sex determination")
      ->default_val(kDefaultSdChrName);

  auto* const par_bed_x_opt = sub->add_option(
      cli_opt_name::kParBedX,
      params->par_bed_x,
      "Path to uncompressed or gzip/BGZF-compressed (.gz) chromosome X pseudoautosomal region BED file");
  par_bed_x_opt->check(cli::BedFileValidator());

  auto* const par_bed_y_opt = sub->add_option(
      cli_opt_name::kParBedY,
      params->par_bed_y,
      "Path to uncompressed or gzip/BGZF-compressed (.gz) chromosome Y pseudoautosomal region BED file");
  par_bed_y_opt->check(cli::BedFileValidator());
}

/**
 * @brief Helper function to add CLI options for outputting SHAP values.
 * @param sub Subcommand application pointer where CLI options are to be added
 * @param params Shared pointer to CLI parameters to store parsed option values
 */
void AddOutputShapValueTsvOption(CLI::App* const sub, const FilterVariantsParamPtr& params) {
  sub->add_option(cli_opt_name::kOutputShapValueTsv,
                  params->shap_value_tsv,
                  "Path to the SHAP value output TSV file, relative to output directory. Intended to be applied to "
                  "target regions due to slow runtime.")
      ->needs(cli_opt_name::kTargetRegions);
}

/**
 * @brief Helper function to add CLI options for outputting SHAP values for SNVs and indels separately. This is intended
 * to be used for germline workflow where separate SHAP value output for SNVs and indels is applicable.
 * @param sub Subcommand application pointer where CLI options are to be added
 * @param params Shared pointer to CLI parameters to store parsed option values
 */
void AddOutputSnvIndelShapValueTsvOptions(CLI::App* const sub, const FilterVariantsParamPtr& params) {
  sub->add_option(cli_opt_name::kOutputSnvShapValueTsv,
                  params->snv_shap_value_tsv,
                  "Path to the SHAP value output TSV file for SNVs, relative to output directory. Intended to be "
                  "applied to target regions due to slow runtime.")
      ->needs(cli_opt_name::kTargetRegions);

  sub->add_option(cli_opt_name::kOutputIndelShapValueTsv,
                  params->indel_shap_value_tsv,
                  "Path to the SHAP value output TSV file for indels, relative to output directory. Intended to be "
                  "applied to target regions due to slow runtime.")
      ->needs(cli_opt_name::kTargetRegions);
}

/**
 * @brief Add an --aligner option that accepts only 'bwa' (auto-select) and 'custom' (explicit --model).
 * @details The tumor-normal-wgs and tumor-only-wgs workflows only ship pre-trained BWA models, so 'bwa'
 *          and 'custom' are the only valid choices. 'giraffe' is intentionally not accepted — it is
 *          rejected at parse time and is not advertised in the help text — to avoid the confusing
 *          combination of '--aligner giraffe' with an explicit '--model'.
 * @param sub Subcommand application pointer where the option is added.
 * @param params CLI parameters shared pointer to store the parsed aligner value.
 * @param workflow_name Workflow label used in the validation error message.
 * @param default_aligner Aligner selected when the user does not pass --aligner. 'bwa' auto-resolves the
 *        model; 'custom' requires an explicit --model.
 */
void AddAlignerOption(CLI::App* const sub,
                      const FilterVariantsParamPtr& params,
                      const std::string_view workflow_name,
                      const AlignerType default_aligner) {
  static constexpr std::array kSupportedAligners = {AlignerType::kBwa, AlignerType::kCustom};
  const auto supported = fmt::format(
      "{{{}, {}}}", enum_util::FormatEnumName(AlignerType::kBwa), enum_util::FormatEnumName(AlignerType::kCustom));
  const std::string workflow{workflow_name};
  sub->add_option_function<std::string>(
         cli_opt_name::kAligner,
         [params, supported, workflow](const std::string& input) {
           const auto parsed = enum_util::ParseEnumName<AlignerType>(input);
           if (!parsed || std::ranges::find(kSupportedAligners, *parsed) == kSupportedAligners.end()) {
             throw CLI::ValidationError(
                 fmt::format("{}: Check {} value in {} FAILED. 'giraffe' is not supported for "
                             "the {} workflow because no pre-trained giraffe model "
                             "is available; use 'custom' with an explicit {} instead.",
                             cli_opt_name::kAligner,
                             input,
                             supported,
                             workflow,
                             cli_opt_name::kModel));
           }
           params->aligner = *parsed;
         },
         fmt::format("Aligner used for BAM alignment. Set to 'bwa' to auto-select the model based on --sample-type, "
                     "or set to 'custom' and provide --model explicitly. Accepted values: {}",
                     supported))
      ->run_callback_for_default()
      ->default_val(enum_util::FormatEnumName(default_aligner));
}

/**
 * @brief Add a workflow-specific sample-type option without mutating shared parameters during option registration.
 * @details All filter_variants subcommands share one FilterVariantsParam. AddEnumOption applies its default callback
 * immediately, which would let the last registered workflow overwrite sample_type even when it is not selected.
 * The active tumor workflow applies its default in ResolveTumorNormalWgsModelPath or ResolveTumorOnlyWgsModelPath.
 * @param sub Subcommand on which to register the option.
 * @param params Shared filter_variants parameters populated for explicitly provided values.
 * @param description CLI help text for the option.
 * @param default_sample_type Default displayed for the selected tumor workflow.
 */
void AddSampleTypeOption(CLI::App* const sub,
                         const FilterVariantsParamPtr& params,
                         const std::string& description,
                         const SampleType default_sample_type) {
  cli::AddOptionalEnumOption(sub, cli_opt_name::kSampleType, params->sample_type, description)
      ->run_callback_for_default(false)
      ->default_val(enum_util::FormatEnumName(default_sample_type));
}

/**
 * @brief Add controls shared by workflows that infer sample type from the input VCF.
 */
void AddSampleTypeDetectionOptions(CLI::App* const sub, const FilterVariantsParamPtr& params) {
  sub->add_option(cli_opt_name::kSampleTypeMaxAf,
                  params->sample_type_max_af,
                  "Maximum allele frequency for SNV error observations used by sample-type auto-detection")
      ->default_val(kDefaultMaxAfForSampleTypeDetection)
      ->check(kCliRangeFraction);

  sub->add_option(cli_opt_name::kSampleTypeMaxFfpePhred,
                  params->sample_type_max_ffpe_phred,
                  "Maximum VCF phred for classifying auto-detected sample type as FFPE")
      ->default_val(kDefaultMaxFfpePhredForSampleTypeDetection)
      ->check(CLI::NonNegativeNumber);
}

/**
 * @brief Helper function to add tumor-normal-wgs specific CLI options.
 * @param app CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
void AddTumorNormalWgsSpecificOptions(const CLI::App* const app, const FilterVariantsParamPtr& params) {
  CLI::App* const sub = app->get_subcommand(enum_util::FormatEnumName(kTumorNormalWgs));
  const auto defaults = SVCConfig(kTumorNormalWgs);

  // Model is auto-resolved when omitted. Explicit --model keeps the custom-model behavior unless --aligner is
  // explicitly provided.
  auto* const model_opt = sub->add_option(cli_opt_name::kModel, params->model)
                              ->description(
                                  "Path to input model file. When omitted, the BWA model is automatically "
                                  "selected based on the inferred or explicitly provided --sample-type");
  model_opt->check(cli::NonEmptyFileValidator());

  AddSampleTypeOption(sub,
                      params,
                      "Sample preparation type for model selection. Ignored when --model is provided. When omitted "
                      "during BWA auto-selection, the sample type is inferred from the input VCF",
                      SampleType::kAutoDetect);

  AddSampleTypeDetectionOptions(sub, params);

  AddAlignerOption(sub, params, enum_util::FormatEnumName(kTumorNormalWgs), AlignerType::kCustom);

  sub->add_option(cli_opt_name::kTumorSampleName, params->tumor_sample_name, "tumor sample name for read groups")
      ->required();

  sub->add_option(
         cli_opt_name::kSnvMinMlScore, params->snv_min_ml_score, "Minimum ML score threshold (inclusive) for SNVs")
      ->default_val(defaults.snv_min_ml_score)
      ->check(kCliRangeFraction);

  sub->add_option(cli_opt_name::kIndelMinMlScore,
                  params->indel_min_ml_score,
                  "Minimum ML score threshold (inclusive) for indels")
      ->default_val(defaults.indel_min_ml_score)
      ->check(kCliRangeFraction);

  sub->add_option(cli_opt_name::kMinTumorSupport,
                  params->min_tumor_support,
                  "Minimum support threshold (inclusive) for variants in the tumor sample")
      ->default_val(defaults.min_tumor_support)
      ->check(CLI::NonNegativeNumber);

  sub->add_option(cli_opt_name::kMaxNormalSupport,
                  params->max_normal_support,
                  "Maximum support threshold (inclusive) for variants in the normal sample")
      ->default_val(defaults.max_normal_support)
      ->check(CLI::NonNegativeNumber);

  sub->add_option(cli_opt_name::kMinTumorAf,
                  params->min_tumor_af,
                  "Minimum allele frequency threshold (inclusive) for variants in the tumor sample")
      ->default_val(defaults.min_tumor_af)
      ->check(kCliRangeFraction);

  // Relaxed thresholds for variants present in the catalog VCF. These require --catalog-vcf and,
  // when unset, fall back to their standard counterparts (no relaxation).
  sub->add_option(cli_opt_name::kCatalogMinTumorAf,
                  params->catalog_min_tumor_af,
                  "Minimum allele frequency threshold (inclusive) for tumor-sample variants present in the "
                  "catalog VCF. Requires --catalog-vcf. Defaults to --min-tumor-af when unset")
      ->needs(cli_opt_name::kCatalogVcf)
      ->check(kCliRangeFraction);

  sub->add_option(cli_opt_name::kCatalogSnvMinMlScore,
                  params->catalog_snv_min_ml_score,
                  "Minimum ML score threshold (inclusive) for SNVs present in the catalog VCF. "
                  "Requires --catalog-vcf. Defaults to --snv-min-ml-score when unset")
      ->needs(cli_opt_name::kCatalogVcf)
      ->check(kCliRangeFraction);

  sub->add_option(cli_opt_name::kCatalogIndelMinMlScore,
                  params->catalog_indel_min_ml_score,
                  "Minimum ML score threshold (inclusive) for indels present in the catalog VCF. "
                  "Requires --catalog-vcf. Defaults to --indel-min-ml-score when unset")
      ->needs(cli_opt_name::kCatalogVcf)
      ->check(kCliRangeFraction);

  sub->add_option(cli_opt_name::kMinDpRatio,
                  params->min_dp_ratio,
                  "Minimum ratio (inclusive) of DP to chromosome's median DP for variants in the tumor sample")
      ->default_val(defaults.min_dp_ratio)
      ->check(CLI::NonNegativeNumber);

  sub->add_option(cli_opt_name::kMaxIndelSize, params->max_indel_size, "Maximum indel size in base pairs (inclusive)")
      ->default_val(defaults.max_indel_size)
      ->check(CLI::PositiveNumber);

  AddOutputShapValueTsvOption(sub, params);
}

/**
 * @brief Helper function to add tumor-only-wgs specific CLI options.
 *
 * Single-sample somatic filtering with a single binary model and SNV/indel score thresholds. A single-sample
 * alt-support gate (`--min-support`) and AF gate (`--min-af`) are applied before ML scoring. The
 * Panel of Normals VCF
 * (`--pon-vcf`) is required, since PON features are part of the tumor-only-wgs scoring schema.
 *
 * @param app CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
void AddTumorOnlyWgsSpecificOptions(const CLI::App* const app, const FilterVariantsParamPtr& params) {
  CLI::App* const sub = app->get_subcommand(enum_util::FormatEnumName(kTumorOnlyWgs));
  const auto defaults = SVCConfig(kTumorOnlyWgs);

  // Model is auto-resolved by default (aligner defaults to 'bwa'). When --aligner is set to 'custom',
  // the model must be provided explicitly via --model.
  auto* const model_opt = sub->add_option(cli_opt_name::kModel, params->model)
                              ->description(
                                  "Path to input model file. Auto-selected by default based on --sample-type. "
                                  "Required when --aligner is set to 'custom'");
  model_opt->check(cli::NonEmptyFileValidator());

  AddSampleTypeOption(sub,
                      params,
                      "Sample preparation type for model selection. Ignored when --model is provided. When omitted "
                      "during BWA auto-selection, the sample type is inferred from the input VCF",
                      SampleType::kAutoDetect);

  AddSampleTypeDetectionOptions(sub, params);

  AddAlignerOption(sub, params, enum_util::FormatEnumName(kTumorOnlyWgs), AlignerType::kBwa);

  sub->add_option(
         cli_opt_name::kSnvMinMlScore, params->snv_min_ml_score, "Minimum ML score threshold (inclusive) for SNVs")
      ->default_val(defaults.snv_min_ml_score)
      ->check(kCliRangeFraction);

  sub->add_option(cli_opt_name::kIndelMinMlScore,
                  params->indel_min_ml_score,
                  "Minimum ML score threshold (inclusive) for indels")
      ->default_val(defaults.indel_min_ml_score)
      ->check(kCliRangeFraction);

  sub->add_option(
         cli_opt_name::kMinSupport, params->min_support, "Minimum alt support threshold (inclusive) for variants")
      ->default_val(defaults.min_support)
      ->check(CLI::NonNegativeNumber);

  sub->add_option(cli_opt_name::kMinAf, params->min_af, "Minimum allele frequency threshold (inclusive) for variants")
      ->default_val(defaults.min_af)
      ->check(kCliRangeFraction);

  sub->add_option(cli_opt_name::kMaxIndelSize, params->max_indel_size, "Maximum indel size in base pairs (inclusive)")
      ->default_val(defaults.max_indel_size)
      ->check(CLI::PositiveNumber);

  AddOutputShapValueTsvOption(sub, params);
}

/**
 * @brief Helper function to add tumor-only-te specific CLI options.
 * @param app CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
void AddTumorOnlyTeSpecificOptions(const CLI::App* const app, const FilterVariantsParamPtr& params) {
  CLI::App* const sub = app->get_subcommand(enum_util::FormatEnumName(kTumorOnlyTe));
  const auto defaults = SVCConfig(kTumorOnlyTe);

  auto* const model_opt = sub->add_option(cli_opt_name::kModel, params->model)
                              ->description("Path to input model file")
                              ->default_val(kDefaultTumorOnlyTeModelPath);
  model_opt->check(cli::NonEmptyFileValidator());

  auto* const blocklist_opt = sub->add_option(cli_opt_name::kBlocklist,
                                              params->block_list,
                                              "Text file listing variants to skip. Variants are represented as "
                                              "`chr_pos_ref_alt`, one per line. Variant position is 1-based.");
  blocklist_opt->check(cli::NonEmptyFileValidator());

  auto* const hotspot_opt = sub->add_option(
      cli_opt_name::kHotspotVcf, params->hotspot_list, "A VCF file containing a list of hotspot variants");
  hotspot_opt->check(cli::NonEmptyFileValidator());

  auto* const forcecall_opt = sub->add_option(cli_opt_name::kForcecallBed,
                                              params->forcecall_list,
                                              "An uncompressed or gzip/BGZF-compressed (.gz) 0-based BED file with "
                                              "forced call positions to include in the output");
  forcecall_opt->check(cli::BedFileValidator());

  sub->add_option(cli_opt_name::kMinAltCounts,
                  params->min_alt_counts,
                  "Minimum number of alt counts for retaining a variant (inclusive)")
      ->default_val(defaults.min_alt_counts)
      ->check(CLI::NonNegativeNumber);

  sub->add_option(cli_opt_name::kMinAf, params->min_af, "Minimum allowed variant allele frequency (inclusive)")
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

  AddOutputShapValueTsvOption(sub, params);
}

/**
 * @brief Helper function to add germline-tagging specific CLI options.
 * @param app CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
void AddGermlineTaggingSpecificOptions(const CLI::App* const app, const FilterVariantsParamPtr& params) {
  CLI::App* const sub = app->get_subcommand(enum_util::FormatEnumName(kGermlineTagging));

  // no default model for germline-tagging, user must provide one
  auto* const opt =
      sub->add_option(cli_opt_name::kModel, params->model)->description("Path to input model file")->required();
  opt->check(cli::NonEmptyFileValidator());

  AddOutputShapValueTsvOption(sub, params);
}

/**
 * @brief Add germline WGS model-selection options.
 * @param sub Subcommand application pointer where CLI options are to be added
 * @param params Shared pointer to CLI parameters to store parsed option values
 */
void AddGermlineModelSelectionOptions(CLI::App* const sub, const FilterVariantsParamPtr& params) {
  cli::AddEnumOption(sub,
                     cli_opt_name::kAligner,
                     params->aligner,
                     "Aligner used for BAM alignment. Set to 'bwa' or 'giraffe' to auto-select models, "
                     "or leave as 'custom' and provide --snv-model and --indel-model explicitly",
                     AlignerType::kCustom);

  cli::AddEnumOption(
      sub,
      cli_opt_name::kLibraryType,
      params->library_type,
      "Library preparation type for config and model selection; auto-detect inspects BAM metadata and read content",
      LibraryType::kAutoDetect);

  cli::AddEnumOption(sub,
                     cli_opt_name::kModelType,
                     params->germline_model_type,
                     "Optional bundled germline WGS model family to use for any model path not supplied manually. "
                     "Defaults to 'large'",
                     GermlineModelType::kLarge);

  auto* const snv_model_opt = sub->add_option(cli_opt_name::kSnvModel, params->snv_model)
                                  ->description("Path to input SNV model file. Required when --aligner is 'custom'");
  snv_model_opt->check(cli::NonEmptyFileValidator());

  auto* const indel_model_opt =
      sub->add_option(cli_opt_name::kIndelModel, params->indel_model)
          ->description("Path to input indel model file. Required when --aligner is 'custom'");
  indel_model_opt->check(cli::NonEmptyFileValidator());
}

/**
 * @brief Helper function to add germline specific CLI options.
 * @param app CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
void AddGermlineSpecificOptions(const CLI::App* const app, const FilterVariantsParamPtr& params) {
  CLI::App* const sub = app->get_subcommand(enum_util::FormatEnumName(kGermlineWgs));

  AddGermlineModelSelectionOptions(sub, params);
  AddOutputSnvIndelShapValueTsvOptions(sub, params);
}

/**
 * @brief Helper function to add mitochondria specific CLI options.
 *  Mitochondria uses AF/depth thresholds instead of ML-based filtering.
 * @param app CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
void AddMitochondriaSpecificOptions(const CLI::App* const app, const FilterVariantsParamPtr& params) {
  CLI::App* const sub = app->get_subcommand(enum_util::FormatEnumName(kGermlineMitochondria));
  const auto defaults = SVCConfig(kGermlineMitochondria);

  auto* const blocklist_bed_opt =
      sub->add_option(cli_opt_name::kBlocklistBed,
                      params->blocklist_bed,
                      "Path to uncompressed or gzip/BGZF-compressed (.gz) BED file of 0-based, half-open regions whose "
                      "variant starts are tagged blocklist");
  blocklist_bed_opt->check(cli::BedFileValidator());

  sub->add_option(cli_opt_name::kMinAf, params->min_af, "Minimum allele frequency threshold for variant to pass")
      ->default_val(defaults.min_af)
      ->check(kCliRangeFraction);

  sub->add_option(cli_opt_name::kMinDp, params->min_dp, "Minimum depth threshold for variant to pass")
      ->default_val(defaults.min_dp)
      ->check(CLI::NonNegativeNumber);

  sub->add_option(cli_opt_name::kMinHomAltAf,
                  params->min_hom_alt_af,
                  "Minimum allele frequency at or above which a variant is genotyped as homozygous ALT (1/1)")
      ->default_val(defaults.min_hom_alt_af)
      ->check(kCliRangeFraction);
}

/**
 * @brief Resolve and render germline WGS profile defaults before the common CLI logs the invocation.
 * @details The root parse-complete callback runs before the module pre-callback. Rendering a local resolved profile
 *          here prevents the `Running:` log from showing duplex registration placeholders for simplex runs.
 */
void SyncGermlineWgsRenderedDefaults(CLI::App* const sub, const FilterVariantsParamPtr& params) {
  params->library_type = ResolveLibraryType(params->library_type, params->bam_files);
  const auto config =
      JsonToConfig(params->config_file.value_or(fs::path{}), GermlineWgsConfigProfile(params->library_type));
  SetRenderedDefault(sub, cli_opt_name::kLibraryType, enum_util::FormatEnumName(params->library_type));
  compute_bam_features::SyncRenderedConfigDefaults(sub, config);
  SetRenderedDefault(sub, cli_opt_name::kNormalizeFeatures, enum_util::FormatEnumName(config.normalize_features));
}

}  // namespace

void filter_variants::DefineOptions(CLI::App* const app, FilterVariantsParamPtr& params) {
  // Add a subcommand for each workflow
  for (const Workflow workflow : kSupportedWorkflows) {
    const std::string name = enum_util::FormatEnumName(workflow);
    const std::string desc = fmt::format("Filter variants for the {} workflow", name);
    CLI::App* const sub = app->add_subcommand(name, desc)->fallthrough();
    // Do not apply force_callback() to subcommand options to avoid overwriting params set by other subcommands
    AddCommonOptions(sub, params);

    // The germline-mitochondria workflow filters on AF/depth only; it does not perform sex determination or
    // ML-feature normalization, so the corresponding options are omitted from that subcommand.
    const bool is_mitochondria = (workflow == kGermlineMitochondria);
    if (!is_mitochondria) {
      AddSexChromosomeOptions(sub, params);
    }

    // Add options for VCF feature extraction
    for (auto* const opt : compute_vcf_features::AddSharedOptions(sub, params)) {
      if (opt->get_name() == cli_opt_name::kPonVcf) {
        // --pon-vcf is required for tumor-only-wgs. It is optional for tumor-normal-wgs, where omitted PON features
        // retain their zero defaults during scoring. It remains hidden for workflows without PON scoring features.
        if (workflow == kTumorOnlyWgs) {
          opt->group(kCliOptGroupVcfFeatureExtraction)->required();
        } else if (workflow == kTumorNormalWgs) {
          opt->group(kCliOptGroupVcfFeatureExtraction);
        } else {
          opt->group("");
        }
      } else {
        opt->group(kCliOptGroupVcfFeatureExtraction);
      }
    }

    // Germline defaults cannot be final until BAM inspection or --library-type resolve the profile. Duplex values are
    // used only to register CLI11 options and are replaced by ApplyConfig before filtering begins.
    const auto defaults = workflow == kGermlineWgs ? SVCConfig(ConfigProfile::kGermlineWgsDuplex) : SVCConfig(workflow);

    // Add options for BAM feature extraction
    for (auto* const opt : compute_bam_features::AddSharedOptions(sub, params, defaults)) {
      opt->group(kCliOptGroupBamFeatureExtraction);
    }

    if (!is_mitochondria) {
      cli::AddEnumOption(sub,
                         cli_opt_name::kNormalizeFeatures,
                         params->normalize_features,
                         "Normalize read-depth related feature values",
                         defaults.normalize_features);
    }
  }
  app->require_subcommand(kMinSubcommands, kMaxSubcommands);

  // Add workflow-specific CLI options for each subcommand
  AddTumorNormalWgsSpecificOptions(app, params);
  AddTumorOnlyWgsSpecificOptions(app, params);
  AddTumorOnlyTeSpecificOptions(app, params);
  AddGermlineTaggingSpecificOptions(app, params);
  AddGermlineSpecificOptions(app, params);
  AddMitochondriaSpecificOptions(app, params);

  auto* const germline_sub = app->get_subcommand(enum_util::FormatEnumName(kGermlineWgs));
  germline_sub->parse_complete_callback(
      [germline_sub, params] { SyncGermlineWgsRenderedDefaults(germline_sub, params); });

  // hide the `germline-tagging` subcommand by assigning it to an empty string group
  app->get_subcommand(enum_util::FormatEnumName(kGermlineTagging))->group("");
}

namespace {

/**
 * @brief Check whether a CLI option was explicitly supplied.
 * @details CLI11's count can include defaulted/shared option state in this command setup. This helper
 * checks the option parse result directly so model, aligner, sample-type, and threshold resolution can
 * distinguish user-provided values from defaults.
 * @param sub Parsed subcommand containing the option.
 * @param option_name CLI option flag to query.
 * @return True when the option exists and has parsed occurrences; otherwise false.
 */
bool OptionProvided(const CLI::App* const sub, const std::string_view option_name) {
  const auto* const option = sub->get_option_no_throw(std::string{option_name});
  return option != nullptr && !option->empty();
}

/**
 * @brief Check whether a CLI option was explicitly supplied on a subcommand.
 * @details Returns false when the option is not registered on the subcommand. This is useful for
 * workflow-specific options because CLI11's count() reports an error for absent option names.
 * @param sub Parsed subcommand to query.
 * @param option_name CLI option flag to query.
 * @return True when the option exists and has parsed occurrences; otherwise false.
 */
bool OptionCounted(const CLI::App* const sub, const std::string_view option_name) {
  return sub->get_option_no_throw(std::string{option_name}) != nullptr && sub->count(std::string{option_name}) > 0;
}

/**
 * @brief Helper function to apply workflow config to tumor-normal-wgs unique CLI parameters.
 * @param sub Subcommand CLI application pointer to check which options were specified via CLI
 * @param params CLI parameters shared pointer to apply defaults to
 */
void ApplyConfigToTumorNormalWgsUniqueOptions(const CLI::App* const sub, const FilterVariantsParamPtr& params) {
  if (sub->count(cli_opt_name::kMinTumorSupport) == 0) {
    params->min_tumor_support = params->config.min_tumor_support;
  }
  if (sub->count(cli_opt_name::kMaxNormalSupport) == 0) {
    params->max_normal_support = params->config.max_normal_support;
  }
  if (sub->count(cli_opt_name::kMinDpRatio) == 0) {
    params->min_dp_ratio = params->config.min_dp_ratio;
  }
}

/**
 * @brief Apply config defaults to ML-score and max-indel-size options shared by the tumor-normal-wgs and
 * tumor-only-wgs subcommands.
 *
 * These options are registered on both subcommands, so they are handled in one place to avoid the two
 * workflow-specific handlers clobbering each other's parsed CLI values. The value is reset to the config default
 * unless the active subcommand explicitly provided it. `sub->count()` is only consulted for subcommands that
 * register the option (tumor-normal-wgs and tumor-only-wgs); for any other workflow the config default is used.
 * @param sub Active subcommand CLI application pointer
 * @param params CLI parameters shared pointer to apply defaults to
 */
void ApplyConfigToSharedMlScoreOptions(const CLI::App* const sub, const FilterVariantsParamPtr& params) {
  const Workflow workflow = params->config.workflow;
  const bool exposes_options = (workflow == kTumorNormalWgs || workflow == kTumorOnlyWgs);
  if (!exposes_options || sub->count(cli_opt_name::kSnvMinMlScore) == 0) {
    params->snv_min_ml_score = params->config.snv_min_ml_score;
  }
  if (!exposes_options || sub->count(cli_opt_name::kIndelMinMlScore) == 0) {
    params->indel_min_ml_score = params->config.indel_min_ml_score;
  }
  if (workflow != kTumorNormalWgs || !OptionCounted(sub, cli_opt_name::kMinTumorAf)) {
    params->min_tumor_af = params->config.min_tumor_af;
  }
  if (!exposes_options || sub->count(cli_opt_name::kMaxIndelSize) == 0) {
    params->max_indel_size = params->config.max_indel_size;
  }
}

/**
 * @brief Helper function to apply workflow config to tumor-only-wgs unique CLI parameters.
 * @param sub Subcommand CLI application pointer to check which options were specified via CLI
 * @param params CLI parameters shared pointer to apply defaults to
 */
void ApplyConfigToTumorOnlyWgsUniqueOptions(const CLI::App* const sub, const FilterVariantsParamPtr& params) {
  if (sub->count(cli_opt_name::kMinSupport) == 0) {
    params->min_support = params->config.min_support;
  }
  if (!OptionCounted(sub, cli_opt_name::kMinAf)) {
    params->min_af = params->config.min_af;
  }
}

/**
 * @brief Helper function to apply workflow config to tumor-only-te unique CLI parameters.
 * @param sub Subcommand CLI application pointer to check which options were specified via CLI
 * @param params CLI parameters shared pointer to apply defaults to
 */
void ApplyConfigToTumorOnlyTeUniqueOptions(const CLI::App* const sub, const FilterVariantsParamPtr& params) {
  if (sub->count(cli_opt_name::kMinAltCounts) == 0) {
    params->min_alt_counts = params->config.min_alt_counts;
  }
  if (sub->count(cli_opt_name::kMinAf) == 0) {
    params->min_af = params->config.min_af;
  }
  if (sub->count(cli_opt_name::kMinPhasedAf) == 0) {
    params->min_phased_allele_freq = params->config.min_phased_af;
  }
  if (sub->count(cli_opt_name::kMaxPhasedAf) == 0) {
    params->max_phased_allele_freq = params->config.max_phased_af;
  }
  if (sub->count(cli_opt_name::kMinWeightedCounts) == 0) {
    params->weighted_counts_threshold = params->config.min_weighted_counts;
  }
  if (sub->count(cli_opt_name::kHotspotMinWeightedCounts) == 0) {
    params->hotspot_weighted_counts_threshold = params->config.hotspot_min_weighted_counts;
  }
  if (sub->count(cli_opt_name::kMinMlScore) == 0) {
    params->ml_threshold = params->config.min_ml_score;
  }
  if (sub->count(cli_opt_name::kHotspotMinMlScore) == 0) {
    params->hotspot_ml_threshold = params->config.hotspot_min_ml_score;
  }
  if (sub->count(cli_opt_name::kPhased) == 0) {
    params->phased = params->config.phased;
  }
}

/**
 * @brief Apply a workflow's static model path when the user did not provide one.
 * @param sub Active workflow subcommand used to detect an explicit model option.
 * @param params CLI parameters to update.
 * @param default_path Built-in model path for the workflow.
 */
void SetNonGermlineDefaultModelPath(const CLI::App* const sub,
                                    const FilterVariantsParamPtr& params,
                                    const std::string_view default_path) {
  if (sub->count(cli_opt_name::kModel) == 0) {
    params->model = default_path;
  }
}

/**
 * @brief Resolve the tumor-normal-wgs model path from sample-type and aligner CLI state.
 * @details When --model is omitted, BWA can auto-resolve the model from the selected sample type. The
 * auto-detect sample type is intentionally left unresolved here so filter_variants can collect metrics
 * during VCF pre-scan and resolve the model afterward. Explicit --model skips this resolution path.
 * @param sub Parsed tumor-normal-wgs subcommand used to detect explicit option values.
 * @param params CLI parameters updated with model_provided, aligner, sample_type, and model state.
 */
void ResolveTumorNormalWgsModelPath(const CLI::App* const sub, const FilterVariantsParamPtr& params) {
  const bool model_provided = OptionProvided(sub, cli_opt_name::kModel);
  const bool aligner_provided = OptionProvided(sub, cli_opt_name::kAligner);
  const bool sample_type_provided = OptionProvided(sub, cli_opt_name::kSampleType);
  params->model_provided = model_provided;
  if (!model_provided) {
    params->model.clear();
  }

  // params->aligner and params->sample_type are shared across all subcommands, and CLI11 fires each
  // subcommand's default callback at definition time — so the last subcommand defined wins the shared
  // default. Reassert this workflow's defaults when the user did not pass the option, making the
  // invoked subcommand authoritative regardless of subcommand definition order.
  if (!aligner_provided) {
    params->aligner = model_provided ? AlignerType::kCustom : AlignerType::kBwa;
  }
  if (!sample_type_provided) {
    params->sample_type = SampleType::kAutoDetect;
  }

  if (params->aligner == AlignerType::kCustom && !model_provided) {
    throw error::Error(fmt::format("{} is required when {} is '{}'",
                                   cli_opt_name::kModel,
                                   cli_opt_name::kAligner,
                                   enum_util::FormatEnumName(AlignerType::kCustom)));
  }

  if (model_provided) {
    return;
  }

  if (params->sample_type == SampleType::kAutoDetect) {
    if (!sample_type_provided && params->aligner != AlignerType::kBwa) {
      params->sample_type = SampleType::kFfpe;
    } else {
      return;
    }
  }

  params->model = ResolveTumorNormalModel(params->sample_type, params->aligner);
}

/**
 * @brief Resolve the tumor-only-wgs model path from sample-type and aligner CLI state.
 * @details The default BWA aligner resolves to a bundled tumor-only model for the selected sample type.
 * Custom aligners require an explicit model path, and explicit --model disables automatic model resolution.
 * @param sub Parsed tumor-only-wgs subcommand used to detect explicit option values.
 * @param params CLI parameters updated with aligner, sample_type, and model state.
 */
void ResolveTumorOnlyWgsModelPath(const CLI::App* const sub, const FilterVariantsParamPtr& params) {
  const bool model_provided = OptionProvided(sub, cli_opt_name::kModel);
  const bool aligner_provided = OptionProvided(sub, cli_opt_name::kAligner);
  const bool sample_type_provided = OptionProvided(sub, cli_opt_name::kSampleType);
  params->model_provided = model_provided;
  if (!model_provided) {
    params->model.clear();
  }

  // params->aligner and params->sample_type are shared across all subcommands, and CLI11 fires each
  // subcommand's default callback at definition time — so the last subcommand defined wins the shared
  // default. Reassert this workflow's defaults when the user did not pass the option, making the
  // invoked subcommand authoritative regardless of subcommand definition order.
  if (!aligner_provided) {
    params->aligner = model_provided ? AlignerType::kCustom : AlignerType::kBwa;
  }
  if (!sample_type_provided) {
    params->sample_type = SampleType::kAutoDetect;
  }

  if (params->aligner == AlignerType::kCustom && !model_provided) {
    throw error::Error(fmt::format("{} is required when {} is '{}'",
                                   cli_opt_name::kModel,
                                   cli_opt_name::kAligner,
                                   enum_util::FormatEnumName(AlignerType::kCustom)));
  }

  if (model_provided) {
    return;
  }

  if (params->sample_type == SampleType::kAutoDetect) {
    return;
  }

  params->model = ResolveTumorOnlyWgsModel(params->sample_type, params->aligner);
}

/**
 * @brief Apply per-model thresholds for tumor-normal-wgs and tumor-only-wgs.
 * @details Runs after top-level config defaults are applied so model-specific thresholds can override
 * config defaults. Explicit CLI values are recorded and preserved. When sample type is auto-detect, this
 * function defers threshold application until sample-type resolution after VCF pre-scan.
 * @param sub Parsed workflow subcommand used to detect explicit threshold and model options.
 * @param params CLI parameters updated with provided flags and resolved thresholds.
 */
void ApplyModelThresholds(const CLI::App* const sub, const FilterVariantsParamPtr& params) {
  params->snv_min_ml_score_provided = OptionCounted(sub, cli_opt_name::kSnvMinMlScore);
  params->indel_min_ml_score_provided = OptionCounted(sub, cli_opt_name::kIndelMinMlScore);
  params->min_af_provided = params->workflow == kTumorOnlyWgs && OptionCounted(sub, cli_opt_name::kMinAf);
  params->min_tumor_af_provided = params->workflow == kTumorNormalWgs && OptionCounted(sub, cli_opt_name::kMinTumorAf);

  const bool auto_resolves_model = (params->workflow == kTumorNormalWgs || params->workflow == kTumorOnlyWgs);
  if (!auto_resolves_model || params->aligner == AlignerType::kCustom) {
    return;
  }
  if (params->sample_type == SampleType::kAutoDetect) {
    return;
  }

  // Only apply when the model was auto-resolved (not explicitly provided via --model)
  if (OptionProvided(sub, cli_opt_name::kModel)) {
    return;
  }

  const auto thresholds = ResolveModelThresholds(params->sample_type, params->aligner, params->config.model_thresholds);
  if (thresholds.has_value()) {
    if (!OptionCounted(sub, cli_opt_name::kSnvMinMlScore) && thresholds->snv_min_ml_score.has_value()) {
      params->snv_min_ml_score = *thresholds->snv_min_ml_score;
    }
    if (!OptionCounted(sub, cli_opt_name::kIndelMinMlScore) && thresholds->indel_min_ml_score.has_value()) {
      params->indel_min_ml_score = *thresholds->indel_min_ml_score;
    }
    if (params->workflow == kTumorOnlyWgs && !params->min_af_provided && thresholds->min_af.has_value()) {
      params->min_af = *thresholds->min_af;
    }
    if (params->workflow == kTumorNormalWgs && !params->min_tumor_af_provided && thresholds->min_tumor_af.has_value()) {
      params->min_tumor_af = *thresholds->min_tumor_af;
    }
  }
}

/**
 * @brief Resolve germline WGS SNV/indel model paths from --library-type, --aligner, and --model-type when needed.
 * @details The default aligner is 'custom', which requires --snv-model and --indel-model to be explicitly provided.
 *          When --aligner is set to 'bwa' or 'giraffe', missing model paths are auto-resolved using --library-type
 *          and --model-type. Any manually supplied model path remains authoritative, allowing partial overrides.
 * @param sub Germline WGS subcommand used to determine which options were explicitly supplied.
 * @param params CLI parameters to update with resolved aligner and model paths.
 * @throws CLI::ValidationError for an unsupported library-type/aligner combination.
 * @throws xoos::error::Error when automatic model resolution lacks a supported aligner.
 */
void ResolveGermlineModelPaths(const CLI::App* const sub, const FilterVariantsParamPtr& params) {
  const bool snv_model_provided = sub->count(cli_opt_name::kSnvModel) > 0;
  const bool indel_model_provided = sub->count(cli_opt_name::kIndelModel) > 0;

  // Explicit model pairs bypass automatic model selection because the user is responsible for model compatibility.
  if (snv_model_provided && indel_model_provided) {
    return;
  }

  if (params->aligner == AlignerType::kCustom) {
    throw error::Error(fmt::format("{} and {} are required when {} is '{}'",
                                   cli_opt_name::kSnvModel,
                                   cli_opt_name::kIndelModel,
                                   cli_opt_name::kAligner,
                                   enum_util::FormatEnumName(AlignerType::kCustom)));
  }

  const auto [snv_path, indel_path] =
      ResolveGermlineWgsModels(params->library_type, params->aligner, params->germline_model_type);
  if (!snv_model_provided) {
    params->snv_model = snv_path;
  }
  if (!indel_model_provided) {
    params->indel_model = indel_path;
  }
}

/**
 * @brief Helper function to set default model paths based on workflow if not provided via CLI.
 * @param sub Subcommand CLI application pointer to check which options were specified via CLI
 * @param params CLI parameters shared pointer to apply defaults to
 */
void SetDefaultModelPaths(const CLI::App* const sub, const FilterVariantsParamPtr& params) {
  switch (params->workflow) {
    case kGermlineWgs:
      ResolveGermlineModelPaths(sub, params);
      break;
    case kTumorOnlyTe:
      SetNonGermlineDefaultModelPath(sub, params, kDefaultTumorOnlyTeModelPath);
      break;
    case kTumorNormalWgs:
      ResolveTumorNormalWgsModelPath(sub, params);
      break;
    case kTumorOnlyWgs:
      ResolveTumorOnlyWgsModelPath(sub, params);
      break;
    case kGermlineMitochondria:
      // Mitochondria workflow doesn't use ML models
      break;
    default:
      SetNonGermlineDefaultModelPath(sub, params, "");
      break;
  }
}

/**
 * @brief Synchronize automatically resolved model settings with the command line rendered for provenance.
 * @details This must run after ApplyConfig, when model paths and thresholds are final, and before PreCallback captures
 * the command line. Otherwise the output VCF would record empty models or generic thresholds instead of used values.
 * @param sub Active workflow subcommand whose rendered defaults are updated.
 * @param params Resolved filtering parameters to reflect in the rendered command line.
 */
void SyncResolvedModelOptions(CLI::App* const sub, const FilterVariantsParamPtr& params) {
  if (params->workflow == kGermlineWgs) {
    SetRenderedDefault(sub, cli_opt_name::kAligner, enum_util::FormatEnumName(params->aligner));
    SetRenderedDefault(sub, cli_opt_name::kLibraryType, enum_util::FormatEnumName(params->library_type));
    SetRenderedDefault(sub, cli_opt_name::kModelType, enum_util::FormatEnumName(params->germline_model_type));
    SetRenderedDefault(sub, cli_opt_name::kSnvModel, params->snv_model.string());
    SetRenderedDefault(sub, cli_opt_name::kIndelModel, params->indel_model.string());
    return;
  }

  const bool auto_resolves_model = (params->workflow == kTumorNormalWgs || params->workflow == kTumorOnlyWgs);
  if (params->workflow == kTumorNormalWgs) {
    SetRenderedDefault(sub, cli_opt_name::kMinTumorAf, fmt::format("{}", params->min_tumor_af));
  } else if (params->workflow == kTumorOnlyWgs) {
    SetRenderedDefault(sub, cli_opt_name::kMinAf, fmt::format("{}", params->min_af));
  }
  if (!auto_resolves_model) {
    return;
  }
  if (OptionProvided(sub, cli_opt_name::kModel)) {
    return;
  }
  if (params->aligner == AlignerType::kCustom) {
    return;
  }
  if (params->sample_type == SampleType::kAutoDetect) {
    return;
  }
  SetRenderedDefault(sub, cli_opt_name::kAligner, enum_util::FormatEnumName(params->aligner));
  SetRenderedDefault(sub, cli_opt_name::kSampleType, enum_util::FormatEnumName(params->sample_type));
  SetRenderedDefault(sub, cli_opt_name::kModel, params->model.string());
  SetRenderedDefault(sub, cli_opt_name::kSnvMinMlScore, fmt::format("{}", params->snv_min_ml_score));
  SetRenderedDefault(sub, cli_opt_name::kIndelMinMlScore, fmt::format("{}", params->indel_min_ml_score));
}

/**
 * @brief Helper function to apply CLI defaults and workflow config presets to CLI parameters if the corresponding CLI
 * options were not specified.
 * @param sub Subcommand CLI application pointer to check which options were specified via CLI
 * @param params CLI parameters shared pointer to apply defaults to
 */
void ApplyConfigToMitochondriaUniqueOptions(const CLI::App* const sub, const FilterVariantsParamPtr& params) {
  if (sub->count(cli_opt_name::kMinAf) == 0) {
    params->min_af = params->config.min_af;
  }
  if (sub->count(cli_opt_name::kMinDp) == 0) {
    params->min_dp = params->config.min_dp;
  }
  if (sub->count(cli_opt_name::kMinHomAltAf) == 0) {
    params->min_hom_alt_af = params->config.min_hom_alt_af;
  }
}

/**
 * @brief Apply the selected profile and derived defaults to filtering parameters.
 * @details Explicit CLI values take precedence. Derived output paths, model selection, feature extraction defaults,
 * and model-specific thresholds are applied in dependency order.
 * @param sub Active workflow subcommand used to identify explicitly provided options.
 * @param params Filtering parameters to update.
 * @throws CLI::ValidationError if a derived output path already exists.
 */
void ApplyConfig(const CLI::App* const sub, const FilterVariantsParamPtr& params) {
  if (sub->count(cli_opt_name::kVcfOutput) == 0) {
    // update VCF output path to be relative to output directory
    params->vcf_output = fs::path(params->output_dir) / params->vcf_output;
    if (fs::exists(params->vcf_output)) {
      throw CLI::ValidationError(fmt::format("Output VCF file '{}' already exists", params->vcf_output.string()));
    }
  }

  if (sub->count(cli_opt_name::kOutputMetrics) == 0) {
    // update metrics output path to be relative to output directory
    params->metrics_output = fs::path(params->output_dir) / params->metrics_output;
  }
  if (fs::exists(params->metrics_output)) {
    throw CLI::ValidationError(fmt::format("Output metrics file '{}' already exists", params->metrics_output.string()));
  }

  // resolve SHAP value TSV paths relative to output directory (order-independent)
  const auto resolve_shap_path = [&params](std::optional<fs::path>& path) {
    if (path.has_value()) {
      path = params->output_dir / path.value();
      if (fs::exists(path.value())) {
        throw CLI::ValidationError(fmt::format("SHAP value output file '{}' already exists", path->string()));
      }
    }
  };
  resolve_shap_path(params->shap_value_tsv);
  resolve_shap_path(params->snv_shap_value_tsv);
  resolve_shap_path(params->indel_shap_value_tsv);

  // set default model paths if not provided via CLI
  SetDefaultModelPaths(sub, params);

  // apply config to compute-bam-features options
  compute_bam_features::ApplyConfig(sub, params);

  // The germline-mitochondria subcommand does not define --normalize-features (it filters on AF/depth only), so guard
  // the lookup to avoid querying an option that does not exist for that workflow.
  if (params->config.workflow != kGermlineMitochondria && sub->count(cli_opt_name::kNormalizeFeatures) == 0) {
    params->normalize_features = params->config.normalize_features;
  }

  // apply config to subcommand unique options
  ApplyConfigToSharedMlScoreOptions(sub, params);
  switch (params->config.workflow) {
    case kTumorNormalWgs:
      ApplyConfigToTumorNormalWgsUniqueOptions(sub, params);
      break;
    case kTumorOnlyWgs:
      ApplyConfigToTumorOnlyWgsUniqueOptions(sub, params);
      break;
    case kTumorOnlyTe:
      ApplyConfigToTumorOnlyTeUniqueOptions(sub, params);
      break;
    case kGermlineMitochondria:
      ApplyConfigToMitochondriaUniqueOptions(sub, params);
      break;
    default:
      break;
  }

  // apply per-model thresholds after top-level config defaults so they take precedence
  ApplyModelThresholds(sub, params);
}
}  // namespace

/**
 * @brief Preprocess the CLI options and validate the parameters.
 * @param app CLI application pointer
 * @param params CLI parameters to be validated
 * @throws CLI::ValidationError if the parameters are invalid
 */
void filter_variants::PreCallback(const CLI::App* const app, const FilterVariantsParamPtr& params) {
  // Check which subcommand was used, set workflow and config accordingly, and apply config defaults as needed
  for (const Workflow workflow : kSupportedWorkflows) {
    const std::string name = enum_util::FormatEnumName(workflow);
    if (app->got_subcommand(name)) {
      params->workflow = workflow;
      if (workflow == kGermlineWgs) {
        params->library_type = ResolveLibraryType(params->library_type, params->bam_files);
        params->config =
            JsonToConfig(params->config_file.value_or(fs::path{}), GermlineWgsConfigProfile(params->library_type));
      } else {
        params->config = JsonToConfig(params->config_file.value_or(fs::path{}), params->workflow);
      }
      // get_subcommand on a const app returns a non-const App*, allowing the resolved selections to be
      // written back into the subcommand's options below.
      CLI::App* const sub = app->get_subcommand(name);
      ApplyConfig(sub, params);
      // Reflect auto-resolved model and thresholds in the options so the captured command line records them.
      SyncResolvedModelOptions(sub, params);
      if (workflow == kGermlineWgs) {
        SyncGermlineWgsRenderedDefaults(sub, params);
        SetRenderedDefault(
            sub, cli_opt_name::kNormalizeFeatures, enum_util::FormatEnumName(params->normalize_features));
      }
      break;
    }
  }

  // Capture the command line after config and model/threshold resolution so the provenance recorded in the
  // output VCF and logs reflects the values actually used for filtering, not just the raw user input.
  params->command_line = cli::GetCommandLineInfo(app);
  if (params->workflow == kGermlineWgs) {
    ReplaceOrAppendCommandLineOption(params->command_line->command_line,
                                     cli_opt_name::kLibraryType,
                                     enum_util::FormatEnumName(params->library_type));
  }

  if (params->decode_yc != YcDecodeMethod::kNone && !IsDuplexProtocol(params->sequencing_protocol)) {
    throw CLI::ValidationError("Duplex sequencing protocol must be used when decoding YC tags");
  }
}

}  // namespace xoos::svc
