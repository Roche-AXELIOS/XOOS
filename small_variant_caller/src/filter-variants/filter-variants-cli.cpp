#include "filter-variants-cli.h"

#ifdef SOMATIC_ENABLE
#include <algorithm>
#endif  // SOMATIC_ENABLE

#include <xoos/cli/enum-option-util.h>
#include <xoos/cli/file-option-util.h>
#include <xoos/cli/thread-count-option-util.h>

#include "core/command-line-info.h"
#include "util/cli-util.h"

namespace xoos::svc {

using cli::AddOptionalEnumOption;
using cli::AddOutputFileOption;
using cli::AddThreadCountOption;

/**
 * @brief Get the configuration based on the CLI parameters and JSON configuration file.
 * @param params CLI parameters
 * @param config_json Path to the JSON configuration file
 * @return SVCConfig object containing the configuration
 */
static SVCConfig GetConfig(const FilterVariantsParamPtr& params, const fs::path& config_json) {
  SVCConfig model_config = JsonToConfig(config_json, params->workflow);
  // Load the default values from the config. If the user overwrites a value it is done below
  params->min_mapq = model_config.min_mapq;
  params->min_bq = model_config.min_bq;
  params->min_allowed_distance_from_end = model_config.min_allowed_distance_from_end;
  params->min_family_size = model_config.min_family_size;
  params->filter_homopolymer = model_config.filter_homopolymer;
  params->min_homopolymer_length = model_config.min_homopolymer_length;
  params->duplex = model_config.duplex;
  params->min_alt_counts = model_config.min_alt_counts;
  params->min_allele_freq_threshold = model_config.min_ctdna_allele_freq_threshold;
  params->min_phased_allele_freq = model_config.min_phased_allele_freq;
  params->max_phased_allele_freq = model_config.max_phased_allele_freq;
  params->weighted_counts_threshold = model_config.weighted_counts_threshold;
  params->hotspot_weighted_counts_threshold = model_config.hotspot_weighted_counts_threshold;
  params->ml_threshold = model_config.ml_threshold;
  params->hotspot_ml_threshold = model_config.hotspot_ml_threshold;
  params->dynamic_thresholds = model_config.dynamic_thresholds;
  params->phased = model_config.phased;
  params->normalize_features = model_config.normalize_features;
  params->somatic_tn_snv_ml_threshold = model_config.somatic_tn_snv_ml_threshold;
  params->somatic_tn_indel_ml_threshold = model_config.somatic_tn_indel_ml_threshold;
  params->tumor_support_threshold = model_config.tumor_support_threshold;
  params->decode_yc = model_config.decode_yc;
  return model_config;
}

const vec<fs::path> kDefaultModelPaths{fs::path("/resources/model-germline-snv.txt.gz"),
                                       fs::path("/resources/model-germline-indel.txt.gz")};
const std::string kDefaultModelPathsStr{"/resources/model-germline-snv.txt.gz /resources/model-germline-indel.txt.gz"};
const std::string kBamFeatureExtractionOptions{"BAM feature extraction options"};
const std::string kVcfFeatureExtractionOptions{"VCF feature extraction options"};

void DefineOptionsFilterVariants(cli::AppPtr app, const FilterVariantsParamPtr& params) {
  app->add_option("--bam-input", params->bam_files, "Input BAM file to be analyzed")
      ->required()
      ->check(kCliIndexedBamFile);
  app->add_option("--vcf-input", params->vcf_file, "Input VCF file, produced by GATK Mutect2/HaplotypeCaller (1-based)")
      ->required()
      ->check(kCliIndexedVcfFile);
  app->add_option("--genome", params->genome, "Path to reference genome (indexed FASTA)")
      ->required()
      ->check(kCliIndexedFastaFile);
  app->add_option("--model", params->model, "Path to LightGBM model file(s)")
      ->default_val(kDefaultModelPaths)
      ->default_str(kDefaultModelPathsStr)
      ->check(kCliNonEmptyFile);
  AddOptionalEnumOption(app, "--workflow", params->workflow, "compute features for the designated workflow")
      ->required();
  // The `--workflow` option must be defined before the `--config` option.
  // Otherwise, `params->workflow` will always have the default value (somatic) within the `GetConfig` function.
  app->add_option_function<fs::path>(
         "--config",
         [&params](const fs::path& value) { params->config = GetConfig(params, value); },
         "Path to config JSON file")
      ->force_callback();
  app->add_option("--output-dir", params->output_dir, "Output directory")->default_val(".");
  AddOutputFileOption(app,
                      "--vcf-output",
                      params->vcf_output,
                      "VCF output location, relative to out-dir",
                      "output.vcf.gz",
                      params->output_dir);
  app->add_option("--target-regions", params->bed_file, "Path to a BED file of target regions")->check(kCliBedFile);
  app->add_option(
         "--interest-regions", params->interest_bed_file, "Path to a BED file of regions for VCF feature `at_interest`")
      ->check(kCliBedFile)
      ->group(kVcfFeatureExtractionOptions);
  app->add_option("--pop-af-vcf",
                  params->pop_af_vcf,
                  "GATK gnomAD population allele frequency VCF (i.e. af-only-gnomad.hg38.vcf.gz)")
      ->check(kCliIndexedVcfFile)
      ->group(kVcfFeatureExtractionOptions);
  app->add_option(
         "--min-mapq", params->min_mapq, "Minimum alignment mapping quality required to support a variant (inclusive)")
      ->check(kCliRangeMapq)
      ->group(kBamFeatureExtractionOptions);
  app->add_option(
         "--min-bq", params->min_bq, "Minimum alignment base quality required to support a variant (inclusive)")
      ->check(kCliRangeBaseq)
      ->group(kBamFeatureExtractionOptions);
  app->add_option("--min-dist",
                  params->min_allowed_distance_from_end,
                  "Minimum distance of variant from fragment alignment end (inclusive)")
      ->check(CLI::NonNegativeNumber)
      ->group(kBamFeatureExtractionOptions);
  app->add_option("--min-family-size",
                  params->min_family_size,
                  "Minimum cluster size required for an alignment to support a variant (inclusive)")
      ->check(CLI::NonNegativeNumber)
      ->group(kBamFeatureExtractionOptions);
  app->add_option("--max-variants-per-read",
                  params->max_read_variant_count,
                  "Max number of variants allowed per read (inclusive); `0` can also turn off this option")
      ->check(CLI::NonNegativeNumber)
      ->group(kBamFeatureExtractionOptions);
  app->add_option("--max-variants-per-read-normalized",
                  params->max_read_variant_count_normalized,
                  "Max number of variants allowed per read, normalized by alignment length (inclusive); `0` can also "
                  "turn off this option")
      ->check(kCliRangeFraction)
      ->group(kBamFeatureExtractionOptions);
  app->add_option(
         "--skip-variants-vcf",
         params->skip_variants_vcf,
         "VCF containing variants not counted by `--max-variants-per-read` or --max-variants-per-read-normalized`")
      ->check(kCliIndexedVcfFile)
      ->group(kBamFeatureExtractionOptions);
#ifdef SOMATIC_ENABLED
  app->add_option("--block-list", params->block_list, "Block list file")->check(kCliNonEmptyFile);
  app->add_option("--hotspot-list", params->hotspot_list, "A VCF file containing a list of hotspot variants")
      ->check(kCliNonEmptyFile);
  app->add_option("--forcecall-list",
                  params->forcecall_list,
                  "A 0-based BED file with forced call positions to be included in the output")
      ->check(kCliNonEmptyFile);
  app->add_option(
         "--min-alt-counts", params->min_alt_counts, "Minimum number of alt counts for retaining a variant (inclusive)")
      ->check(CLI::NonNegativeNumber);
  app->add_option(
         "--af-threshold", params->min_allele_freq_threshold, "Minimum allowed variant allele frequency (inclusive)")
      ->check(kCliRangeFraction);
  app->add_option("--min-phased-af",
                  params->min_phased_allele_freq,
                  "Minimum allowed phased allele frequency for a variant (inclusive)")
      ->check(kCliRangeFraction);
  app->add_option("--max-phased-af",
                  params->max_phased_allele_freq,
                  "Maximum allowed phased allele frequency for a variant (inclusive)")
      ->check(kCliRangeFraction);
  app->add_flag("--phased", params->phased, "Call phased variants in VCF");
  app->add_flag(
      "--dynamic-wc-thresholds", params->dynamic_thresholds, "Determine weighted count thresholds dynamically.");
  app->add_option("--wc-threshold",
                  params->weighted_counts_threshold,
                  "Threshold for weighted counts (inclusive). Variants must score at least this much to be included")
      ->check(CLI::NonNegativeNumber);
  app->add_option(
         "--hotspot-wc-threshold",
         params->hotspot_weighted_counts_threshold,
         "Threshold for hotspot weighted counts (inclusive). Variants must score at least this much to be included")
      ->check(CLI::NonNegativeNumber);
  app->add_option("--ml-threshold",
                  params->ml_threshold,
                  "Threshold for ML score (inclusive). Variants must score at least this much to be included")
      ->check(kCliRangeFraction);
  app->add_option("--somatic-tn-snv-ml-threshold",
                  params->somatic_tn_snv_ml_threshold,
                  "ML score threshold for Somatic TN for SNVs (inclusive). Variants must score at least this much to "
                  "be included")
      ->check(kCliRangeFraction);
  app->add_option("--somatic-tn-indel-ml-threshold",
                  params->somatic_tn_indel_ml_threshold,
                  "ML score threshold for Somatic TN for InDels (inclusive). Variants must score at least this much to "
                  "be included")
      ->check(kCliRangeFraction);
  app->add_option("--tumor-support-threshold",
                  params->tumor_support_threshold,
                  "Tumor support threshold for Somatic TN (inclusive). Variants must have at least this much support "
                  "to be included")
      ->check(CLI::NonNegativeNumber);
  app->add_option(
         "--hotspot-ml-threshold",
         params->hotspot_ml_threshold,
         "Threshold for ML score in hotspots (inclusive). Variants must score at least this much to be included")
      ->check(kCliRangeFraction);
  app->add_option_function<std::string>(
         "--sample-type",
         [&params](const std::string& value) {
           if (value == "FFPE") {
             params->min_allele_freq_threshold = std::max(params->min_allele_freq_threshold, kDefaultMinFFPEAf);
           }
         },
         "Sample type, either FFPE or ctDNA.")
      ->default_val(kDefaultSampleType)
      ->force_callback();
#endif  // SOMATIC_ENABLED
  AddThreadCountOption(app, "--threads", params->threads);
  app->add_flag("--filter-homopolymer,!--no-filter-homopolymer",
                params->filter_homopolymer,
                "skip variant adjacent to homopolymer that spans beyond a read's 3' end")
      ->group(kBamFeatureExtractionOptions);
  app->add_option("--min-homopolymer-length",
                  params->min_homopolymer_length,
                  "Minimum length of homopolymer in reference. (inclusive)")
      ->check(CLI::NonNegativeNumber)
      ->group(kBamFeatureExtractionOptions);
  app->add_flag("--duplex,!--no-duplex", params->duplex, "input BAM contains duplex data")
      ->group(kBamFeatureExtractionOptions);
  app->add_option("--max-bam-region-size-per-thread",
                  params->max_bam_region_size_per_thread,
                  "Maximum BAM region size per thread during feature extraction (inclusive)")
      ->default_val(kDefaultMaxBamRegionSizePerThread)
      ->check(CLI::PositiveNumber);
  app->add_option("--max-vcf-region-size-per-thread",
                  params->max_vcf_region_size_per_thread,
                  "Maximum VCF region size per thread during variant filtration (inclusive)")
      ->default_val(kDefaultMaxVcfRegionSizePerThread)
      ->check(CLI::PositiveNumber);
  // First, `run_callback_default()` ensures the default is computed in `transform()`,
  // then `force_callback()` ensures `transform()` is run if the number of threads is specified.
  app->add_option("--filter-variant-window-size",
                  params->filter_variant_window_size,
                  "Window size for filtering variants, the default is 0 which means 8 * --threads")
      ->transform([&threads = params->threads](const std::string& value) -> std::string {
        auto requested_window_size = std::stoul(value);
        return std::to_string(requested_window_size == 0 ? threads * 8 : requested_window_size);
      })
      ->force_callback()
      ->run_callback_for_default()
      ->default_val(0)
      ->check(CLI::PositiveNumber);
  app->add_flag("--normalize,!--no-normalize", params->normalize_features, "Normalize read count features");
  app->add_option("--sd-chr-name", params->sd_chr_name, "autosome name for sex determination")
      ->default_val(kDefaultChr1Name);
  app->add_option("--par-bed-x", params->par_bed_x, "Path to chromosome X pseudoautosommal region BED file")
      ->check(kCliNonEmptyFile);
  app->add_option("--par-bed-y", params->par_bed_y, "Path to chromosome Y pseudoautosommal region BED file")
      ->check(kCliNonEmptyFile);
  app->add_flag("--decode-yc,!--no-decode-yc", params->decode_yc, "decode YC tags within input BAM file(s)")
      ->needs("--duplex")
      ->group(kBamFeatureExtractionOptions);
  cli::AddEnumOption(app,
                     "--min-base-type",
                     params->min_base_type,
                     "minimum base type in duplex reads for variant support",
                     yc_decode::BaseType::kSimplex)
      ->needs("--decode-yc")
      ->group(kBamFeatureExtractionOptions);
  AddWarnAsErrorOption(app);
}

/**
 * @brief Preprocess the CLI options and validate the parameters.
 * @param app CLI application pointer
 * @param params FilterVariantsParamPtr containing the parameters
 */
void FilterVariantsPreCallback(cli::ConstAppPtr app, const FilterVariantsParamPtr& params) {
  params->command_line = GetCommandLineInfo(app);
  using enum Workflow;
  switch (params->workflow) {
    case kGermlineMultiSample:
    case kGermline: {
      break;
    }
    default: {
      throw CLI::ValidationError("Workflow not Supported for Filtering");
    }
  }
}

}  // namespace xoos::svc
