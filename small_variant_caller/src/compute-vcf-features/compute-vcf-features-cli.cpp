#include "compute-vcf-features-cli.h"

#include <xoos/cli/enum-option-util.h>
#include <xoos/cli/thread-count-option-util.h>

#include "util/cli-util.h"

namespace xoos::svc {

using cli::AddOptionalEnumOption;
using cli::AddThreadCountOption;

static SVCConfig GetConfig(const ComputeVcfFeaturesParamPtr& params, const fs::path& config_json) {
  SVCConfig model_config = JsonToConfig(config_json, params->workflow);
  // CLI parameters to be added in the future can be updated here with the config settings.
  return model_config;
}

void DefineOptionsComputeVcfFeatures(cli::AppPtr app, const ComputeVcfFeaturesParamPtr& params) {
  app->add_option("--vcf-input", params->vcf_file, "Path to input VCF file")->required()->check(kCliIndexedVcfFile);
  app->add_option("--genome", params->genome, "Path to reference genome (indexed FASTA)")
      ->required()
      ->check(kCliIndexedFastaFile);
  AddOptionalEnumOption(app, "--workflow", params->workflow, "Compute features for the designated workflow")
      ->required();
  // The `--workflow` option must be defined before the `--config` option.
  // Otherwise, `params->workflow` will always have the default value (somatic) within the `GetConfig` function.
  app->add_option_function<fs::path>(
         "--config",
         [&params](const fs::path& value) { params->config = GetConfig(params, value); },
         "Path to config JSON file")
      ->force_callback();
  app->add_option("--output-file", params->output_file, "Output features file name")
      ->default_val(kDefaultVcfFeaturesFileName);
  app->add_option_function<fs::path>(
         "--target-regions",
         [&params](const fs::path& value) { params->target_regions = GetChromIntervalMap(value); },
         "Path to a BED file of target regions")
      ->check(kCliBedFile);
  app->add_option_function<fs::path>(
         "--interest-regions",
         [&params](const fs::path& value) { params->interest_regions = GetChromIntervalMap(value); },
         "Path to a BED file of regions for VCF feature `at_interest_region`")
      ->check(kCliBedFile);
  AddThreadCountOption(app, "--threads", params->threads);
  app->add_option("--pop-af-vcf",
                  params->pop_af_vcf,
                  "Path to GATK gnomAD population allele frequency VCF (i.e. af-only-gnomad.hg38.vcf.gz)")
      ->check(kCliIndexedVcfFile);
  app->add_option("--output-bed", params->output_bed, "Path to output BED file for extracted features");
  app->add_option("--left-pad", params->left_pad, "left-padding to variant start position for output BED file")
      ->default_val(kDefaultLeftPad)
      ->check(CLI::NonNegativeNumber);
  app->add_option("--right-pad", params->right_pad, "right-padding to variant end position for output BED file")
      ->default_val(kDefaultRightPad)
      ->check(CLI::NonNegativeNumber);
  app->add_option("--collapse-distance",
                  params->collapse_dist,
                  "distance to collapse nearby variants after padding is applied for output BED file")
      ->default_val(kDefaultCollapsableDist)
      ->check(CLI::NonNegativeNumber);
  AddWarnAsErrorOption(app);
}

}  // namespace xoos::svc
