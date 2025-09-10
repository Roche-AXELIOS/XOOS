#include "compute-bam-features-cli.h"

#include <xoos/cli/enum-option-util.h>
#include <xoos/cli/thread-count-option-util.h>

#include "util/cli-util.h"

namespace xoos::svc {

using cli::AddOptionalEnumOption;
using cli::AddThreadCountOption;

static SVCConfig GetConfig(const ComputeBamFeaturesCliParamsPtr& params, const fs::path& config_json) {
  SVCConfig model_config = JsonToConfig(config_json, params->workflow);
  // Load the default values from the config. If the user overwrites a value it is done below
  params->min_mapq = model_config.min_mapq;
  params->min_bq = model_config.min_bq;
  params->min_allowed_distance_from_end = model_config.min_allowed_distance_from_end;
  params->min_family_size = model_config.min_family_size;
  params->filter_homopolymer = model_config.filter_homopolymer;
  params->min_homopolymer_length = model_config.min_homopolymer_length;
  params->duplex = model_config.duplex;
  params->decode_yc = model_config.decode_yc;
  return model_config;
}

void DefineOptions(cli::AppPtr app, const ComputeBamFeaturesCliParamsPtr& params) {
  app->add_option("--bam-input", params->bam_input, "Input BAM file(s) to be analyzed")
      ->required()
      ->check(kCliIndexedBamFile);
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
      ->default_val(kDefaultBamFeaturesFileName);
  app->add_option_function<fs::path>(
         "--target-regions",
         [&params](const fs::path& value) { params->bed_regions = GetBedRegions(value); },
         "Path to a BED file of target regions")
      ->check(kCliBedFile);
  app->add_option(
         "--min-mapq", params->min_mapq, "Minimum alignment mapping quality required to support a variant (inclusive)")
      ->check(kCliRangeMapq);
  app->add_option(
         "--min-bq", params->min_bq, "Minimum alignment base quality required to support a variant (inclusive)")
      ->check(kCliRangeBaseq);
  app->add_option("--min-dist",
                  params->min_allowed_distance_from_end,
                  "Minimum distance of variant from fragment alignment end (inclusive)")
      ->check(CLI::NonNegativeNumber);
  app->add_option("--min-family-size",
                  params->min_family_size,
                  "Minimum cluster size required for an alignment to support a variant (inclusive)")
      ->check(CLI::NonNegativeNumber);
  app->add_option("--max-variants-per-read",
                  params->max_read_variant_count,
                  "Max number of variants allowed per read (inclusive); `0` can also turn off this option")
      ->check(CLI::NonNegativeNumber);
  app->add_option("--max-variants-per-read-normalized",
                  params->max_read_variant_count_normalized,
                  "Max number of variants allowed per read, normalized by alignment length (inclusive); `0` can also "
                  "turn off this option")
      ->check(kCliRangeFraction);
  app->add_option(
         "--skip-variants-vcf",
         params->skip_variants_vcf,
         "VCF containing variants not counted by `--max-variants-per-read` or --max-variants-per-read-normalized`")
      ->check(kCliIndexedVcfFile);
  AddThreadCountOption(app, "--threads", params->threads);
  app->add_flag("--filter-homopolymer,!--no-filter-homopolymer",
                params->filter_homopolymer,
                "skip variant adjacent to homopolymer that spans beyond a read's 3' end");
  app->add_option(
         "--min-homopolymer-length", params->min_homopolymer_length, "Minimum length of homopolymer in reference")
      ->check(CLI::NonNegativeNumber);
  app->add_flag("--duplex,!--no-duplex", params->duplex, "input BAM contains duplex data");
#ifdef SOMATIC_ENABLED
  app->add_option("--tumor-read-group", params->tumor_read_group, "tumor sample's read group name");
#endif  // SOMATIC_ENABLED
  app->add_option("--max-region-size-per-thread", params->max_region_size_per_thread, "Maximum region size per thread")
      ->default_val(kDefaultMaxBamRegionSizePerThread)
      ->check(CLI::PositiveNumber);
  app->add_flag("--decode-yc,!--no-decode-yc", params->decode_yc, "decode YC tags within input BAM file(s)")
      ->needs("--duplex");
  cli::AddEnumOption(app,
                     "--min-base-type",
                     params->min_base_type,
                     "minimum base type in duplex reads for variant support",
                     yc_decode::BaseType::kSimplex)
      ->needs("--decode-yc");
  AddWarnAsErrorOption(app);
}

}  // namespace xoos::svc
