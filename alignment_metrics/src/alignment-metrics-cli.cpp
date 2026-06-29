#include "alignment-metrics-cli.h"

#include <string>

#include <fmt/format.h>

#include <xoos/cli/bam-option-util.h>
#include <xoos/cli/cli.h>
#include <xoos/cli/enum-option-util.h>
#include <xoos/cli/fasta-option-util.h>
#include <xoos/cli/thread-count-option-util.h>
#include <xoos/enum/enum-util.h>
#include <xoos/types/int.h>

#include "alignment-metrics-options.h"
#include "run-alignment-metrics.h"

namespace xoos::alignment_metrics {

using cli::AddBamFileOption;
using cli::AddFastaFileOption;
using cli::AddThreadCountOption;

void AddInputOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts, const bool requires_reference) {
  AddBamFileOption(app,
                   "-b,--bam-input",
                   opts->bam_input,
                   opts->bam_index,
                   "Path to the BAM file which the metrics are to be calculated. The BAM file must be "
                   "coordinate-sorted and indexed.")
      ->required()
      ->group(kOptGroupNameInputOptions);
  app->add_option(
         "-i,--bed-input",
         opts->bed_input,
         "Path to a BED file which specifies regions of interest for the metrics calculation. If specified, only reads "
         "overlapping with a region of interest and positions contained within a region of interest will contribute to "
         "the metrics. If not specified, all contigs present in the BAM header will contribute to the metrics.")
      ->group(kOptGroupNameInputOptions)
      ->check(CLI::ExistingFile);
  AddFastaFileOption(
      app,
      "-r,--reference",
      opts->reference_input,
      "Path to the reference FASTA file used to compare the reads against. The FASTA file must be indexed with a "
      "`.fai` file. This is required for calculating accuracy metrics as well as homopolymer coverage metrics.")
      ->required(requires_reference)
      ->group(kOptGroupNameInputOptions);
}

void AddOutputOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts) {
  app->add_option(
         "-o,--output-dir", opts->out_dir, "Path to the output directory where the metrics will be written to.")
      ->default_val(".")
      ->group(kOptGroupNameOutputOptions);
}

void AddReadFilterOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts) {
  app->add_option("-e,--exclude-flags",
                  opts->exclude_flags,
                  "Exclude alignments from coverage metrics, accuracy metrics, and post-filter read metrics if any "
                  "bits in their FLAG field match the specified decimal integer. Refer to "
                  "https://broadinstitute.github.io/picard/explain-flags.html for flag generation. Default is "
                  "supplementary, secondary, unmapped, and duplicate (3332)")
      ->default_val(kDefaultExcludeFlag)
      ->group(kOptGroupNameReadFilterOptions);
  app->add_option("-Q,--min-mapq",
                  opts->min_mapq,
                  "Minimum mapping quality. Reads with a mapping quality less than the threshold are not considered "
                  "for coverage metrics, accuracy metrics, and post-filter read metrics.")
      ->default_val(std::to_string(kDefaultMinMapq))
      ->check(CLI::Range(kDefaultMinMapq, kMaxBaseQuality))
      ->group(kOptGroupNameReadFilterOptions);
}

void AddReadTrimmingOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts) {
  app->add_option(
         "--trim-leading-bases", opts->trim_leading_bases, "Number of bases to trim from the start of each read")
      ->default_val(0)
      ->check(CLI::Range(0u, kMaxTrimLeadingBases))
      ->group(kOptGroupNameReadTrimmingOptions);
  app->add_option(
         "--trim-trailing-bases", opts->trim_trailing_bases, "Number of bases to trim from the end of each read")
      ->default_val(0)
      ->check(CLI::Range(0u, kMaxTrimTrailingBases))
      ->group(kOptGroupNameReadTrimmingOptions);
}

void AddBaseFilterOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts) {
  app->add_option("-q,--min-baseq",
                  opts->min_baseq,
                  "Minimum base quality. Bases with a base quality less than the threshold are not considered for "
                  "accuracy metrics.")
      ->default_val(std::to_string(kDefaultMinBaseq))
      ->check(CLI::Range(0u, 255u))
      ->group(kOptGroupNameBaseFilterOptions);
  app->add_flag("--disable-base-type-decoding",
                opts->disable_base_type_decoding,
                "Skip decoding base types from the YC tag in the BAM file. This can speed up processing for simplex "
                "datasets where base type decoding is not necessary.")
      ->group(kOptGroupNameBaseFilterOptions);
  cli::AddOptionalEnumOption(
      app,
      "--min-base-type",
      opts->min_base_type,
      "Minimum duplex base type. For duplex datasets, bases with a base type worse than the threshold are not "
      "considered for accuracy and post-filter coverage metrics. For simplex datasets, it should always be set to "
      "`simplex`.")
      ->default_val(enum_util::FormatEnumName(kDefaultMinBaseType))
      ->group(kOptGroupNameBaseFilterOptions)
      ->check([&opts](const std::string& input) -> std::string {
        const auto base_type = cli::ParseEnumNameOrThrow<yc_decode::BaseType>("min-base-type", input);
        if (opts->disable_base_type_decoding && base_type == yc_decode::BaseType::kConcordant) {
          throw CLI::ValidationError(
              "When --disable-base-type-decoding is enabled, --min-base-type must not be set to `concordant` as no "
              "base type "
              "information will be available and all bases will be treated as simplex. Setting --min-base-type to "
              "`concordant` in this case will cause all bases to be filtered out.");
        }
        return "";
      });
}

void AddHpMetricsOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts) {
  const auto hp_metrics_flag =
      app->add_flag("--enable-hp-metrics", opts->calculate_hp_metrics, "Enable calculation of homopolymer metrics.")
          ->group(kOptGroupNameHpMetricsOptions)
          ->default_val(false)
          ->needs("--reference");

  // Need to cast to u16 for range check because CLI11 treats u8 as
  // unsigned char and it won't render correctly
  app->add_option(
         "--min-hp-length", opts->min_hp_length, "Minimum homopolymer length to consider for homopolymer metrics.")
      ->default_val(static_cast<u16>(kHpMinLength))
      ->check(CLI::Range(static_cast<u16>(kHpMinLength), static_cast<u16>(kHpMaxLength)))
      ->group(kOptGroupNameHpMetricsOptions)
      ->needs(hp_metrics_flag);

  app->add_option(
         "--max-hp-length", opts->max_hp_length, "Maximum homopolymer length to consider for homopolymer metrics.")
      ->default_val(static_cast<u16>(kHpMaxLength))
      ->check(CLI::Range(static_cast<u16>(kHpMinLength), static_cast<u16>(kHpMaxLength)))
      ->group(kOptGroupNameHpMetricsOptions)
      ->needs(hp_metrics_flag)
      ->check([&opts](const std::string& input) -> std::string {
        // Parse the input value first. This step should not throw because
        // the CLI::Range check above should have already validated the input.
        const auto max_value = static_cast<u16>(std::stoul(input));
        // Check constraint only if HP metrics are enabled and both values are set
        if (opts->calculate_hp_metrics && opts->min_hp_length > max_value) {
          throw CLI::ValidationError(fmt::format(
              "min-hp-length ({}) must be less than or equal to max-hp-length ({})", opts->min_hp_length, max_value));
        }
        // Return empty string if validation passes
        return "";
      });

  app->add_flag("--hp-allow-heterogeneous-insertions",
                opts->hp_allow_heterogeneous_insertions,
                "Allow insertions that are heterogeneous (e.g. an insertion of T into an AAAAAA homopolymer) when "
                "evaluating homopolymer errors.")
      ->default_val(false)
      ->group(kOptGroupNameHpMetricsOptions)
      ->needs(hp_metrics_flag);
  app->add_option(
         "--hp-subsampling-fraction",
         opts->hp_subsampling_fraction,
         "Consider only a specified fraction of all reference homopolymer regions to reduce runtime and memory usage.")
      ->default_val(kDefaultHpSubsamplingFraction)
      ->check(CLI::Range(0.0, 1.0))
      ->group(kOptGroupNameHpMetricsOptions);
  app->add_option(
         "--hp-subsampling-seed", opts->hp_subsampling_seed, "The seed used for subsampling homopolymer regions.")
      ->group(kOptGroupNameHpMetricsOptions);
}

void AddHpMaskingOption(cli::AppPtr app, const CoverageMetricsOptsPtr& opts) {
  app->add_flag("--disable-discordant-hp-masking",
                opts->disable_hp_quality_modification,
                "Disables masking of homopolymer regions with discordant duplex bases (only for duplex reads)")
      ->group(kOptGroupNameHpMaskingOptions);
  app->add_option("--min-baseq-for-hp-masking",
                  opts->base_quality_threshold_for_hp_masking,
                  "Base quality threshold for masking homopolymers for non-duplex dataset or when "
                  "--disable-base-type-decoding is set.")
      ->default_val(std::to_string(kDefaultBaseQualityThresholdForHpMasking))
      ->check(CLI::Range(kDefaultMinBaseq, kMaxBaseQuality))
      ->excludes("--disable-discordant-hp-masking")
      ->group(kOptGroupNameHpMaskingOptions);
}

void AddAccuracyMetricsOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts) {
  app->add_option("--min-depth",
                  opts->min_depth,
                  "Minimum depth threshold for errors at the base to contribute to the accuracy metrics.")
      ->default_val(std::to_string(kDefaultMinDepth))
      ->group(kOptGroupNameAccuracyMetricsOptions);
  app->add_option("--max-alt-allele-fraction",
                  opts->max_alt_allele_fraction,
                  "Threshold for the alternate allele fraction at a position to distinguish between germline variants "
                  "and sequencing errors. "
                  "When the alt allele fraction at a position exceeds this threshold, those errors are considered "
                  "likely germline variants "
                  "and do not contribute to accuracy metrics. For HP metrics, positions with depth less than the "
                  "threshold are not filtered out but are not considered for the alt allele fraction filter.")
      ->default_val(std::to_string(kDefaultMaxAltAlleleFraction))
      ->check(CLI::Range(0.0, 1.0))
      ->group(kOptGroupNameAccuracyMetricsOptions);
  app->add_option(
         "--max-cluster-size-bin",
         opts->max_cluster_size_bin,
         "Maximum cluster size for computing errors by cluster size. Errors originating in a consensus read generated "
         "from a cluster of size larger than the threshold will be put into the max-cluster-size+ bin.")
      ->default_val(kDefaultMaxClusterSizeBin)
      ->check(CLI::Range(kMinClusterSize, kMaxClusterSize))
      ->group(kOptGroupNameAccuracyMetricsOptions);
}

static void AddCoverageMetricsOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts) {
  app->add_option("--coverage-cutoffs",
                  opts->coverage_cutoff,
                  "A list of cutoff coverage values to be included in the coverage distribution summary metric.")
      ->delimiter(',')
      ->default_val(kDefaultCoverageCutoffs)
      ->group(kOptGroupNameCoverageMetricsOptions);
  app->add_flag("--exclude-uncovered-positions",
                opts->exclude_empty_positions,
                "Exclude positions with no coverage when calculating coverage metrics and histograms. If this flag is "
                "turned on, the 0 bin in coverage histogram will only include reference positions where all supporting "
                "reads have a deletion.")
      ->group(kOptGroupNameCoverageMetricsOptions);
  app->add_option(
         "--max-coverage-bin",
         opts->max_coverage_bin,
         "Maximum depth bin for coverage histograms. Positions with coverage greater than the threshold will only "
         "contribute to the `max-coverage+` bin.")
      ->default_val(std::to_string(kDefaultMaxCoverageBin))
      ->group(kOptGroupNameCoverageMetricsOptions);
}

void AddReadMetricsOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts) {
  app->add_option("--max-read-length-bin",
                  opts->max_read_length_bin,
                  "Maximum read length to consider for read length histograms. Reads longer than the threshold will "
                  "only contribute to the max-read-length+ bin.")
      ->default_val(kDefaultMaxReadLengthBin)
      ->group(kOptGroupNameReadMetricsOptions);
}

static void AddSummaryStatsOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts) {
  app->add_option("--summary-stats-percentiles",
                  opts->summary_stats_percentiles,
                  "A list of percentiles to compute for coverage and read length summary statistics.")
      ->delimiter(',')
      ->default_val(kDefaultSummaryStatsPercentiles)
      ->group(kOptGroupNameSummaryStatsOptions);
}

void AddPerformanceOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts, const bool ref_required) {
  app->add_option("--region-size", opts->region_size, "The size of each region that are processed in parallel.")
      ->default_val(kDefaultRegionSize)
      ->check(CLI::Range(u32{1}, std::numeric_limits<u32>::max()))
      ->group(kOptGroupNamePerformanceOptions);
  if (ref_required) {
    app->add_option("--reference-padding",
                    opts->reference_padding,
                    "The number of extra reference bases to fetch on the 3' side of the region. This is to ensure all "
                    "alignments that are fetched have the necessary context even if they extend beyond the region.")
        ->default_val(kDefaultMaxReferenceBufferLength)
        ->group(kOptGroupNamePerformanceOptions);
  }
  AddThreadCountOption(app, "--threads", opts->threads)->default_val(1)->group(kOptGroupNamePerformanceOptions);
}

void AddTEMetricsOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts) {
  app->add_flag("--enable-te-metrics",
                opts->enable_te_metrics,
                "Enable calculation of target enrichment read metrics. Requires a BED file specifying target regions.")
      ->group(kOptGroupNameTEMetricsOptions)
      ->needs("--bed-input");
}

void AddLegacyOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts) {
  app->add_flag("--ignore-family", opts->ignore_family, "Ignore family information when calculating metrics")
      ->group("Legacy Options");
}

void DefineCoverageMetricsOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts) {
  AddInputOptions(app, opts, false);
  AddOutputOptions(app, opts);
  AddReadFilterOptions(app, opts);
  AddBaseFilterOptions(app, opts);
  AddCoverageMetricsOptions(app, opts);
  AddSummaryStatsOptions(app, opts);
  AddHpMetricsOptions(app, opts);
  AddHpMaskingOption(app, opts);
  AddReadTrimmingOptions(app, opts);
  AddTEMetricsOptions(app, opts);
  AddPerformanceOptions(app, opts, true);
}

void DefineAccuracyMetricsOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts) {
  AddInputOptions(app, opts, true);
  AddOutputOptions(app, opts);
  AddReadFilterOptions(app, opts);
  AddBaseFilterOptions(app, opts);
  AddAccuracyMetricsOptions(app, opts);
  AddHpMaskingOption(app, opts);
  AddHpMetricsOptions(app, opts);
  AddPerformanceOptions(app, opts, true);
  AddReadTrimmingOptions(app, opts);
  AddLegacyOptions(app, opts);
}

void DefineReadMetricsOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts) {
  AddInputOptions(app, opts, false);
  // We don't need the reference for read metrics, so we can remove it
  app->remove_option(app->get_option("--reference"));
  AddOutputOptions(app, opts);
  AddReadFilterOptions(app, opts);
  AddReadMetricsOptions(app, opts);
  AddSummaryStatsOptions(app, opts);
  AddPerformanceOptions(app, opts, false);
  AddReadTrimmingOptions(app, opts);
  AddTEMetricsOptions(app, opts);
  AddLegacyOptions(app, opts);
}

void DefineAllMetricsOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts) {
  AddInputOptions(app, opts, true);
  AddOutputOptions(app, opts);
  AddReadFilterOptions(app, opts);
  AddBaseFilterOptions(app, opts);
  AddHpMaskingOption(app, opts);
  AddCoverageMetricsOptions(app, opts);
  AddAccuracyMetricsOptions(app, opts);
  AddReadMetricsOptions(app, opts);
  AddHpMetricsOptions(app, opts);
  AddSummaryStatsOptions(app, opts);
  AddPerformanceOptions(app, opts, true);
  AddReadTrimmingOptions(app, opts);
  AddTEMetricsOptions(app, opts);
  AddLegacyOptions(app, opts);
}

// This adds version and commandline to the comments field of the option
// This handles parent or subcommands appropriately if available.
static void AddVersionAndCommandLineComment(const cli::ConstAppPtr& app,
                                            const CoverageMetricsOptsPtr& opts,
                                            const std::string& version) {
  const std::string cli_args = cli::RenderFullCli(app);
  io::AddVersionAndCommandLineComment(opts->comments, version, cli_args);
}

static void AddSubcommands(const std::shared_ptr<CLI::App>& app,
                           CoverageMetricsOptsPtr& opts,
                           const std::string& version) {
  // add coverage metrics subcommand
  cli::AddSubcommand<AlignmentMetricsOptions>(
      app.get(),
      kSubcommandNameCoverageMetrics,
      DefineCoverageMetricsOptions,
      opts,
      RunAlignmentMetrics,
      kSubcommandDescriptionCoverageMetrics,
      [&version](const cli::ConstAppPtr& app, const CoverageMetricsOptsPtr& options) {
        options->metric_types.has_accuracy_metrics = false;
        options->metric_types.has_coverage_metrics = true;
        options->metric_types.has_read_metrics = false;
        AddVersionAndCommandLineComment(app, options, version);
      })
      ->fallthrough();

  // add accuracy metrics subcommand
  cli::AddSubcommand<AlignmentMetricsOptions>(
      app.get(),
      kSubcommandNameAccuracyMetrics,
      DefineAccuracyMetricsOptions,
      opts,
      RunAlignmentMetrics,
      kSubcommandDescriptionAccuracyMetrics,
      [&version](const cli::ConstAppPtr& app, const CoverageMetricsOptsPtr& options) {
        options->metric_types.has_accuracy_metrics = true;
        options->metric_types.has_coverage_metrics = false;
        options->metric_types.has_read_metrics = false;
        AddVersionAndCommandLineComment(app, options, version);
      })
      ->fallthrough();

  // add read metrics subcommand
  cli::AddSubcommand<AlignmentMetricsOptions>(
      app.get(),
      kSubcommandNameReadMetrics,
      DefineReadMetricsOptions,
      opts,
      RunAlignmentMetrics,
      kSubcommandDescriptionReadMetrics,
      [&version](const cli::ConstAppPtr& app, const CoverageMetricsOptsPtr& options) {
        options->metric_types.has_accuracy_metrics = false;
        options->metric_types.has_coverage_metrics = false;
        options->metric_types.has_read_metrics = true;
        AddVersionAndCommandLineComment(app, options, version);
      })
      ->fallthrough();

  // add all metrics subcommand
  cli::AddSubcommand<AlignmentMetricsOptions>(
      app.get(),
      kSubcommandNameAlignmentMetrics,
      // define options
      DefineAllMetricsOptions,
      // shared_ptr
      opts,
      // subcommand main
      RunAlignmentMetrics,
      kSubcommandDescriptionAlignmentMetrics,
      [&version](const cli::ConstAppPtr& app, const CoverageMetricsOptsPtr& options) {
        options->metric_types.has_accuracy_metrics = true;
        options->metric_types.has_coverage_metrics = true;
        options->metric_types.has_read_metrics = true;
        AddVersionAndCommandLineComment(app, options, version);
      })
      ->fallthrough();
}

int RunAlignmentMetricsApp(int argc, char** argv, const std::string& program_name, const std::string& version) {
  auto opts = std::make_shared<AlignmentMetricsOptions>();
  const auto app = cli::SetupDefaultCli(program_name, version);
  app->require_subcommand(1);
  AddSubcommands(app, opts, version);
  return cli::RunCli(app.get(), argc, argv);
}

}  // namespace xoos::alignment_metrics
