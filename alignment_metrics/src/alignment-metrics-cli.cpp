#include "alignment-metrics-cli.h"

#include <string>

#include <fmt/format.h>

#include <xoos/cli/bam-option-util.h>
#include <xoos/cli/cli.h>
#include <xoos/cli/enum-option-util.h>
#include <xoos/cli/fasta-option-util.h>
#include <xoos/cli/thread-count-option-util.h>
#include <xoos/cli/validators/bed-validator.h>
#include <xoos/cli/validators/filename-component-validator.h>
#include <xoos/cli/validators/indexed-bam-validator.h>
#include <xoos/cli/validators/indexed-fasta-validator.h>
#include <xoos/enum/enum-util.h>
#include <xoos/types/int.h>

#include "CLI/CLI.hpp"
#include "alignment-metrics-options.h"
#include "run-alignment-metrics.h"

namespace xoos::alignment_metrics {

namespace {
void ValidateGcBiasOptionConstraints(const AlignmentMetricsOptsPtr& opts) {
  if (opts->metric_types.has_gc_bias_metrics) {
    if (opts->region_size > 0 && opts->gc_scan_window_size >= opts->region_size) {
      throw CLI::ValidationError("--gc-window-size must be smaller than --region-size");
    }
    if (opts->reference_padding == 0) {
      throw CLI::ValidationError("--reference-padding must be greater than 0 when GC bias metrics are enabled");
    }
    if (opts->reference_padding > 0 && opts->gc_scan_window_size >= opts->reference_padding) {
      throw CLI::ValidationError("--gc-window-size must be smaller than --reference-padding");
    }
  }
}
}  // namespace

using cli::AddBamFileOption;
using cli::AddFastaFileOption;
using cli::AddThreadCountOption;

void AddInputOptions(CLI::App* const app, const AlignmentMetricsOptsPtr& opts) {
  AddBamFileOption(app,
                   "--bam",
                   opts->bam_input,
                   opts->bam_index,
                   "Path to the BAM file which the metrics are to be calculated. The BAM file must be "
                   "coordinate-sorted and indexed.")
      ->required()
      ->check(cli::IndexedBamFileValidator())
      ->group(kOptGroupNameInputOptions);
  AddFastaFileOption(
      app,
      "--reference",
      opts->reference_input,
      "Path to the reference FASTA file used to compare the reads against. The FASTA file must be indexed with a "
      "`.fai` file.")
      ->required()
      ->check(cli::IndexedFastaFileValidator())
      ->group(kOptGroupNameInputOptions);
}

void AddOutputOptions(CLI::App* const app, const AlignmentMetricsOptsPtr& opts) {
  app->add_option("--output-dir", opts->out_dir, "Path to the output directory where the metrics will be written to.")
      ->default_val(".")
      ->group(kOptGroupNameOutputOptions);
  app->add_option("--output-prefix",
                  opts->output_prefix,
                  "Prefix prepended to all output filenames (e.g. a sample ID), joined with an underscore. "
                  "If unset, default filenames are used.")
      ->check(cli::FilenameComponentValidator())
      ->group(kOptGroupNameOutputOptions);
  app->add_flag("--overwrite",
                opts->overwrite,
                "Overwrite existing output files. Without this flag, the tool errors out if any output file already "
                "exists.")
      ->group(kOptGroupNameOutputOptions);
}

void AddRegionOptions(CLI::App* const app, const AlignmentMetricsOptsPtr& opts) {
  app->add_option(
         "--bed",
         opts->bed_input,
         "Path to a BED file which specifies regions of interest for the metrics calculation. If specified, only reads "
         "overlapping with a region of interest and positions contained within a region of interest will contribute to "
         "the metrics. If not specified, all contigs present in the BAM header will contribute to the metrics. The "
         "file may be gzip/bgzip compressed when named with a `.gz` extension.")
      ->group(kOptGroupNameRegionOptions)
      ->check(cli::BedFileValidator());
  app->add_flag(
         "--stratify-by-label",
         opts->stratify_by_label,
         "Interpret the 4th column of the BED file as a comma-separated list of region labels and stratify the output "
         "by "
         "label. Metrics for each label are written to a separate subdirectory of the output directory (e.g. "
         "<output-dir>/<label>/{accuracy,coverage,read_metrics}). Every region in the BED file must have a 4th column. "
         "At most 16 distinct labels are supported. Requires --bed.")
      ->group(kOptGroupNameRegionOptions)
      ->needs("--bed");
}

void AddReadFilteringOptions(CLI::App* const app, const AlignmentMetricsOptsPtr& opts) {
  app->add_option("--exclude-flags",
                  opts->exclude_flags,
                  "Exclude alignments from coverage metrics, accuracy metrics, and post-filter read metrics if any "
                  "bits in their FLAG field match the specified decimal integer. Refer to "
                  "https://broadinstitute.github.io/picard/explain-flags.html for flag generation. Default is "
                  "supplementary, secondary, unmapped, and duplicate (3332)")
      ->default_val(kDefaultExcludeFlag)
      ->group(kOptGroupNameReadFilteringOptions);
  app->add_option("--min-mapping-quality",
                  opts->min_mapq,
                  "Minimum mapping quality. Reads with a mapping quality less than the threshold are not considered "
                  "for coverage metrics, accuracy metrics, and post-filter read metrics.")
      ->default_val(std::to_string(kDefaultMinMapq))
      ->check(CLI::Range(kMinBaseQuality, kMaxBaseQuality))
      ->group(kOptGroupNameReadFilteringOptions);
  app->add_option(
         "--trim-start-length", opts->trim_leading_bases, "Number of bases to trim from the start of each read")
      ->default_val(0)
      ->check(CLI::Range(0U, kMaxTrimLeadingBases))
      ->group(kOptGroupNameReadFilteringOptions);
  app->add_option("--trim-end-length", opts->trim_trailing_bases, "Number of bases to trim from the end of each read")
      ->default_val(0)
      ->check(CLI::Range(0U, kMaxTrimTrailingBases))
      ->group(kOptGroupNameReadFilteringOptions);
}

void AddBaseFilteringOptions(CLI::App* const app, const AlignmentMetricsOptsPtr& opts) {
  app->add_option("--min-base-quality",
                  opts->min_baseq,
                  "Minimum base quality. Bases with a base quality less than the threshold are not considered for "
                  "accuracy metrics.")
      ->default_val(std::to_string(kDefaultMinBaseq))
      ->check(CLI::Range(0U, kMaxBaseQuality))
      ->group(kOptGroupNameBaseFilteringOptions);
  app->add_flag("--skip-base-type-decoding",
                opts->disable_base_type_decoding,
                "Skip decoding base types from the YC tag in the BAM file. This can speed up processing for simplex "
                "datasets where base type decoding is not necessary.")
      ->group(kOptGroupNameBaseFilteringOptions);
  app->add_flag(
         "--include-ambiguous-reference-bases",
         opts->include_ambiguous_reference_bases,
         "Include IUPAC ambiguous reference bases (R, Y, S, W, K, M, B, D, H, V) in coverage and accuracy metrics. "
         "When enabled, query bases matched against an ambiguous reference that is not N are not counted as errors "
         "when the query base is consistent with the ambiguity code. Only N/n reference positions are excluded. "
         "By default, all non-ACGT reference positions are excluded.")
      ->group(kOptGroupNameBaseFilteringOptions);
  cli::AddOptionalEnumOption(
      app,
      "--min-base-type",
      opts->min_base_type,
      "Minimum duplex base type. For duplex datasets, bases with a base type worse than the threshold are not "
      "considered for accuracy and post-filter coverage metrics. For simplex datasets, it should always be set to "
      "`simplex`.")
      ->default_val(enum_util::FormatEnumName(kDefaultMinBaseType))
      ->group(kOptGroupNameBaseFilteringOptions)
      ->check([&opts](const std::string& input) -> std::string {
        const auto base_type = cli::ParseEnumNameOrThrow<yc_decode::BaseType>("min-base-type", input);
        if (opts->disable_base_type_decoding && base_type == yc_decode::BaseType::kConcordant) {
          throw CLI::ValidationError(
              "When --skip-base-type-decoding is enabled, --min-base-type must not be set to `concordant` as no "
              "base type "
              "information will be available and all bases will be treated as simplex. Setting --min-base-type to "
              "`concordant` in this case will cause all bases to be filtered out.");
        }
        return "";
      });
}

void AddHomopolymerMetricOptions(CLI::App* const app, const AlignmentMetricsOptsPtr& opts) {
  auto* const hp_metrics_flag =
      app->add_flag(
             "--enable-homopolymer-metrics", opts->calculate_hp_metrics, "Enable calculation of homopolymer metrics.")
          ->group(kOptGroupNameHomopolymerMetricOptions)
          ->default_val(false)
          ->needs("--reference");

  // Need to cast to u16 for range check because CLI11 treats u8 as
  // unsigned char and it won't render correctly
  app->add_option("--min-homopolymer-length",
                  opts->min_hp_length,
                  "Minimum homopolymer length to consider for homopolymer metrics.")
      ->default_val(static_cast<u16>(kHpMinLength))
      ->check(CLI::Range(static_cast<u16>(kHpMinLength), static_cast<u16>(kHpMaxLength)))
      ->group(kOptGroupNameHomopolymerMetricOptions)
      ->needs(hp_metrics_flag);

  app->add_option("--max-homopolymer-length",
                  opts->max_hp_length,
                  "Maximum homopolymer length to consider for homopolymer metrics.")
      ->default_val(static_cast<u16>(kHpMaxLength))
      ->check(CLI::Range(static_cast<u16>(kHpMinLength), static_cast<u16>(kHpMaxLength)))
      ->group(kOptGroupNameHomopolymerMetricOptions)
      ->needs(hp_metrics_flag)
      ->check([&opts](const std::string& input) -> std::string {
        // Parse the input value first. This step should not throw because
        // the CLI::Range check above should have already validated the input.
        const auto max_value = static_cast<u16>(std::stoul(input));
        // Check constraint only if HP metrics are enabled and both values are set
        if (opts->calculate_hp_metrics && opts->min_hp_length > max_value) {
          throw CLI::ValidationError(
              fmt::format("min-homopolymer-length ({}) must be less than or equal to max-homopolymer-length ({})",
                          opts->min_hp_length,
                          max_value));
        }
        // Return empty string if validation passes
        return "";
      });

  app->add_flag("--homopolymer-allow-heterogeneous-insertions",
                opts->hp_allow_heterogeneous_insertions,
                "Allow insertions that are heterogeneous (e.g. an insertion of T into an AAAAAA homopolymer) when "
                "evaluating homopolymer errors.")
      ->default_val(false)
      ->group(kOptGroupNameHomopolymerMetricOptions)
      ->needs(hp_metrics_flag);
}

void AddHpMaskingOption(CLI::App* const app, const AlignmentMetricsOptsPtr& opts) {
  app->add_flag("--skip-discordant-homopolymer-masking",
                opts->disable_hp_quality_modification,
                "Disables masking of homopolymer regions with discordant duplex bases (only for duplex reads)")
      ->group(kOptGroupNameHpMaskingOptions);
  app->add_option("--min-homopolymer-masking-base-quality",
                  opts->base_quality_threshold_for_hp_masking,
                  "Base quality threshold for masking homopolymers for non-duplex dataset or when "
                  "--skip-base-type-decoding is set.")
      ->default_val(std::to_string(kDefaultBaseQualityThresholdForHpMasking))
      ->check(CLI::Range(kMinBaseQuality, kMaxBaseQuality))
      ->excludes("--skip-discordant-homopolymer-masking")
      ->group(kOptGroupNameHpMaskingOptions);
}

void AddAccuracyMetricsOptions(CLI::App* const app, const AlignmentMetricsOptsPtr& opts) {
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

void AddCoverageMetricsOptions(CLI::App* const app, const AlignmentMetricsOptsPtr& opts) {
  app->add_option("--coverage-distribution-cutoffs",
                  opts->coverage_cutoff,
                  "A list of cutoff coverage values to be included in the coverage distribution summary metric.")
      ->delimiter(',')
      ->default_val(kDefaultCoverageCutoffs)
      ->group(kOptGroupNameCoverageMetricsOptions);
  app->add_flag("--skip-uncovered-positions",
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
  app->add_flag("--enable-gc-bias-metrics",
                opts->metric_types.has_gc_bias_metrics,
                "Enable calculation of GC bias coverage metrics. Requires a reference sequence.")
      ->needs("--reference")
      ->group(kOptGroupNameCoverageMetricsOptions);
  app->add_option("--gc-window-size",
                  opts->gc_scan_window_size,
                  "The size of the sliding window used to compute GC content for GC bias metrics. "
                  "Must be smaller than --region-size and --reference-padding.")
      ->default_val(kDefaultGcWindowSize)
      ->check(CLI::Range(u32{1}, kMaxGcWindowSize))
      ->group(kOptGroupNameCoverageMetricsOptions);
}

void AddGcBiasOptions(CLI::App* const app, const AlignmentMetricsOptsPtr& opts) {
  app->add_option("--gc-window-size",
                  opts->gc_scan_window_size,
                  "The size of the sliding window used to compute GC content for GC bias metrics. "
                  "Must be smaller than --region-size and --reference-padding.")
      ->default_val(kDefaultGcWindowSize)
      ->check(CLI::Range(u32{1}, kMaxGcWindowSize))
      ->group(kOptGroupNameGcBiasMetricsOptions);
}

void AddTargetEnrichmentMetricsOptions(CLI::App* const app, const AlignmentMetricsOptsPtr& opts) {
  app->add_flag("--enable-target-enrichment-metrics",
                opts->enable_te_metrics,
                "Enable calculation of target enrichment read metrics. Requires a BED file specifying target regions.")
      ->group(kOptGroupNameTargetEnrichmentMetricsOptions)
      ->needs("--bed");
}

void AddReadMetricsOptions(CLI::App* const app, const AlignmentMetricsOptsPtr& opts) {
  app->add_option("--max-read-length-bin",
                  opts->max_read_length_bin,
                  "Maximum read length to consider for read length histograms. Reads longer than the threshold will "
                  "only contribute to the max-read-length+ bin.")
      ->default_val(kDefaultMaxReadLengthBin)
      ->group(kOptGroupNameReadMetricsOptions);
}

void AddGeneralMetricsOptions(CLI::App* const app, const AlignmentMetricsOptsPtr& opts) {
  app->add_option("--summary-stats-percentiles",
                  opts->summary_stats_percentiles,
                  "A list of percentiles to compute for coverage and read length summary statistics.")
      ->delimiter(',')
      ->default_val(kDefaultSummaryStatsPercentiles)
      ->group(kOptGroupNameGeneralMetricsOptions);
}

void AddPerformanceOptions(CLI::App* const app, const AlignmentMetricsOptsPtr& opts, const bool ref_required) {
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

void AddHiddenOptions(CLI::App* const app,
                      const AlignmentMetricsOptsPtr& opts,
                      const bool include_hp_subsampling_options) {
  app->add_flag("--ignore-family", opts->ignore_family, "Ignore family information when calculating metrics")
      ->group("");
  if (include_hp_subsampling_options) {
    app->add_option("--homopolymer-subsample-fraction",
                    opts->hp_subsampling_fraction,
                    "Consider only a specified fraction of all reference homopolymer regions to reduce runtime and "
                    "memory usage.")
        ->default_val(kDefaultHpSubsamplingFraction)
        ->check(CLI::Range(0.0, 1.0))
        ->group("");
    app->add_option("--homopolymer-subsample-seed",
                    opts->hp_subsampling_seed,
                    "The seed used for subsampling homopolymer regions.")
        ->group("");
  }
}

void DefineCoverageMetricsOptions(CLI::App* const app, const AlignmentMetricsOptsPtr& opts) {
  AddInputOptions(app, opts);
  AddOutputOptions(app, opts);
  AddRegionOptions(app, opts);
  AddReadFilteringOptions(app, opts);
  AddBaseFilteringOptions(app, opts);
  AddCoverageMetricsOptions(app, opts);
  AddHomopolymerMetricOptions(app, opts);
  AddHpMaskingOption(app, opts);
  AddTargetEnrichmentMetricsOptions(app, opts);
  AddGeneralMetricsOptions(app, opts);
  AddPerformanceOptions(app, opts, true);
  AddHiddenOptions(app, opts, true);
}

void DefineAccuracyMetricsOptions(CLI::App* const app, const AlignmentMetricsOptsPtr& opts) {
  AddInputOptions(app, opts);
  AddOutputOptions(app, opts);
  AddRegionOptions(app, opts);
  AddReadFilteringOptions(app, opts);
  AddBaseFilteringOptions(app, opts);
  AddAccuracyMetricsOptions(app, opts);
  AddHomopolymerMetricOptions(app, opts);
  AddHpMaskingOption(app, opts);
  // No TE metrics for accuracy
  AddGeneralMetricsOptions(app, opts);
  AddPerformanceOptions(app, opts, true);
  AddHiddenOptions(app, opts, true);
}

void DefineReadMetricsOptions(CLI::App* const app, const AlignmentMetricsOptsPtr& opts) {
  AddInputOptions(app, opts);
  // We don't need the reference for read metrics, so we can remove it
  app->remove_option(app->get_option("--reference"));
  AddOutputOptions(app, opts);
  AddRegionOptions(app, opts);
  AddReadFilteringOptions(app, opts);
  // No base filtering, HP, or HP masking for read metrics
  AddReadMetricsOptions(app, opts);
  AddTargetEnrichmentMetricsOptions(app, opts);
  AddGeneralMetricsOptions(app, opts);
  AddPerformanceOptions(app, opts, false);
  AddHiddenOptions(app, opts, false);
}

void DefineGcBiasMetricsOptions(CLI::App* const app, const AlignmentMetricsOptsPtr& opts) {
  AddInputOptions(app, opts);
  AddOutputOptions(app, opts);
  AddRegionOptions(app, opts);
  AddReadFilteringOptions(app, opts);
  AddGcBiasOptions(app, opts);
  AddPerformanceOptions(app, opts, true);
}

void DefineAllMetricsOptions(CLI::App* const app, const AlignmentMetricsOptsPtr& opts) {
  AddInputOptions(app, opts);
  AddOutputOptions(app, opts);
  AddRegionOptions(app, opts);
  AddReadFilteringOptions(app, opts);
  AddBaseFilteringOptions(app, opts);
  AddCoverageMetricsOptions(app, opts);
  AddAccuracyMetricsOptions(app, opts);
  AddReadMetricsOptions(app, opts);
  AddHomopolymerMetricOptions(app, opts);
  AddHpMaskingOption(app, opts);
  AddTargetEnrichmentMetricsOptions(app, opts);
  AddGeneralMetricsOptions(app, opts);
  AddPerformanceOptions(app, opts, true);
  AddHiddenOptions(app, opts, true);
}

namespace {
void AddSubcommands(const std::shared_ptr<CLI::App>& app, AlignmentMetricsOptsPtr& opts) {
  // add coverage metrics subcommand
  cli::AddSubcommand<AlignmentMetricsOptions>(app.get(),
                                              kSubcommandNameCoverageMetrics,
                                              DefineCoverageMetricsOptions,
                                              opts,
                                              RunAlignmentMetrics,
                                              kSubcommandDescriptionCoverageMetrics,
                                              [](const CLI::App* const app, const AlignmentMetricsOptsPtr& options) {
                                                options->metric_types.has_accuracy_metrics = false;
                                                options->metric_types.has_coverage_metrics = true;
                                                options->metric_types.has_read_metrics = false;
                                                options->command_line_info = cli::GetCommandLineInfo(app);
                                                ValidateGcBiasOptionConstraints(options);
                                              })
      ->fallthrough();

  // add accuracy metrics subcommand
  cli::AddSubcommand<AlignmentMetricsOptions>(app.get(),
                                              kSubcommandNameAccuracyMetrics,
                                              DefineAccuracyMetricsOptions,
                                              opts,
                                              RunAlignmentMetrics,
                                              kSubcommandDescriptionAccuracyMetrics,
                                              [](const CLI::App* const app, const AlignmentMetricsOptsPtr& options) {
                                                options->metric_types.has_accuracy_metrics = true;
                                                options->metric_types.has_coverage_metrics = false;
                                                options->metric_types.has_read_metrics = false;
                                                options->command_line_info = cli::GetCommandLineInfo(app);
                                              })
      ->fallthrough();

  // add read metrics subcommand
  cli::AddSubcommand<AlignmentMetricsOptions>(app.get(),
                                              kSubcommandNameReadMetrics,
                                              DefineReadMetricsOptions,
                                              opts,
                                              RunAlignmentMetrics,
                                              kSubcommandDescriptionReadMetrics,
                                              [](const CLI::App* const app, const AlignmentMetricsOptsPtr& options) {
                                                options->metric_types.has_accuracy_metrics = false;
                                                options->metric_types.has_coverage_metrics = false;
                                                options->metric_types.has_read_metrics = true;
                                                options->command_line_info = cli::GetCommandLineInfo(app);
                                              })
      ->fallthrough();

  // add gc bias metrics subcommand
  cli::AddSubcommand<AlignmentMetricsOptions>(app.get(),
                                              kSubcommandNameGcBiasMetrics,
                                              DefineGcBiasMetricsOptions,
                                              opts,
                                              RunAlignmentMetrics,
                                              kSubcommandDescriptionGcBiasMetrics,
                                              [](const CLI::App* const app, const AlignmentMetricsOptsPtr& options) {
                                                options->metric_types.has_accuracy_metrics = false;
                                                options->metric_types.has_coverage_metrics = false;
                                                options->metric_types.has_read_metrics = false;
                                                options->metric_types.has_gc_bias_metrics = true;
                                                options->command_line_info = cli::GetCommandLineInfo(app);
                                                ValidateGcBiasOptionConstraints(options);
                                              })
      ->fallthrough();

  // add all metrics subcommand
  cli::AddSubcommand<AlignmentMetricsOptions>(app.get(),
                                              kSubcommandNameAlignmentMetrics,
                                              DefineAllMetricsOptions,
                                              opts,
                                              RunAlignmentMetrics,
                                              kSubcommandDescriptionAlignmentMetrics,
                                              [](const CLI::App* const app, const AlignmentMetricsOptsPtr& options) {
                                                options->metric_types.has_accuracy_metrics = true;
                                                options->metric_types.has_coverage_metrics = true;
                                                options->metric_types.has_read_metrics = true;
                                                options->command_line_info = cli::GetCommandLineInfo(app);
                                                ValidateGcBiasOptionConstraints(options);
                                              })
      ->fallthrough();
}
}  // namespace

int RunAlignmentMetricsApp(int argc, char** argv, const std::string& program_name, const std::string& version) {
  auto opts = std::make_shared<AlignmentMetricsOptions>();
  const auto app = cli::SetupDefaultCli(program_name, version);
  app->require_subcommand(1);
  AddSubcommands(app, opts);
  return cli::RunCli(app.get(), argc, argv);
}

}  // namespace xoos::alignment_metrics
