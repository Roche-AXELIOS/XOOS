#include "util/cli-option-util.h"

#include <xoos/cli/enum-option-util.h>
#include <xoos/cli/thread-count-option-util.h>
#include <xoos/cli/validators/file-extension-validator.h>
#include <xoos/enum/enum-util.h>
#include <xoos/types/int.h>
#include <xoos/util/container-functions.h>
#include <xoos/util/file-functions.h>
#include <xoos/util/hash.h>

#include "CLI/CLI.hpp"
#include "consensus/qscore-calculator.h"
#include "core/read-collapser-options.h"

namespace xoos::read_collapser {

using cli::AddThreadCountOption;
using cli::ParseEnumNameOrThrow;
using enum_util::FormatEnumName;

size_t PresetHash::operator()(const ReadCollapserPresets& presets) const {
  return util::hash::Hash(presets);
}

static void UpdateDefaultsFromPreset(cli::AppPtr app,
                                     const PresetsMap& preset_map,
                                     const ReadCollapserPresets& preset) {
  const auto it = preset_map.find(preset);
  if (it != preset_map.end()) {
    for (const auto& [key, value] : it->second) {
      const auto opt = app->get_option(key);
      if (opt != nullptr) {
        opt->default_val(value);
      }
    }
  }
}

void AddPresetOption(cli::AppPtr app, const PresetsMap& preset_map, const std::string& group_name) {
  // Build the string list of available presets for the help menu
  auto presets = vec<ReadCollapserPresets>{};
  for (const auto& [preset, _] : preset_map) {
    presets.push_back(preset);
  }
  const std::string desc = fmt::format("value in {}", FormatEnumName(presets));
  const std::string main_desc =
      "Specify a preset to set default parameters. To see defaults for each preset, run the preset with `--help`. Any "
      "parameters provided that conflict with the preset configurations will take precedence over the preset defaults.";
  const std::string full_description = fmt::format("{}\n{}", desc, main_desc);

  // Define the preset option
  app->add_option_function<std::string>(
         "--preset",
         [app, preset_map](const std::string& value) {
           const auto preset = ParseEnumNameOrThrow<ReadCollapserPresets>("preset", value);
           UpdateDefaultsFromPreset(app, preset_map, preset);
         },
         full_description)
      ->trigger_on_parse()
      ->group(group_name);
}

void AddCommonInputOptions(cli::AppPtr app, const ReadCollapserOptionsPtr& options, const std::string& group_name) {
  app->add_option("--bam-input", options->bam_input, "Input BAM file, BAM must have an index.")
      ->check(cli::FileExtensionValidator({".bam"}))
      ->required()
      ->group(group_name);
  app->add_option("--bed-input", options->bed_input, "Regions to handle in BED format.")
      ->group(kOptGroupNameInputOptions);
  app->add_option("--padding",
                  options->padding,
                  "Number of bases of padding to apply to the regions specified by `--bed-input`.")
      ->check(CLI::NonNegativeNumber)
      ->default_val(kDefaultPadding)
      ->group(group_name);
}

void AddCommonOutputOptions(cli::AppPtr app, const ReadCollapserOptionsPtr& options, const std::string& group_name) {
  app->add_option("--output-dir", options->output_dir, "Directory for output file(s) and metrics.")
      ->default_val("output")
      ->group(group_name);
}

void AddCommonReadFilteringOptions(cli::AppPtr app,
                                   const ReadCollapserOptionsPtr& options,
                                   const std::string& group_name) {
  app->add_option("--min-mapq", options->min_mapq, "Minimum mapping quality.")
      ->check(CLI::Range(0u, kU8Limit))
      ->default_val(0)
      ->group(group_name);
  app->add_option(
         "--max-discordant-duplex-error-percentage",
         options->max_discordant_duplex_error_percentage,
         "Maximum percentage of discordant duplex bases allowed in a duplex read before the read is discarded.")
      ->check(CLI::Range(0.0, 100.0))
      ->group(group_name);
  app->add_option<bool>(
         "--exclude-partial-reads", options->exclude_partial_reads, "Exclude partial reads from analysis.")
      ->expected(0)
      ->default_val("false")
      ->group(group_name);
}

void AddConsensusReadFilterOptions(cli::AppPtr app,
                                   const ReadCollapserOptionsPtr& options,
                                   const std::string& group_name) {
  app->add_option("--exclude-flags",
                  options->exclude_flags,
                  "Exclude flags to exclude reads from analysis. "
                  "(BAM_FSUPPLEMENTARY | BAM_FSECONDARY).  See "
                  "https://broadinstitute.github.io/picard/explain-flags.html for how to "
                  "generate flags.")
      ->default_val(kDefaultExcludeFlagConsensus)
      ->group(group_name);
}

void AddMarkdupReadFilterOptions(cli::AppPtr app,
                                 const ReadCollapserOptionsPtr& options,
                                 const std::string& group_name) {
  app->add_option("--exclude-flags",
                  options->exclude_flags,
                  "Exclude flags to exclude reads from analysis. "
                  "(BAM_FSECONDARY).  See "
                  "https://broadinstitute.github.io/picard/explain-flags.html for how to "
                  "generate flags.")
      ->default_val(kDefaultExcludeFlagMarkdup)
      ->group(group_name);
}

void AddCommonClusterOptions(cli::AppPtr app, const ReadCollapserOptionsPtr& options, const std::string& group_name) {
  app->add_option<bool>("--cluster-by-umi", options->cluster_by_umi, "Split clusters by UMI.")
      ->expected(0)
      ->default_val("false")
      ->group(group_name);
  app->add_option<bool>("--cluster-by-strand", options->cluster_by_strand, "Split clusters by strand.")
      ->expected(0)
      ->default_val("false")
      ->group(group_name);
  app->add_option<bool>("--make-clusters-of-partial-reads-only",
                        options->make_clusters_of_partial_reads_only,
                        "Make clusters from remaining unassigned partial reads.")
      ->expected(0)
      ->default_val("false")
      ->group(group_name);
  app->add_option(
         "--wiggle-room", options->wiggle_room, "Maximum distance between read alignments to be in the same cluster.")
      ->check(CLI::Range(0u, kU8Limit))
      ->default_val(kDefaultWiggleRoom)
      ->group(group_name);
  app->add_option("--wiggle-room-partial",
                  options->wiggle_room_partial,
                  "Maximum distance between read alignments to be in the same cluster for partial reads.")
      ->check(CLI::Range(0u, kU8Limit))
      ->default_val(kDefaultWiggleRoomPartial)
      ->group(group_name);
}

void AddCommonPerformanceOptions(cli::AppPtr app,
                                 const ReadCollapserOptionsPtr& options,
                                 const std::string& group_name) {
  AddThreadCountOption(app, "--threads", options->threads)->check(CLI::NonNegativeNumber)->group(group_name);
  app->add_option("--region-size", options->region_size, "The size of each region (in bases) to process in parallel.")
      ->check(CLI::PositiveNumber)
      ->default_val(kDefaultRegionSize)
      ->group(group_name);
  app->add_option("--batch-size",
                  options->batch_size,
                  "Number of reads to process in a batch. Used in combination with `--region-size` to limit the number "
                  "of reads being processed simultaneously.")
      ->check(CLI::PositiveNumber)
      ->default_val(kDefaultBatchSize)
      ->group(group_name);
}

void AddConsensusOptions(cli::AppPtr app, const ReadCollapserOptionsPtr& options, const std::string& group_name) {
  app->add_option("--min-trim-read-support",
                  options->min_trim_read_support,
                  "The minimum number of reads required for a specific base position at the ends of the consensus "
                  "sequence. Any bases at "
                  "the start or end of the consensus sequence with read support below this threshold will be trimmed. "
                  "Must be <= --min-cluster-size.")
      ->check(CLI::PositiveNumber)
      ->default_val(kMinReadsPerCluster)
      ->group(group_name);
  app->add_option("--min-same-strand-cluster-size",
                  options->min_same_strand_cluster_size,
                  "Consensus sequences generated from same-strand clusters will be discarded if the effective cluster "
                  "size is below this threshold. The effective cluster size is calculated by the total number of bases "
                  "in the consensus matrix divided by the number of covered positions. "
                  "Defaults to the value of --min-cluster-size if not specified.")
      ->check(CLI::PositiveNumber)
      ->group(group_name);
  app->add_option("--min-mixed-strand-cluster-size",
                  options->min_mixed_strand_cluster_size,
                  "Consensus sequences generated from mixed-strand clusters will be discarded if the effective cluster "
                  "size is below this threshold. The effective cluster size is calculated by the total number of bases "
                  "in the consensus matrix divided by the number of covered positions. "
                  "Defaults to the value of --min-cluster-size if not specified.")
      ->check(CLI::PositiveNumber)
      ->group(group_name);
  app->add_option("--consensus-threshold",
                  options->consensus_threshold,
                  "The fraction of reads that must agree on a base call for it to be considered a strong consensus.")
      ->check(CLI::Range(0.0, 1.0))
      ->default_val(kDefaultConsensusThreshold)
      ->group(group_name);
  app->add_option("--consensus-gap-threshold",
                  options->consensus_gap_threshold,
                  "The fraction of reads at a given position that must support a gap for a gap to be called in the "
                  "consensus sequence.")
      ->check(CLI::Range(0.0, 1.0))
      ->default_val(kDefaultConsensusGapThreshold)
      ->group(group_name);
  app->add_option<bool>("--include-softclips",
                        options->include_softclips,
                        "Include soft-clipped segments for consensus generation. Only supported for data with UMIs. "
                        "Needs --cluster-by-umi.")
      ->expected(0)
      ->default_val("false")
      ->group(group_name);
  app->add_option<bool>("--enable-legacy-qscore-model",
                        options->enable_legacy_qscore_model,
                        "Use legacy Q-score model for consensus quality scoring.")
      ->expected(0)
      ->default_val("false")
      ->group(group_name);
  cli::AddEnumOption(app,
                     "--duplex-library-type",
                     options->duplex_library_type,
                     "Specifies the method for decoding duplex reads into their constituent single-strand reads (R1 "
                     "and R2). This mode is used to inform downstream consensus calling by treating R1 and R2 as "
                     "separate entities from a single molecule. None implies a simplex library.",
                     HDDeconvolutionType::kNone)
      ->group(group_name);
  app->add_option("--min-consensus-read-length",
                  options->min_consensus_read_length,
                  "Discard consensus reads shorter than the specified length. If not specified, no consensus reads are "
                  "discarded based on length.")
      ->check(CLI::PositiveNumber)
      ->group(group_name);
}

void AddMarkdupOutputOptions(cli::AppPtr app, const ReadCollapserOptionsPtr& options, const std::string& group_name) {
  app->add_option<bool>(
         "--remove-duplicates", options->remove_duplicates, "Remove duplicate reads from the output BAM file(s).")
      ->expected(0)
      ->default_val("false")
      ->group(group_name);
  app->add_option<bool>("--exclude-cluster-tags",
                        options->exclude_cluster_tags,
                        "Excludes cluster ID and size from the output BAM file(s).")
      ->expected(0)
      ->group(group_name);
  app->add_option<bool>(
         "--merge-output", options->merge_output, "Merge output BAM files into a single position-sorted BAM file.")
      ->expected(0)
      ->default_val("false")
      ->group(group_name);
}

void AddConsensusOutputOptions(cli::AppPtr app, const ReadCollapserOptionsPtr& options, const std::string& group_name) {
  app->add_option(
         "--compression-level", options->compression_level, "The level of compression used for the FASTQ output, 1-9.")
      ->check(CLI::Range(kMinCompressionLevel, kMaxCompressionLevel))
      ->default_val(kMinCompressionLevel)
      ->group(group_name);
  app->add_option<bool>("--output-cluster-bam", options->output_cluster_bam, "Output cluster BAM file(s).")
      ->expected(0)
      ->default_val("false")
      ->group(group_name);
}

void AddConsensusClusterOptions(cli::AppPtr app,
                                const ReadCollapserOptionsPtr& options,
                                const std::string& group_name) {
  app->add_option(
         "--max-cluster-size", options->max_cluster_size, "The maximum reads in a cluster before downsampling.")
      ->check(CLI::Range(1, QscoreCalculator::kMaxReadsInCluster))
      ->default_val(kMaxReadsPerCluster)
      ->group(group_name);
  app->add_option("--min-cluster-size",
                  options->min_cluster_size,
                  "The minimum reads in a cluster. Clusters with fewer reads will be discarded.")
      ->check(CLI::PositiveNumber)
      ->default_val(kMinReadsPerCluster)
      ->group(group_name);
}

void AddConsensusDebugOptions(cli::AppPtr app, const ReadCollapserOptionsPtr& options, const std::string& group_name) {
  app->add_option<bool>("--include-per-base-read-support-tags",
                        options->include_per_base_read_support_tags,
                        "If enabled, includes per-base read support information in the FASTQ output files.")
      ->expected(0)
      ->default_val("false")
      ->group(group_name);
  app->add_option<bool>("--include-per-base-majority-count-tags",
                        options->include_per_base_majority_count_tags,
                        "If enabled, includes per-base majority count information in the FASTQ output files.")
      ->expected(0)
      ->default_val("false")
      ->group(group_name);
}

void ValidateConsensusOptions(const cli::ConstAppPtr& app, const ReadCollapserOptionsPtr& options) {
  // If not explicitly provided (by user or preset), default strand cluster sizes to the value of --min-cluster-size
  if (options->min_same_strand_cluster_size == 0) {
    options->min_same_strand_cluster_size = options->min_cluster_size;
  }
  if (options->min_mixed_strand_cluster_size == 0) {
    options->min_mixed_strand_cluster_size = options->min_cluster_size;
  }
  // `--min-trim-read-support` must be <= `--min-cluster-size`
  if (options->min_trim_read_support > options->min_cluster_size) {
    throw CLI::ValidationError("--min-trim-read-support must be <= --min-cluster-size.");
  }
  // `--include-softclips` requires `--cluster-by-umi`
  // We need to check this in a callback because the default value of `--cluster-by-umi` may be changed by presets
  if (options->include_softclips && !options->cluster_by_umi) {
    throw CLI::ValidationError(
        "--include-softclips is only supported for data with UMIs and --cluster-by-umi must be enabled.");
  }
  // Manually set exclude_flags to default value for consensus to avoid it
  // being overwritten by markdup default value
  if (app->get_option("--exclude-flags")->count() == 0) {
    options->exclude_flags = kDefaultExcludeFlagConsensus;
  }
}

void ValidateMarkdupOptions(const cli::ConstAppPtr& app, const ReadCollapserOptionsPtr& options) {
  // Manually set exclude_flags to default value for markdup to avoid it
  // being overwritten by consensus default value
  if (app->get_option("--exclude-flags")->count() == 0) {
    options->exclude_flags = kDefaultExcludeFlagMarkdup;
  }
}

}  // namespace xoos::read_collapser
