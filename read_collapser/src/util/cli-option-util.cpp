#include "util/cli-option-util.h"

#include <xoos/cli/enum-option-util.h>
#include <xoos/cli/thread-count-option-util.h>
#include <xoos/cli/validators/bed-validator.h>
#include <xoos/cli/validators/indexed-bam-validator.h>
#include <xoos/enum/enum-util.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/util/container-functions.h>
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

namespace {

void UpdateDefaultsFromPreset(CLI::App* const app, const PresetsMap& preset_map, const ReadCollapserPresets& preset) {
  const auto it = preset_map.find(preset);
  if (it != preset_map.end()) {
    for (const auto& [key, value] : it->second) {
      auto* const opt = app->get_option(key);
      if (opt != nullptr) {
        opt->default_val(value);
      }
    }
  }
}

}  // namespace

void AddPresetOption(CLI::App* const app, const PresetsMap& preset_map, const std::string& group_name) {
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

void AddCommonInputOptions(CLI::App* const app, const ReadCollapserOptionsPtr& options, const std::string& group_name) {
  app->add_option("--bam", options->bam_input, "Input BAM file, BAM must have an index.")
      ->check(cli::IndexedBamFileValidator())
      ->required()
      ->group(group_name);
}

void AddCommonRegionOptions(CLI::App* const app,
                            const ReadCollapserOptionsPtr& options,
                            const std::string& group_name) {
  app->add_option("--bed", options->bed_input, "Regions to handle in BED format.")
      ->check(cli::BedFileValidator())
      ->group(group_name);
  app->add_option(
         "--padding", options->padding, "Number of bases of padding to apply to the regions specified by `--bed`.")
      ->check(CLI::NonNegativeNumber)
      ->default_val(kDefaultPadding)
      ->group(group_name);
}

void AddCommonHiddenOptions(CLI::App* const app, const ReadCollapserOptionsPtr& options) {
  app->add_option(
         "--ignore-read-name-parsing-errors",
         options->ignore_read_name_parsing_errors,
         "If enabled, parsing errors encountered when parsing read names be ignored and "
         "the read will be treated as if it is a full read with no UMIs. This is useful for allowing "
         "processing of BAM files that do not follow the expected read name format without having to preprocess "
         "the BAM to fix the read names. Use with caution as this can lead to unintended consequences if the "
         "input BAM does in fact have UMIs and/or SIDs in the read names but they are not in the expected format.")
      ->expected(0)
      ->default_val("false")
      ->group(kOptGroupNameHiddenOptions);
}

void AddCommonOutputOptions(CLI::App* const app,
                            const ReadCollapserOptionsPtr& options,
                            const std::string& group_name) {
  app->add_option("--output-dir", options->output_dir, "Directory for output file(s) and metrics.")
      ->default_val(".")
      ->group(group_name);
  app->add_option("--output-prefix",
                  options->output_prefix,
                  "Prefix prepended to every output filename, joined with an underscore "
                  "(e.g. 'SAMPLE' -> 'SAMPLE_<filename>'). Defaults to no prefix.")
      ->group(group_name);
  app->add_option<bool>("--overwrite", options->overwrite, "Overwrite existing output files.")
      ->expected(0)
      ->default_val("false")
      ->group(group_name);
}

void AddCommonReadFilteringOptions(CLI::App* const app,
                                   const ReadCollapserOptionsPtr& options,
                                   const std::string& group_name) {
  app->add_option("--min-mapping-quality", options->min_mapq, "Minimum mapping quality.")
      ->check(CLI::Range(0U, kU8Limit))
      ->default_val(0)
      ->group(group_name);
  app->add_option(
         "--max-discordant-duplex-error-percentage",
         options->max_discordant_duplex_error_percentage,
         "Maximum percentage of discordant duplex bases allowed in a hairpin duplex record before it is discarded.")
      ->check(CLI::Range(0.0, 100.0))
      ->group(group_name);
  app->add_option<bool>(
         "--exclude-partial-reads", options->exclude_partial_reads, "Exclude partial reads from analysis.")
      ->expected(0)
      ->default_val("false")
      ->group(group_name);
}

void AddConsensusReadFilterOptions(CLI::App* const app,
                                   const ReadCollapserOptionsPtr& options,
                                   const std::string& group_name) {
  app->add_option("--exclude-flags",
                  options->exclude_flags,
                  "Exclude flags to exclude BAM records from analysis. "
                  "(BAM_FSUPPLEMENTARY | BAM_FSECONDARY).  See "
                  "https://broadinstitute.github.io/picard/explain-flags.html for how to "
                  "generate flags.")
      ->default_val(kDefaultExcludeFlagConsensus)
      ->group(group_name);
}

void AddMarkdupReadFilterOptions(CLI::App* const app,
                                 const ReadCollapserOptionsPtr& options,
                                 const std::string& group_name) {
  app->add_option("--exclude-flags",
                  options->exclude_flags,
                  "Exclude flags to exclude BAM records from analysis. "
                  "(BAM_FSECONDARY).  See "
                  "https://broadinstitute.github.io/picard/explain-flags.html for how to "
                  "generate flags.")
      ->default_val(kDefaultExcludeFlagMarkdup)
      ->group(group_name);
}

void AddCommonClusterOptions(CLI::App* const app,
                             const ReadCollapserOptionsPtr& options,
                             const std::string& group_name) {
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
      ->check(CLI::Range(0U, kU8Limit))
      ->default_val(kDefaultWiggleRoom)
      ->group(group_name);
  app->add_option("--wiggle-room-partial",
                  options->wiggle_room_partial,
                  "Maximum distance between read alignments to be in the same cluster for partial reads.")
      ->check(CLI::Range(0U, kU8Limit))
      ->default_val(kDefaultWiggleRoomPartial)
      ->group(group_name);
}

void AddCommonPerformanceOptions(CLI::App* const app,
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

void AddConsensusOptions(CLI::App* const app, const ReadCollapserOptionsPtr& options, const std::string& group_name) {
  app->add_option("--min-trim-read-support",
                  options->consensus_options.min_trim_read_support,
                  "The minimum number of reads required for a specific base position at the ends of the consensus "
                  "sequence. Any bases at "
                  "the start or end of the consensus sequence with read support below this threshold will be trimmed. "
                  "Must be <= --min-cluster-size.")
      ->check(CLI::PositiveNumber)
      ->default_val(kMinReadsPerCluster)
      ->group(group_name);
  app->add_option("--min-same-strand-cluster-size",
                  options->consensus_options.min_same_strand_cluster_size,
                  "Consensus sequences generated from same-strand clusters will be discarded if the effective cluster "
                  "size is below this threshold. The effective cluster size is calculated by the total number of bases "
                  "in the consensus matrix divided by the number of covered positions. "
                  "Defaults to the value of --min-cluster-size if not specified.")
      ->check(CLI::PositiveNumber)
      ->group(group_name);
  app->add_option("--min-mixed-strand-cluster-size",
                  options->consensus_options.min_mixed_strand_cluster_size,
                  "Consensus sequences generated from mixed-strand clusters will be discarded if the effective cluster "
                  "size is below this threshold. The effective cluster size is calculated by the total number of bases "
                  "in the consensus matrix divided by the number of covered positions. "
                  "Defaults to the value of --min-cluster-size if not specified.")
      ->check(CLI::PositiveNumber)
      ->group(group_name);
  app->add_option("--min-consensus-fraction",
                  options->consensus_options.consensus_threshold,
                  "The fraction of reads that must agree on a base call for it to be considered a strong consensus.")
      ->check(CLI::Range(0.0, 1.0))
      ->default_val(kDefaultConsensusThreshold)
      ->group(group_name);
  app->add_option("--min-consensus-gap-fraction",
                  options->consensus_options.consensus_gap_threshold,
                  "The fraction of reads at a given position that must support a gap for a gap to be called in the "
                  "consensus sequence.")
      ->check(CLI::Range(0.0, 1.0))
      ->default_val(kDefaultConsensusGapThreshold)
      ->group(group_name);
  app->add_option<bool>("--skip-softclip-consensus",
                        options->consensus_options.skip_softclips,
                        "Skip soft-clip consensus generation. By default, soft-clipped segments are included in "
                        "consensus generation.")
      ->expected(0)
      ->default_val("false")
      ->group(group_name);
  app->add_option<bool>("--enable-legacy-qscore-model",
                        options->consensus_options.enable_legacy_qscore_model,
                        "Use legacy Q-score model for consensus quality scoring.")
      ->expected(0)
      ->default_val("false")
      ->group(group_name);
  cli::AddEnumOption(app,
                     "--duplex-library-type",
                     options->consensus_options.duplex_library_type,
                     "Specifies the method for decoding duplex reads into their constituent single-strand reads (R1 "
                     "and R2). This mode is used to inform downstream consensus calling by treating R1 and R2 as "
                     "separate entities from a single molecule. None implies a simplex library.",
                     HDDeconvolutionType::kNone)
      ->group(group_name);
  app->add_option("--min-consensus-read-length",
                  options->consensus_options.min_consensus_read_length,
                  "Discard consensus reads shorter than the specified length. If not specified, no consensus reads are "
                  "discarded based on length.")
      ->check(CLI::PositiveNumber)
      ->group(group_name);
}

void AddMarkdupOutputOptions(CLI::App* const app,
                             const ReadCollapserOptionsPtr& options,
                             const std::string& group_name) {
  app->add_option<bool>("--remove-duplicates",
                        options->markdup_options.remove_duplicates,
                        "Remove duplicate reads from the output BAM file(s).")
      ->expected(0)
      ->default_val("false")
      ->group(group_name);
  app->add_option<bool>("--exclude-cluster-tags",
                        options->markdup_options.exclude_cluster_tags,
                        "Excludes cluster ID and size from the output BAM file(s).")
      ->expected(0)
      ->group(group_name);
  app->add_option<bool>(
         "--merge-output", options->merge_output, "Merge output BAM files into a single position-sorted BAM file.")
      ->expected(0)
      ->default_val("false")
      ->group(group_name);
  app->add_option<bool>("--mark-supplementary-alignments",
                        options->markdup_options.mark_supplementary_alignments,
                        "Mark supplementary alignments of duplicate reads as duplicates. "
                        "Requires an additional read/write pass over the output, which increases runtime.")
      ->expected(0)
      ->default_val("false")
      ->group(group_name);
  app->add_option<bool>("--keep-pg-tag",
                        options->markdup_options.keep_pg_tag,
                        "Keep the PG auxiliary tag on output BAM records. By default, PG tags are removed to avoid "
                        "confusion when re-marking duplicates on BAMs produced by other tools.")
      ->expected(0)
      ->default_val("false")
      ->group(group_name);
}

void AddConsensusOutputOptions(CLI::App* const app,
                               const ReadCollapserOptionsPtr& options,
                               const std::string& group_name) {
  app->add_option("--compression-level",
                  options->consensus_options.compression_level,
                  "The level of compression used for the FASTQ output, 1-9.")
      ->check(CLI::Range(kMinCompressionLevel, kMaxCompressionLevel))
      ->default_val(kMinCompressionLevel)
      ->group(group_name);
  app->add_option<bool>(
         "--write-cluster-bam", options->consensus_options.output_cluster_bam, "Output cluster BAM file(s).")
      ->expected(0)
      ->default_val("false")
      ->group(group_name);
  app->add_option<bool>(
         "--include-consensus-aux-tags",
         options->consensus_options.include_consensus_aux_tags,
         "If enabled, adds duplex consensus auxiliary tags to the FASTQ output: per-base single-strand depth (ad/bd) "
         "and errors (ae/be) as arrays, plus per-strand and combined max depth (aD/bD/cD), min depth "
         "(aM/bM/cM), and error rate (aE/bE/cE).")
      ->expected(0)
      ->default_val("false")
      ->group(group_name);
}

void AddConsensusClusterOptions(CLI::App* const app,
                                const ReadCollapserOptionsPtr& options,
                                const std::string& group_name) {
  app->add_option("--max-cluster-size",
                  options->consensus_options.max_cluster_size,
                  "The maximum reads in a cluster before downsampling.")
      ->check(CLI::Range(1, QscoreCalculator::kMaxReadsInCluster))
      ->default_val(kMaxReadsPerCluster)
      ->group(group_name);
  app->add_option("--min-cluster-size",
                  options->consensus_options.min_cluster_size,
                  "The minimum reads in a cluster. Clusters with fewer reads will be discarded.")
      ->check(CLI::PositiveNumber)
      ->default_val(kMinReadsPerCluster)
      ->group(group_name);
}

namespace {

void WarnIfIgnoringReadNameParsingErrors(const ReadCollapserOptionsPtr& options) {
  if (options->ignore_read_name_parsing_errors) {
    Logging::Warn(
        "--ignore-read-name-parsing-errors is enabled. Errors encountered when parsing read names will be ignored "
        "and the read will be treated as if it has both SIDs present and no UMIs. This can lead to unintended "
        "consequences "
        "if the input BAM does in fact have UMIs and/or SIDs in the read names but they are not in the expected "
        "format.");
    if (options->cluster_by_umi) {
      Logging::Warn(
          "--cluster-by-umi is enabled in combination with --ignore-read-name-parsing-errors. "
          "If an error occurs when parsing a read name, the read will be treated as if it has no UMIs and discarded. "
          "This can lead to unintended consequences. UMI clustering is not possible if the read name does not follow "
          "the expected format. "
          "If you wish to process these reads anyway, you can disable UMI clustering and the reads will be clustered "
          "based on alignment positions alone.");
    }
  }
}

}  // namespace

void ValidateConsensusOptions(const CLI::App* const app, const ReadCollapserOptionsPtr& options) {
  // If not explicitly provided (by user or preset), default strand cluster sizes to the value of --min-cluster-size
  if (options->consensus_options.min_same_strand_cluster_size == 0) {
    options->consensus_options.min_same_strand_cluster_size = options->consensus_options.min_cluster_size;
  }
  if (options->consensus_options.min_mixed_strand_cluster_size == 0) {
    options->consensus_options.min_mixed_strand_cluster_size = options->consensus_options.min_cluster_size;
  }
  // `--min-trim-read-support` must be <= `--min-cluster-size`
  if (options->consensus_options.min_trim_read_support > options->consensus_options.min_cluster_size) {
    throw CLI::ValidationError("--min-trim-read-support must be <= --min-cluster-size.");
  }
  // Manually set exclude_flags to default value for consensus to avoid it
  // being overwritten by markdup default value
  if (app->get_option("--exclude-flags")->count() == 0) {
    options->exclude_flags = kDefaultExcludeFlagConsensus;
  }
  // Print warning if the user specifies `--ignore-read-name-parsing-errors` since this can lead to unintended
  // consequences if the input BAM does in fact have UMIs and/or SIDs in the read names but they are not in the expected
  // format
  WarnIfIgnoringReadNameParsingErrors(options);
}

void ValidateMarkdupOptions(const CLI::App* const app, const ReadCollapserOptionsPtr& options) {
  // Manually set exclude_flags to default value for markdup to avoid it
  // being overwritten by consensus default value
  if (app->get_option("--exclude-flags")->count() == 0) {
    options->exclude_flags = kDefaultExcludeFlagMarkdup;
  }
  options->markdup_options.cluster_rescued_secondaries = true;
  // Print warning if the user specifies `--ignore-read-name-parsing-errors` since this can lead to unintended
  // consequences if the input BAM does in fact have UMIs and/or SIDs in the read names but they are not in the expected
  // format
  WarnIfIgnoringReadNameParsingErrors(options);
}

}  // namespace xoos::read_collapser
