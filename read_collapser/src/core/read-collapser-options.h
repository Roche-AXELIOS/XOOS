#pragma once

#include <memory>
#include <optional>
#include <string>

#include <htslib/sam.h>

#include <CLI/Formatter.hpp>

#include <xoos/cli/cli.h>
#include <xoos/io/metadata-util.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::read_collapser {

enum class HDDeconvolutionType : u8 {
  kParentParent,
  kParentDaughter,
  kNone,
};

enum class ReadCollapserPresets : u8 {
  kWgsDuplex,
  kWgsDuplexCfdna,
  kWgsSimplex,
  kTeDuplex,
  kTeSimplex,
  kRnaBulk,
  kNone,
};

// This enum defines the mode of operation for read collapser, specifically whether it is handling duplicate marking or
// consensus. The set of output files and the metric output both differ between modes (e.g. consensus metrics will be
// NA for duplicate marking mode).
enum class ReadCollapserMode : u8 {
  kMarkDuplicate,
  kConsensus
};

// Options specific to the `markdup` subcommand.
struct MarkdupOptions {
  bool cluster_rescued_secondaries{false};
  bool remove_duplicates{false};
  bool exclude_cluster_tags{false};
  bool mark_supplementary_alignments{false};
  bool keep_pg_tag{false};
};

// Options specific to the `consensus` subcommand.
struct ConsensusOptions {
  u8 compression_level{};
  bool output_cluster_bam{false};
  u32 max_cluster_size{};
  u8 min_cluster_size{};
  u8 min_trim_read_support{};
  u8 min_same_strand_cluster_size{};
  u8 min_mixed_strand_cluster_size{};
  f64 consensus_threshold{};
  f64 consensus_gap_threshold{};
  bool skip_softclips{false};
  bool enable_legacy_qscore_model{false};
  HDDeconvolutionType duplex_library_type{};
  std::optional<u32> min_consensus_read_length;

  // When enabled, emit per-consensus-read duplex auxiliary tags in the FASTQ output: per-base single-strand depth/error
  // arrays and per-strand and combined depth/error-rate summaries. This lets downstream pipelines built for that tag
  // schema consume Read Collapser output.
  bool include_consensus_aux_tags{false};
};

struct ReadCollapserOptions {
  // Common options
  fs::path bam_input;
  std::optional<fs::path> bed_input;
  s32 padding{};

  fs::path output_dir;
  // Prepended to every output filename, joined with an underscore, e.g. prefix "SAMPLE" ->
  // "SAMPLE_output.bam". Empty (the default) leaves filenames unprefixed.
  std::string output_prefix;
  bool overwrite{false};
  bool merge_output{false};

  u8 min_mapq{};
  std::optional<f64> max_discordant_duplex_error_percentage;
  u16 exclude_flags{};
  bool exclude_partial_reads{false};

  bool cluster_by_umi{false};
  bool cluster_by_strand{false};
  bool make_clusters_of_partial_reads_only{false};
  u8 wiggle_room{};
  u8 wiggle_room_partial{};
  bool ignore_read_name_parsing_errors{false};

  size_t threads{};
  u32 region_size{};
  u32 batch_size{};

  // Mode-specific options
  MarkdupOptions markdup_options;
  ConsensusOptions consensus_options;

  // Program information
  io::CommandLineInfo command_line_info{};
};

using ReadCollapserOptionsPtr = std::shared_ptr<ReadCollapserOptions>;

}  // namespace xoos::read_collapser
