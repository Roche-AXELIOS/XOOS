#pragma once

#include <memory>
#include <optional>
#include <string>

#include <htslib/sam.h>

#include <CLI/Formatter.hpp>

#include <xoos/cli/cli.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::read_collapser {

enum class HDDeconvolutionType {
  kParentParent,
  kParentDaughter,
  kNone,
};

enum class ReadCollapserPresets {
  kWgsDuplex,
  kWgsDuplexMrd,
  kWgsDuplexCfdna,
  kWgsSimplex,
  kTeDuplex,
  kTeSimplex,
  kRnaBulk,
  kNone,
};

struct ReadCollapserOptions {
  // Common options
  fs::path bam_input{};
  std::optional<fs::path> bed_input{};
  s32 padding{};

  fs::path output_dir{};
  bool merge_output{false};

  u8 min_mapq{};
  std::optional<f64> max_discordant_duplex_error_percentage{};
  u16 exclude_flags{};
  bool exclude_partial_reads{false};

  bool cluster_by_umi{false};
  bool cluster_by_strand{false};
  bool make_clusters_of_partial_reads_only{false};
  u8 wiggle_room{};
  u8 wiggle_room_partial{};

  size_t threads{};
  u32 region_size{};
  u32 batch_size{};

  // Markdup options
  bool remove_duplicates{false};
  bool exclude_cluster_tags{false};

  // Consensus options
  u8 compression_level{};
  bool output_cluster_bam{false};
  u32 max_cluster_size{};
  u8 min_cluster_size{};
  u8 min_trim_read_support{};
  u8 min_same_strand_cluster_size{};
  u8 min_mixed_strand_cluster_size{};
  f64 consensus_threshold{};
  f64 consensus_gap_threshold{};
  bool include_softclips{false};
  bool enable_legacy_qscore_model{false};
  HDDeconvolutionType duplex_library_type{};
  std::optional<u32> min_consensus_read_length{};

  // Consensus debug options
  bool include_per_base_read_support_tags{false};
  bool include_per_base_majority_count_tags{false};

  // Program information
  std::string command_line{};
  std::string program_name{};
  std::string version{};
};

using ReadCollapserOptionsPtr = std::shared_ptr<ReadCollapserOptions>;

}  // namespace xoos::read_collapser
