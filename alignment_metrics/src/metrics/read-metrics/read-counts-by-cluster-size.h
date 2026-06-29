#pragma once
#include <string>
#include <vector>

#include <xoos/io/metadata-util.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

namespace xoos::alignment_metrics {

/**
 * @brief Read counts for a single cluster size stratified by cluster characteristics.
 *
 * Tracks the number of reads belonging to a cluster of a specific size, broken down by
 * strand composition (mixed, forward-only, reverse-only) and read completeness (full,
 * partial, or mixed).
 */
struct ReadCountByClusterSize {
  u64 mixed_strand_cluster_reads_passing_filter{};
  u64 forward_cluster_reads_passing_filter{};
  u64 reverse_cluster_reads_passing_filter{};
  u64 full_cluster_reads_passing_filter{};
  u64 mixed_full_and_partial_cluster_reads_passing_filter{};
  u64 partial_cluster_reads_passing_filter{};
  u32 cluster_size{};

  bool operator==(const ReadCountByClusterSize& other) const = default;
};

/**
 * @brief Container for read counts across a range of cluster sizes.
 *
 * Aggregates ReadCountByClusterSize statistics for cluster sizes from 1 up to max_cluster_size_bin.
 */
struct ReadCountsByClusterSize {
  ReadCountsByClusterSize() = default;
  explicit ReadCountsByClusterSize(u32 max_cluster_size_bin);

  u32 max_cluster_size_bin{};
  std::vector<ReadCountByClusterSize> read_counts_by_cluster_size{};

  static std::vector<std::string> GetHeaders();
  void WriteTsv(const fs::path& output_path, const io::Comments& comments) const;

  ReadCountsByClusterSize& operator+=(const ReadCountsByClusterSize& other);
  ReadCountByClusterSize& operator[](size_t cluster_size_index);
};
}  // namespace xoos::alignment_metrics
