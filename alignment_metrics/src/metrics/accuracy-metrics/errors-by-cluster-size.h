#pragma once
#include <string>
#include <vector>

#include <ankerl/unordered_dense.h>

#include <xoos/io/metadata-util.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "core/error-counts.h"

namespace xoos::alignment_metrics {

// Error counts for a single cluster size
struct ErrorByClusterSize {
  u64 total_bases_mixed_strand_cluster{};
  u64 total_bases_forward_cluster{};
  u64 total_bases_reverse_cluster{};
  u64 total_bases_full_cluster{};
  u64 total_bases_mixed_full_and_partial_cluster{};
  u64 total_bases_partial_cluster{};
  u64 substitutions_mixed_strand_cluster{};
  u64 substitutions_forward_cluster{};
  u64 substitutions_reverse_cluster{};
  u64 substitutions_full_cluster{};
  u64 substitutions_mixed_full_and_partial_cluster{};
  u64 substitutions_partial_cluster{};
  u64 insertions_mixed_strand_cluster{};
  u64 insertions_forward_cluster{};
  u64 insertions_reverse_cluster{};
  u64 insertions_full_cluster{};
  u64 insertions_mixed_full_and_partial_cluster{};
  u64 insertions_partial_cluster{};
  u64 deletions_mixed_strand_cluster{};
  u64 deletions_forward_cluster{};
  u64 deletions_reverse_cluster{};
  u64 deletions_full_cluster{};
  u64 deletions_mixed_full_and_partial_cluster{};
  u64 deletions_partial_cluster{};
  u32 cluster_size{};
  bool Populated() const;

  bool operator==(const ErrorByClusterSize& other) const = default;
};

// ErrorsByClusterSize is a container for ErrorByClusterSize for a range of cluster sizes up to max_cluster_size_bin
struct ErrorsByClusterSize {
  ErrorsByClusterSize() = default;
  explicit ErrorsByClusterSize(u32 max_cluster_size_bin);
  // key is cluster size
  u32 max_cluster_size_bin{};
  std::vector<ErrorByClusterSize> errors_by_cluster_size{};

  static std::vector<std::string> GetHeaders();
  void WriteTsv(const fs::path& output_path, const io::Comments& comments) const;

  ErrorsByClusterSize& operator+=(const ErrorsByClusterSize& other);
  ErrorsByClusterSize& operator+=(
      const std::vector<CountsByErrorTypeAndReadProperty<u32>>& combined_error_stats_cluster_size);
};
}  // namespace xoos::alignment_metrics
