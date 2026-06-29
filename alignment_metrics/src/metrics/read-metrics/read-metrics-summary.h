#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include <xoos/types/int.h>
// this is needed for compilation in consensus-accuracy-metrics, but flagged as not needed
#include <xoos/io/metadata-util.h>
#include <xoos/types/fs.h>
#include <xoos/types/optional.h>  // pragma: keep

#include "metadata/dataset-metadata.h"

namespace xoos::alignment_metrics {

/**
 * @brief Comprehensive summary statistics for read alignments and sequencing quality.
 *
 * Aggregates high-level metrics about reads, alignments, and bases across the entire dataset.
 * Tracks mapping statistics (mapped/unmapped, forward/reverse), alignment flags (secondary,
 * supplementary, duplicate), quality filtering results, and base composition. For post-consensus
 * datasets with cluster information, also tracks read stratification by cluster characteristics
 * (strand composition, read completeness).
 */
struct ReadMetricsSummary {
  u64 alignments_plus_unmapped_reads{};
  u64 total_reads{};
  u64 unmapped_reads{};
  u64 mapped_reads{};
  u64 forward_reads{};
  u64 reverse_reads{};
  std::optional<u64> full_length_reads;
  std::optional<u64> partial_length_reads;
  std::optional<u64> mixed_strand_cluster_reads;
  std::optional<u64> forward_cluster_reads;
  std::optional<u64> reverse_cluster_reads;
  std::optional<u64> full_cluster_reads;
  std::optional<u64> mixed_full_and_partial_cluster_reads;
  std::optional<u64> partial_cluster_reads;
  u64 secondary_alignments{};
  u64 supplementary_alignments{};
  u64 duplicate_reads{};
  u64 mapq_zero_reads{};
  u64 reads_passing_filter{};
  u64 forward_reads_passing_filter{};
  u64 reverse_reads_passing_filter{};
  std::optional<u64> full_length_reads_passing_filter;
  std::optional<u64> partial_length_reads_passing_filter;
  std::optional<u64> mixed_strand_cluster_reads_passing_filter;
  std::optional<u64> forward_cluster_reads_passing_filter;
  std::optional<u64> reverse_cluster_reads_passing_filter;
  std::optional<u64> full_cluster_reads_passing_filter;
  std::optional<u64> mixed_full_and_partial_cluster_reads_passing_filter;
  std::optional<u64> partial_cluster_reads_passing_filter;
  u64 total_bases{};
  u64 unmapped_bases{};
  u64 aligned_bases{};
  u64 soft_clipped_bases{};
  u64 gc_bases{};
  u64 total_bases_passing_filter{};
  u64 aligned_bases_passing_filter{};
  u64 soft_clipped_bases_passing_filter{};
  u64 gc_bases_passing_filter{};

  ReadMetricsSummary() = default;
  explicit ReadMetricsSummary(const DatasetMetadata& dataset_metadata);

  ReadMetricsSummary& operator+=(const ReadMetricsSummary& obj);
  static std::vector<std::string> GetHeaders();
  void WriteTsv(const fs::path& output_path, const io::Comments& comments) const;
};

}  // namespace xoos::alignment_metrics
