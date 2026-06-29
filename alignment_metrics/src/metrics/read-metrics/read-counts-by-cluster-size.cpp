#include "metrics/read-metrics/read-counts-by-cluster-size.h"

#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>

#include <csv.hpp>

#include <xoos/types/int.h>

#include "metrics/metrics-names.h"

namespace xoos::alignment_metrics {
ReadCountsByClusterSize::ReadCountsByClusterSize(const u32 max_cluster_size_bin)
    : max_cluster_size_bin(max_cluster_size_bin) {
  // +1 to for max+ indexes
  read_counts_by_cluster_size.resize(max_cluster_size_bin + 1);

  for (u32 i = 0; i < max_cluster_size_bin + 1; ++i) {
    // we store `cluster_size` at index `cluster_size - 1` since cluster_size starts from 1
    read_counts_by_cluster_size[i].cluster_size = i + 1;
  }
}

std::vector<std::string> ReadCountsByClusterSize::GetHeaders() {
  return {kNameClusterSize,
          kNameMixedStrandClusterReadsPassingFilter,
          kNameForwardClusterReadsPassingFilter,
          kNameReverseClusterReadsPassingFilter,
          kNameFullClusterReadsPassingFilter,
          kNameMixedFullAndPartialClusterReadsPassingFilter,
          kNamePartialClusterReadsPassingFilter};
}

// Writes TSV file for the errors by pass for each cluster size
void ReadCountsByClusterSize::WriteTsv(const fs::path& output_path, const io::Comments& comments) const {
  std::ofstream ofstream(output_path);
  auto writer = csv::make_tsv_writer_buffered(ofstream);
  io::WriteTsvComments(writer, comments);
  writer << GetHeaders();
  for (const auto& read_count_by_cluster_size : read_counts_by_cluster_size) {
    // if the cluster size is 0 then we skip it
    if (read_count_by_cluster_size.cluster_size == 0) {
      continue;
    }
    auto cluster_size = read_count_by_cluster_size.cluster_size <= max_cluster_size_bin
                            ? std::to_string(read_count_by_cluster_size.cluster_size)
                            : std::to_string(max_cluster_size_bin) + "+";

    // write the data to the tsv file
    writer << std::make_tuple(cluster_size,
                              read_count_by_cluster_size.mixed_strand_cluster_reads_passing_filter,
                              read_count_by_cluster_size.forward_cluster_reads_passing_filter,
                              read_count_by_cluster_size.reverse_cluster_reads_passing_filter,
                              read_count_by_cluster_size.full_cluster_reads_passing_filter,
                              read_count_by_cluster_size.mixed_full_and_partial_cluster_reads_passing_filter,
                              read_count_by_cluster_size.partial_cluster_reads_passing_filter);
  }
}

ReadCountsByClusterSize& ReadCountsByClusterSize::operator+=(const ReadCountsByClusterSize& other) {
  for (const auto& value : other.read_counts_by_cluster_size) {
    auto cluster_size = value.cluster_size;
    // if cluster_size is more than the max_cluster_size_bin then store it into max_cluster_size_bin+
    if (value.cluster_size > max_cluster_size_bin) {
      cluster_size = max_cluster_size_bin + 1;
    }
    // we store `cluster_size` at index `cluster_size - 1` since cluster_size starts from 1
    auto index = cluster_size - 1;
    auto& read_count_by_cluster_size = read_counts_by_cluster_size[index];

    read_count_by_cluster_size.cluster_size = value.cluster_size;
    read_count_by_cluster_size.mixed_strand_cluster_reads_passing_filter +=
        value.mixed_strand_cluster_reads_passing_filter;
    read_count_by_cluster_size.forward_cluster_reads_passing_filter += value.forward_cluster_reads_passing_filter;
    read_count_by_cluster_size.reverse_cluster_reads_passing_filter += value.reverse_cluster_reads_passing_filter;
    read_count_by_cluster_size.full_cluster_reads_passing_filter += value.full_cluster_reads_passing_filter;
    read_count_by_cluster_size.mixed_full_and_partial_cluster_reads_passing_filter +=
        value.mixed_full_and_partial_cluster_reads_passing_filter;
    read_count_by_cluster_size.partial_cluster_reads_passing_filter += value.partial_cluster_reads_passing_filter;
  }
  return *this;
}

ReadCountByClusterSize& ReadCountsByClusterSize::operator[](size_t cluster_size_index) {
  return read_counts_by_cluster_size[cluster_size_index];
}
}  // namespace xoos::alignment_metrics
