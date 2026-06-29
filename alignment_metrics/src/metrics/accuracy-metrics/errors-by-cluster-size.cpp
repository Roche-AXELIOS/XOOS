#include "metrics/accuracy-metrics/errors-by-cluster-size.h"

#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>

#include <ankerl/unordered_dense.h>

#include <csv.hpp>

#include <xoos/types/int.h>
#include <xoos/util/math.h>

#include "core/error-counts.h"
#include "metrics/metrics-names.h"
#include "util/format-util.h"

namespace xoos::alignment_metrics {
ErrorsByClusterSize::ErrorsByClusterSize(const u32 max_cluster_size_bin) : max_cluster_size_bin(max_cluster_size_bin) {
  // +1 to for max+ indexes
  errors_by_cluster_size.resize(max_cluster_size_bin + 1);

  for (u32 i = 0; i < max_cluster_size_bin + 1; i++) {
    // we store `cluster_size` at index `cluster_size - 1` since cluster_size starts from 1
    errors_by_cluster_size[i].cluster_size = i + 1;
  }
}

std::vector<std::string> ErrorsByClusterSize::GetHeaders() {
  return {kNameClusterSize,
          kNameTotalBases,
          kNameOverallPhred,
          kNameTotalBasesMixedStrandCluster,
          kNameTotalBasesForwardCluster,
          kNameTotalBasesReverseCluster,
          kNameTotalBasesFullCluster,
          kNameTotalBasesMixedFullAndPartialCluster,
          kNameTotalBasesPartialCluster,
          kNameSubstitutionsTotal,
          kNameSubstitutionsPhred,
          kNameSubstitutionsMixedStrandCluster,
          kNameSubstitutionsForwardCluster,
          kNameSubstitutionsReverseCluster,
          kNameSubstitutionsFullCluster,
          kNameSubstitutionsMixedFullAndPartialCluster,
          kNameSubstitutionsPartialCluster,
          kNameInsertionsTotal,
          kNameInsertionsPhred,
          kNameInsertionsMixedStrandCluster,
          kNameInsertionsForwardCluster,
          kNameInsertionsReverseCluster,
          kNameInsertionsFullCluster,
          kNameInsertionsMixedFullAndPartialCluster,
          kNameInsertionsPartialCluster,
          kNameDeletionsTotal,
          kNameDeletionsPhred,
          kNameDeletionsMixedStrandCluster,
          kNameDeletionsForwardCluster,
          kNameDeletionsReverseCluster,
          kNameDeletionsFullCluster,
          kNameDeletionsMixedFullAndPartialCluster,
          kNameDeletionsPartialCluster};
}

// Writes TSV file for the errors by pass for each cluster size
void ErrorsByClusterSize::WriteTsv(const fs::path& output_path, const io::Comments& comments) const {
  std::ofstream ofstream(output_path);
  auto writer = csv::make_tsv_writer_buffered(ofstream);
  io::WriteTsvComments(writer, comments);
  writer << GetHeaders();
  for (const auto& error_by_cluster_size : errors_by_cluster_size) {
    // if the cluster size is 0 then we skip it
    if (error_by_cluster_size.cluster_size == 0) {
      continue;
    }
    auto cluster_size = error_by_cluster_size.cluster_size <= max_cluster_size_bin
                            ? std::to_string(error_by_cluster_size.cluster_size)
                            : std::to_string(max_cluster_size_bin) + "+";

    auto total_bases = error_by_cluster_size.total_bases_mixed_strand_cluster +
                       error_by_cluster_size.total_bases_forward_cluster +
                       error_by_cluster_size.total_bases_reverse_cluster;
    auto substitutions_total = error_by_cluster_size.substitutions_mixed_strand_cluster +
                               error_by_cluster_size.substitutions_forward_cluster +
                               error_by_cluster_size.substitutions_reverse_cluster;
    auto insertions_total = error_by_cluster_size.insertions_mixed_strand_cluster +
                            error_by_cluster_size.insertions_forward_cluster +
                            error_by_cluster_size.insertions_reverse_cluster;
    auto deletions_total = error_by_cluster_size.deletions_mixed_strand_cluster +
                           error_by_cluster_size.deletions_forward_cluster +
                           error_by_cluster_size.deletions_reverse_cluster;
    auto overall_phred =
        total_bases > 0
            ? ToStringWithPrecision(
                  math::ErrorRateToPhred(static_cast<f64>(substitutions_total + insertions_total + deletions_total) /
                                         static_cast<f64>(total_bases)),
                  2)
            : kNotApplicable;
    auto substitutions_phred =
        total_bases > 0
            ? ToStringWithPrecision(
                  math::ErrorRateToPhred(static_cast<f64>(substitutions_total) / static_cast<f64>(total_bases)), 2)
            : kNotApplicable;
    auto insertions_phred =
        total_bases > 0
            ? ToStringWithPrecision(
                  math::ErrorRateToPhred(static_cast<f64>(insertions_total) / static_cast<f64>(total_bases)), 2)
            : kNotApplicable;
    auto deletions_phred =
        total_bases > 0
            ? ToStringWithPrecision(
                  math::ErrorRateToPhred(static_cast<f64>(deletions_total) / static_cast<f64>(total_bases)), 2)
            : kNotApplicable;

    // write the data to the tsv file
    writer << std::make_tuple(cluster_size,
                              total_bases,
                              overall_phred,
                              error_by_cluster_size.total_bases_mixed_strand_cluster,
                              error_by_cluster_size.total_bases_forward_cluster,
                              error_by_cluster_size.total_bases_reverse_cluster,
                              error_by_cluster_size.total_bases_full_cluster,
                              error_by_cluster_size.total_bases_mixed_full_and_partial_cluster,
                              error_by_cluster_size.total_bases_partial_cluster,
                              substitutions_total,
                              substitutions_phred,
                              error_by_cluster_size.substitutions_mixed_strand_cluster,
                              error_by_cluster_size.substitutions_forward_cluster,
                              error_by_cluster_size.substitutions_reverse_cluster,
                              error_by_cluster_size.substitutions_full_cluster,
                              error_by_cluster_size.substitutions_mixed_full_and_partial_cluster,
                              error_by_cluster_size.substitutions_partial_cluster,
                              insertions_total,
                              insertions_phred,
                              error_by_cluster_size.insertions_mixed_strand_cluster,
                              error_by_cluster_size.insertions_forward_cluster,
                              error_by_cluster_size.insertions_reverse_cluster,
                              error_by_cluster_size.insertions_full_cluster,
                              error_by_cluster_size.insertions_mixed_full_and_partial_cluster,
                              error_by_cluster_size.insertions_partial_cluster,
                              deletions_total,
                              deletions_phred,
                              error_by_cluster_size.deletions_mixed_strand_cluster,
                              error_by_cluster_size.deletions_forward_cluster,
                              error_by_cluster_size.deletions_reverse_cluster,
                              error_by_cluster_size.deletions_full_cluster,
                              error_by_cluster_size.deletions_mixed_full_and_partial_cluster,
                              error_by_cluster_size.deletions_partial_cluster);
  }
}

ErrorsByClusterSize& ErrorsByClusterSize::operator+=(const ErrorsByClusterSize& other) {
  for (const auto& value : other.errors_by_cluster_size) {
    auto cluster_size = value.cluster_size;
    // if cluster_size is more than the max_cluster_size_bin then store it into max_cluster_size_bin+
    if (value.cluster_size > max_cluster_size_bin) {
      cluster_size = max_cluster_size_bin + 1;
    }
    // we store `cluster_size` at index `cluster_size - 1` since cluster_size starts from 1
    const auto index = cluster_size - 1;
    auto& error_by_cluster_size = errors_by_cluster_size.at(index);

    error_by_cluster_size.cluster_size = value.cluster_size;
    error_by_cluster_size.total_bases_mixed_strand_cluster += value.total_bases_mixed_strand_cluster;
    error_by_cluster_size.total_bases_forward_cluster += value.total_bases_forward_cluster;
    error_by_cluster_size.total_bases_reverse_cluster += value.total_bases_reverse_cluster;
    error_by_cluster_size.total_bases_full_cluster += value.total_bases_full_cluster;
    error_by_cluster_size.total_bases_mixed_full_and_partial_cluster +=
        value.total_bases_mixed_full_and_partial_cluster;
    error_by_cluster_size.total_bases_partial_cluster += value.total_bases_partial_cluster;

    error_by_cluster_size.substitutions_mixed_strand_cluster += value.substitutions_mixed_strand_cluster;
    error_by_cluster_size.substitutions_forward_cluster += value.substitutions_forward_cluster;
    error_by_cluster_size.substitutions_reverse_cluster += value.substitutions_reverse_cluster;
    error_by_cluster_size.substitutions_full_cluster += value.substitutions_full_cluster;
    error_by_cluster_size.substitutions_mixed_full_and_partial_cluster +=
        value.substitutions_mixed_full_and_partial_cluster;
    error_by_cluster_size.substitutions_partial_cluster += value.substitutions_partial_cluster;

    error_by_cluster_size.insertions_mixed_strand_cluster += value.insertions_mixed_strand_cluster;
    error_by_cluster_size.insertions_forward_cluster += value.insertions_forward_cluster;
    error_by_cluster_size.insertions_reverse_cluster += value.insertions_reverse_cluster;
    error_by_cluster_size.insertions_full_cluster += value.insertions_full_cluster;
    error_by_cluster_size.insertions_mixed_full_and_partial_cluster += value.insertions_mixed_full_and_partial_cluster;
    error_by_cluster_size.insertions_partial_cluster += value.insertions_partial_cluster;

    error_by_cluster_size.deletions_mixed_strand_cluster += value.deletions_mixed_strand_cluster;
    error_by_cluster_size.deletions_forward_cluster += value.deletions_forward_cluster;
    error_by_cluster_size.deletions_reverse_cluster += value.deletions_reverse_cluster;
    error_by_cluster_size.deletions_full_cluster += value.deletions_full_cluster;
    error_by_cluster_size.deletions_mixed_full_and_partial_cluster += value.deletions_mixed_full_and_partial_cluster;
    error_by_cluster_size.deletions_partial_cluster += value.deletions_partial_cluster;
  }
  return *this;
}

ErrorsByClusterSize& ErrorsByClusterSize::operator+=(
    const std::vector<CountsByErrorTypeAndReadProperty<u32>>& combined_error_stats_cluster_size) {
  for (size_t cluster_size_index = 0; cluster_size_index < combined_error_stats_cluster_size.size();
       ++cluster_size_index) {
    const auto& error_counts = combined_error_stats_cluster_size[cluster_size_index];
    const u16 cluster_size_index_capped =
        cluster_size_index >= errors_by_cluster_size.size() ? errors_by_cluster_size.size() - 1 : cluster_size_index;
    auto& error_by_cluster_size = errors_by_cluster_size[cluster_size_index_capped];
    // Add total aligned base count for the current cluster size
    const auto& total_for_cluster_size = error_counts.totals;
    error_by_cluster_size.total_bases_mixed_strand_cluster += total_for_cluster_size.mixed_strand_cluster;
    error_by_cluster_size.total_bases_forward_cluster += total_for_cluster_size.forward_cluster;
    error_by_cluster_size.total_bases_reverse_cluster += total_for_cluster_size.reverse_cluster;
    error_by_cluster_size.total_bases_full_cluster += total_for_cluster_size.full_cluster;
    error_by_cluster_size.total_bases_mixed_full_and_partial_cluster +=
        total_for_cluster_size.mixed_full_and_partial_cluster;
    error_by_cluster_size.total_bases_partial_cluster += total_for_cluster_size.partial_cluster;

    // Add substitutions for the current cluster size
    const auto& substitutions_for_cluster_size = error_counts.substitutions;
    error_by_cluster_size.substitutions_mixed_strand_cluster += substitutions_for_cluster_size.mixed_strand_cluster;
    error_by_cluster_size.substitutions_forward_cluster += substitutions_for_cluster_size.forward_cluster;
    error_by_cluster_size.substitutions_reverse_cluster += substitutions_for_cluster_size.reverse_cluster;
    error_by_cluster_size.substitutions_full_cluster += substitutions_for_cluster_size.full_cluster;
    error_by_cluster_size.substitutions_mixed_full_and_partial_cluster +=
        substitutions_for_cluster_size.mixed_full_and_partial_cluster;
    error_by_cluster_size.substitutions_partial_cluster += substitutions_for_cluster_size.partial_cluster;

    // Add insertions for the current cluster size
    const auto& insertions_for_cluster_size = error_counts.insertions;
    error_by_cluster_size.insertions_mixed_strand_cluster += insertions_for_cluster_size.mixed_strand_cluster;
    error_by_cluster_size.insertions_forward_cluster += insertions_for_cluster_size.forward_cluster;
    error_by_cluster_size.insertions_reverse_cluster += insertions_for_cluster_size.reverse_cluster;
    error_by_cluster_size.insertions_full_cluster += insertions_for_cluster_size.full_cluster;
    error_by_cluster_size.insertions_mixed_full_and_partial_cluster +=
        insertions_for_cluster_size.mixed_full_and_partial_cluster;
    error_by_cluster_size.insertions_partial_cluster += insertions_for_cluster_size.partial_cluster;

    // Add deletions for the current cluster size
    const auto& deletions_for_cluster_size = error_counts.deletions;
    error_by_cluster_size.deletions_mixed_strand_cluster += deletions_for_cluster_size.mixed_strand_cluster;
    error_by_cluster_size.deletions_forward_cluster += deletions_for_cluster_size.forward_cluster;
    error_by_cluster_size.deletions_reverse_cluster += deletions_for_cluster_size.reverse_cluster;
    error_by_cluster_size.deletions_full_cluster += deletions_for_cluster_size.full_cluster;
    error_by_cluster_size.deletions_mixed_full_and_partial_cluster +=
        deletions_for_cluster_size.mixed_full_and_partial_cluster;
    error_by_cluster_size.deletions_partial_cluster += deletions_for_cluster_size.partial_cluster;
  }

  return *this;
}

// Returns true if any of the base or deletion counts is greater than 0 to indicate that the error by pass is populated
// with a single base in the reference
bool ErrorByClusterSize::Populated() const {
  return total_bases_forward_cluster + total_bases_reverse_cluster + total_bases_mixed_strand_cluster > 0;
}
}  // namespace xoos::alignment_metrics
