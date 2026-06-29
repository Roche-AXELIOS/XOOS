#include "metrics/accuracy-metrics/errors-by-substitution-type.h"

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <csv.hpp>

#include "core/error-counts.h"
#include "core/substitution-lookup.h"
#include "metrics/metrics-names.h"

namespace xoos::alignment_metrics {

ErrorsBySubstitutionType::ErrorsBySubstitutionType(const DatasetMetadata& dataset_metadata)
    : dataset_metadata(dataset_metadata) {
  // Initialize all members of the substitutions Index
  errors_by_substitution_type.resize(SubstitutionLookup::kMaxSubstitutionIndex);

  // For each substitution index, we set the ref and base values
  for (size_t index = 0; index < errors_by_substitution_type.size(); ++index) {
    const auto [ref, base] = SubstitutionLookup::GetSubstitutionForIndex(index);
    errors_by_substitution_type[index].ref = ref;
    errors_by_substitution_type[index].base = base;
  }
}

std::vector<std::string> ErrorsBySubstitutionType::GetHeaders() const {
  std::vector columns = {kNameSubstitutionType, kNameSubstitutionsTotal};
  if (dataset_metadata.has_cluster_info) {
    columns.emplace_back(kNameSubstitutionsMixedStrandCluster);
    columns.emplace_back(kNameSubstitutionsForwardCluster);
    columns.emplace_back(kNameSubstitutionsReverseCluster);
    columns.emplace_back(kNameSubstitutionsFullCluster);
    columns.emplace_back(kNameSubstitutionsMixedFullAndPartialCluster);
    columns.emplace_back(kNameSubstitutionsPartialCluster);
  }
  if (dataset_metadata.has_read_type_info) {
    columns.emplace_back(kNameSubstitutionsFullRead);
    columns.emplace_back(kNameSubstitutionsPartialRead);
  }
  if (dataset_metadata.has_strand_info) {
    columns.emplace_back(kNameSubstitutionsForwardStrand);
    columns.emplace_back(kNameSubstitutionsReverseStrand);
  }
  return columns;
}

void ErrorsBySubstitutionType::WriteTsv(const fs::path& output_path, const io::Comments& comments) const {
  std::ofstream ofstream(output_path);
  auto writer = csv::make_tsv_writer_buffered(ofstream);
  io::WriteTsvComments(writer, comments);
  writer << GetHeaders();
  // We write the total first
  std::vector total_row = {kNameSubstitutionTypeAny, std::to_string(substitutions_total)};
  // Only output columns that are applicable to the dataset
  if (dataset_metadata.has_cluster_info) {
    total_row.emplace_back(std::to_string(substitutions_mixed_strand_cluster));
    total_row.emplace_back(std::to_string(substitutions_forward_cluster));
    total_row.emplace_back(std::to_string(substitutions_reverse_cluster));
    total_row.emplace_back(std::to_string(substitutions_full_cluster));
    total_row.emplace_back(std::to_string(substitutions_mixed_full_and_partial_cluster));
    total_row.emplace_back(std::to_string(substitutions_partial_cluster));
  }
  if (dataset_metadata.has_read_type_info) {
    total_row.emplace_back(std::to_string(substitutions_full_read));
    total_row.emplace_back(std::to_string(substitutions_partial_read));
  }
  if (dataset_metadata.has_strand_info) {
    total_row.emplace_back(std::to_string(substitutions_forward_read));
    total_row.emplace_back(std::to_string(substitutions_reverse_read));
  }
  writer << total_row;
  for (const auto& error_by_substitution_type : errors_by_substitution_type) {
    // if there are no substitutions of this type, skip it
    if (error_by_substitution_type.substitutions_total == 0) {
      continue;
    }
    std::vector row = {error_by_substitution_type.GetSubstitutionString(),
                       std::to_string(error_by_substitution_type.substitutions_total)};
    if (dataset_metadata.has_cluster_info) {
      row.emplace_back(std::to_string(error_by_substitution_type.substitutions_mixed_strand_cluster));
      row.emplace_back(std::to_string(error_by_substitution_type.substitutions_forward_cluster));
      row.emplace_back(std::to_string(error_by_substitution_type.substitutions_reverse_cluster));
      row.emplace_back(std::to_string(error_by_substitution_type.substitutions_full_cluster));
      row.emplace_back(std::to_string(error_by_substitution_type.substitutions_mixed_full_and_partial_cluster));
      row.emplace_back(std::to_string(error_by_substitution_type.substitutions_partial_cluster));
    }
    if (dataset_metadata.has_read_type_info) {
      row.emplace_back(std::to_string(error_by_substitution_type.substitutions_full_read));
      row.emplace_back(std::to_string(error_by_substitution_type.substitutions_partial_read));
    }
    if (dataset_metadata.has_strand_info) {
      row.emplace_back(std::to_string(error_by_substitution_type.substitutions_forward_read));
      row.emplace_back(std::to_string(error_by_substitution_type.substitutions_reverse_read));
    }
    writer << row;
  }
}

ErrorsBySubstitutionType& ErrorsBySubstitutionType::operator+=(const ErrorsBySubstitutionType& other) {
  substitutions_total += other.substitutions_total;
  substitutions_mixed_strand_cluster += other.substitutions_mixed_strand_cluster;
  substitutions_forward_cluster += other.substitutions_forward_cluster;
  substitutions_reverse_cluster += other.substitutions_reverse_cluster;
  substitutions_full_cluster += other.substitutions_full_cluster;
  substitutions_mixed_full_and_partial_cluster += other.substitutions_mixed_full_and_partial_cluster;
  substitutions_partial_cluster += other.substitutions_partial_cluster;
  substitutions_forward_read += other.substitutions_forward_read;
  substitutions_reverse_read += other.substitutions_reverse_read;
  substitutions_full_read += other.substitutions_full_read;
  substitutions_partial_read += other.substitutions_partial_read;

  // sum all fields in the vector
  for (size_t i = 0; i < errors_by_substitution_type.size(); i++) {
    auto& error_by_substitution_type = errors_by_substitution_type[i];
    auto other_error_by_substitution_type = other.errors_by_substitution_type[i];
    error_by_substitution_type.substitutions_total += other_error_by_substitution_type.substitutions_total;
    error_by_substitution_type.substitutions_mixed_strand_cluster +=
        other_error_by_substitution_type.substitutions_mixed_strand_cluster;
    error_by_substitution_type.substitutions_forward_cluster +=
        other_error_by_substitution_type.substitutions_forward_cluster;
    error_by_substitution_type.substitutions_reverse_cluster +=
        other_error_by_substitution_type.substitutions_reverse_cluster;
    error_by_substitution_type.substitutions_full_cluster +=
        other_error_by_substitution_type.substitutions_full_cluster;
    error_by_substitution_type.substitutions_mixed_full_and_partial_cluster +=
        other_error_by_substitution_type.substitutions_mixed_full_and_partial_cluster;
    error_by_substitution_type.substitutions_partial_cluster +=
        other_error_by_substitution_type.substitutions_partial_cluster;
    error_by_substitution_type.substitutions_forward_read +=
        other_error_by_substitution_type.substitutions_forward_read;
    error_by_substitution_type.substitutions_reverse_read +=
        other_error_by_substitution_type.substitutions_reverse_read;
    error_by_substitution_type.substitutions_full_read += other_error_by_substitution_type.substitutions_full_read;
    error_by_substitution_type.substitutions_partial_read +=
        other_error_by_substitution_type.substitutions_partial_read;
  }

  return *this;
}

ErrorsBySubstitutionType& ErrorsBySubstitutionType::operator+=(
    const ankerl::unordered_dense::map<u16, CountsByReadProperty<u32>>& errors_by_substitution_type_map) {
  for (const auto& [substitution_index, base_error_stats] : errors_by_substitution_type_map) {
    if (substitution_index < errors_by_substitution_type.size()) {
      auto& error_by_substitution_type = errors_by_substitution_type[substitution_index];
      error_by_substitution_type.substitutions_total += base_error_stats.total;
      error_by_substitution_type.substitutions_mixed_strand_cluster += base_error_stats.mixed_strand_cluster;
      error_by_substitution_type.substitutions_forward_cluster += base_error_stats.forward_cluster;
      error_by_substitution_type.substitutions_reverse_cluster += base_error_stats.reverse_cluster;
      error_by_substitution_type.substitutions_full_cluster += base_error_stats.full_cluster;
      error_by_substitution_type.substitutions_mixed_full_and_partial_cluster +=
          base_error_stats.mixed_full_and_partial_cluster;
      error_by_substitution_type.substitutions_partial_cluster += base_error_stats.partial_cluster;
      error_by_substitution_type.substitutions_forward_read += base_error_stats.forward_read;
      error_by_substitution_type.substitutions_reverse_read += base_error_stats.reverse_read;
      error_by_substitution_type.substitutions_full_read += base_error_stats.full_read;
      error_by_substitution_type.substitutions_partial_read += base_error_stats.partial_read;
      // Update the global counts
      substitutions_total += base_error_stats.total;
      substitutions_mixed_strand_cluster += base_error_stats.mixed_strand_cluster;
      substitutions_forward_cluster += base_error_stats.forward_cluster;
      substitutions_reverse_cluster += base_error_stats.reverse_cluster;
      substitutions_full_cluster += base_error_stats.full_cluster;
      substitutions_mixed_full_and_partial_cluster += base_error_stats.mixed_full_and_partial_cluster;
      substitutions_partial_cluster += base_error_stats.partial_cluster;
      substitutions_forward_read += base_error_stats.forward_read;
      substitutions_reverse_read += base_error_stats.reverse_read;
      substitutions_full_read += base_error_stats.full_read;
      substitutions_partial_read += base_error_stats.partial_read;
    }
  }
  return *this;
}

// Returns the substitution string for the base for rendering with 'reference->substitution'
std::string SubstitutionError::GetSubstitutionString() const {
  return std::string{ref} + "->" + base;
}
}  // namespace xoos::alignment_metrics
