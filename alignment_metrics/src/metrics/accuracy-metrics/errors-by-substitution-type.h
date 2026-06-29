#pragma once
#include <string>

#include <ankerl/unordered_dense.h>

#include <xoos/io/metadata-util.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "core/error-counts.h"
#include "metadata/dataset-metadata.h"

namespace xoos::alignment_metrics {

/**
 * A struct to hold counts of a specific substitution error type (e.g., A->C),
 * categorized by various read and cluster properties from which the error was observed.
 */
struct SubstitutionError {
  SubstitutionError() = default;
  char base{};
  char ref{};
  u64 substitutions_total{};
  u64 substitutions_mixed_strand_cluster{};
  u64 substitutions_forward_cluster{};
  u64 substitutions_reverse_cluster{};
  u64 substitutions_full_cluster{};
  u64 substitutions_mixed_full_and_partial_cluster{};
  u64 substitutions_partial_cluster{};
  u64 substitutions_forward_read{};
  u64 substitutions_reverse_read{};
  u64 substitutions_full_read{};
  u64 substitutions_partial_read{};

  // Gets the substitution string for the base for rendering
  std::string GetSubstitutionString() const;
  bool operator==(const SubstitutionError& other) const = default;
};

/**
 * Container for counts of substitutions of each type (e.g., A->C, G->T)
 * further categorized by various read and cluster properties.
 */
struct ErrorsBySubstitutionType {
  ErrorsBySubstitutionType() = default;
  explicit ErrorsBySubstitutionType(const DatasetMetadata& dataset_metadata);

  // Dataset metadata used to determine which columns to omit
  DatasetMetadata dataset_metadata{};

  u64 substitutions_total{};
  u64 substitutions_mixed_strand_cluster{};
  u64 substitutions_forward_cluster{};
  u64 substitutions_reverse_cluster{};
  u64 substitutions_full_cluster{};
  u64 substitutions_mixed_full_and_partial_cluster{};
  u64 substitutions_partial_cluster{};
  u64 substitutions_forward_read{};
  u64 substitutions_reverse_read{};
  u64 substitutions_full_read{};
  u64 substitutions_partial_read{};

  /**
   * A vector of substitution errors, where each entry corresponds to a specific substitution type (e.g., A->C, G->T).
   * Each entry contains counts of substitutions categorized by various read and cluster properties.
   */
  std::vector<SubstitutionError> errors_by_substitution_type{};

  std::vector<std::string> GetHeaders() const;
  void WriteTsv(const fs::path& output_path, const io::Comments& comments) const;

  ErrorsBySubstitutionType& operator+=(const ErrorsBySubstitutionType& other);
  ErrorsBySubstitutionType& operator+=(
      const ankerl::unordered_dense::map<u16, CountsByReadProperty<u32>>& errors_by_substitution_type_map);
};
}  // namespace xoos::alignment_metrics
