#pragma once
#include <string>
#include <vector>

#include <xoos/io/metadata-util.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "core/error-counts.h"
#include "metadata/dataset-metadata.h"

namespace xoos::alignment_metrics {

/**
 * Container for counts of errors of each type
 * categorized by the read type from which the error was observed
 */
struct ErrorsByReadType {
  explicit ErrorsByReadType(const DatasetMetadata& dataset_metadata = DatasetMetadata())
      : dataset_metadata(dataset_metadata) {
  }

  // Dataset metadata used to determine which rows to omit
  DatasetMetadata dataset_metadata;

  u64 substitutions_total{};
  u64 substitutions_forward_read{};
  u64 substitutions_reverse_read{};
  u64 substitutions_full_read{};
  u64 substitutions_partial_read{};
  u64 insertions_total{};
  u64 insertions_forward_read{};
  u64 insertions_reverse_read{};
  u64 insertions_full_read{};
  u64 insertions_partial_read{};
  u64 deletions_total{};
  u64 deletions_forward_read{};
  u64 deletions_reverse_read{};
  u64 deletions_full_read{};
  u64 deletions_partial_read{};

  static std::vector<std::string> GetHeaders();
  void WriteTsv(const fs::path& output_path, const io::Comments& comments) const;
  ErrorsByReadType& operator+=(const ErrorsByReadType& other);
  ErrorsByReadType& operator+=(const CountsByErrorTypeAndReadProperty<u32>& combined_stats);
};
}  // namespace xoos::alignment_metrics
