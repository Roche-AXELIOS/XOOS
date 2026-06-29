#pragma once
#include <string>
#include <vector>

#include <xoos/io/metadata-util.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "core/error-counts.h"

namespace xoos::alignment_metrics {

/**
 * This struct stores a summary of base-level accuracy metrics.
 *
 * It includes the total number of aligned bases, the number of errors (substitutions, insertions, and deletions)
 * across all reads regardless of cluster size or read type.
 */
struct BaseLevelAccuracySummary {
  u64 total_bases{};
  u64 substitutions{};
  u64 insertions{};
  u64 deletions{};
  u64 inserted_bases{};
  u64 deleted_bases{};

  u64 total_positions{};
  u64 positions_skipped_due_to_low_depth{};
  u64 positions_skipped_due_to_high_alt_allele_fraction{};

  static std::vector<std::string> GetHeaders();
  void WriteTsv(const fs::path& output_path, const io::Comments& comments) const;

  BaseLevelAccuracySummary& operator+=(const BaseLevelAccuracySummary& other);
  BaseLevelAccuracySummary& operator+=(const CountsByErrorTypeAndReadProperty<u32>& combined_stats);
};

}  // namespace xoos::alignment_metrics
