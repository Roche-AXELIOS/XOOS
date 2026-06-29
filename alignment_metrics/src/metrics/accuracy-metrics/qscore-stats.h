#pragma once
#include <string>

#include <ankerl/unordered_dense.h>

#include <xoos/io/metadata-util.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "core/error-counts.h"
#include "metadata/dataset-metadata.h"

namespace xoos::alignment_metrics {

// Maximum possible qscore value in the Phred+33 encoding
constexpr u8 kMaxQscoreValue = 93;

/**
 * @brief Statistics tracking sequencing accuracy by base quality score.
 *
 * Aggregates match and mismatch counts for each possible quality score value (0-93 in Phred+33 encoding).
 */
struct QscoreStats {
  // for each qscore, the count of matches and mismatches with that qscore
  vec<MismatchCounts<u64>> mismatches_by_qscore;

  DatasetMetadata dataset_metadata{};

  QscoreStats() = default;
  explicit QscoreStats(const DatasetMetadata& dataset_metadata);

  vec<std::string> GetHeaders() const;
  void WriteTsv(const fs::path& output_path, const io::Comments& comments) const;

  QscoreStats& operator+=(const vec<MismatchCounts<u64>>& other);
  QscoreStats& operator+=(const ankerl::unordered_dense::map<u8, MismatchCounts<u32>>& other);
  QscoreStats& operator+=(const QscoreStats& other);
};

}  // namespace xoos::alignment_metrics
