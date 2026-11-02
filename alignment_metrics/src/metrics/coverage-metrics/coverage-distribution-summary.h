#pragma once

#include <map>
#include <string>

#include <xoos/histogram/histogram-summary.h>
#include <xoos/io/metadata-util.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "metadata/dataset-metadata.h"
#include "metrics/coverage-metrics/coverage-histograms.h"

namespace xoos::alignment_metrics {

/**
 * @brief Summary of coverage distribution statistics across genomic regions.
 *
 * Aggregates coverage histogram data into position counts for various coverage cutoff thresholds.
 * For each histogram and coverage cutoff, the struct calculates the number of positions with
 * depth of coverage greater than or equal to that cutoff. This enables analysis of coverage
 * uniformity and identification of under-covered regions.
 *
 * For every dataset type, tracks:
 * - Any coverage: all bases
 * - Post-filter coverage: high-confidence filtered bases
 */
struct CoverageDistributionSummary {
  // Constructor for manually creating a distribution summary for testing
  CoverageDistributionSummary() = default;

  /**
   * @brief Creates a distribution summary from coverage histograms.
   *
   * For each histogram and coverage cutoff, sums the number of positions with depth of coverage
   * greater than or equal to the cutoff. The dataset_metadata determines which histograms are populated.
   *
   * @param coverage_histograms The coverage histograms to summarize
   */
  explicit CoverageDistributionSummary(const CoverageHistograms& coverage_histograms);

  DatasetMetadata dataset_metadata;
  vec<u64> coverage_cutoffs;

  u64 total_positions{};
  u64 positions_with_any_coverage{};
  u64 positions_with_no_coverage{};
  u64 positions_with_post_filter_coverage{};
  u64 positions_with_no_post_filter_coverage{};
  std::map<u64, u64> any_coverage_cutoff_map;
  std::map<u64, u64> post_filter_coverage_cutoff_map;

  static std::vector<std::string> GetHeaders();
  void WriteTsv(const fs::path& output_path, const io::TsvMetadata& tsv_metadata) const;

  bool operator==(const CoverageDistributionSummary& other) const = default;
};

}  // namespace xoos::alignment_metrics
