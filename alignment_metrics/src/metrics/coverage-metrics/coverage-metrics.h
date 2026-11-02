#pragma once

#include <optional>

#include <xoos/types/int.h>

#include "alignment-metrics-options.h"
#include "metadata/dataset-metadata.h"
#include "metrics/coverage-metrics/coverage-histograms.h"
#include "metrics/coverage-metrics/coverage-uniformity-metrics.h"

namespace xoos::alignment_metrics {

/**
 * @brief Container for all coverage-related metrics and histograms.
 *
 * Tracks depth of coverage distributions across genomic regions. Includes separate histograms for
 * homopolymer regions when HP analysis is enabled, and coverage uniformity metrics for target enrichment experiments.
 */
class CoverageMetrics {
 public:
  CoverageMetrics() = default;
  explicit CoverageMetrics(const AlignmentMetricsOptions& options, const DatasetMetadata& dataset_metadata);

  CoverageHistograms coverage_histograms;
  std::optional<CoverageHistograms> hp_coverage_histograms;
  std::optional<CoverageUniformityMetrics> coverage_uniformity_metrics;

  /**
   * Increment the histogram bin for 0-depth positions by the given number of bases.
   */
  void AddEmptyHistogramData(u64 num_bases);
  /**
   * Update the coverage histograms with the given depth values
   */
  void AddHistogramData(u32 all_depth, u32 concordant_duplex_depth, u32 filtered_depth, bool is_part_of_hp = false);
  CoverageMetrics& operator+=(const CoverageMetrics& other);
};

}  // namespace xoos::alignment_metrics
