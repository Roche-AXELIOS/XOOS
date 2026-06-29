#pragma once

#include <xoos/histogram/histogram-summary.h>
#include <xoos/io/metadata-util.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

namespace xoos::alignment_metrics {

/**
 * A struct to hold coverage uniformity metrics aggregated across target regions.
 * It includes a histogram of target region counts by mean coverage and summary statistics
 * like mean coverage, median coverage, fold 80 base penalty, etc.
 */
struct CoverageUniformityMetrics {
  std::map<u64, u32> region_counts_by_mean_coverage{};

  void WriteMeanCoverageHistogram(const fs::path& output_path, const io::Comments& comments) const;

  void WriteCoverageUniformitySummaryTsv(const fs::path& output_path,
                                         const io::Comments& comments,
                                         const histogram::Histogram<u64>& filtered_coverage_histogram) const;

  CoverageUniformityMetrics& operator+=(const CoverageUniformityMetrics& other);
};

/**
 * A struct to hold coverage uniformity summary metrics with metrics like mean coverage, median coverage, fold 80 base
 * penalty, etc.
 */
struct CoverageUniformitySummary {
  u64 mean_coverage{0};
  u64 median_coverage{0};
  u32 target_region_count{0};
  u32 target_regions_with_no_coverage{0};
  // nullopt if mean_coverage is 0 or if the 20th percentile coverage is 0
  std::optional<f64> fold_80_base_penalty{};
  f64 pct_bases_at_0_5x_to_2x_mean_coverage{0.0};

  /**
   * Construct a CoverageUniformitySummary from the filtered coverage histogram and mean coverage histogram (mean
   * coverage by target region).
   */
  CoverageUniformitySummary(const histogram::Histogram<u64>& filtered_coverage_histogram,
                            const std::map<u64, u32>& mean_coverage_histogram);

  void WriteTsv(const fs::path& output_path, const io::Comments& comments) const;
};

}  // namespace xoos::alignment_metrics
