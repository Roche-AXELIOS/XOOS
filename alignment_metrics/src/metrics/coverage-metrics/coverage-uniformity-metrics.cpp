#include "metrics/coverage-metrics/coverage-uniformity-metrics.h"

#include <tuple>

#include <csv.hpp>

#include <xoos/histogram/histogram-summary.h>
#include <xoos/io/metadata-util.h>
#include <xoos/util/format-functions.h>
#include <xoos/util/missing-value.h>

#include "metrics/metrics-names.h"

namespace xoos::alignment_metrics {

using util::kNA;
using util::ToStringWithPrecision;

namespace {
// Percentile used to compute the median coverage.
constexpr u64 kMedianPercentile = 50;
}  // namespace

CoverageUniformityMetrics& CoverageUniformityMetrics::operator+=(const CoverageUniformityMetrics& other) {
  for (const auto& [mean_coverage_bin, count] : other.region_counts_by_mean_coverage) {
    region_counts_by_mean_coverage[mean_coverage_bin] += count;
  }
  return *this;
}

void CoverageUniformityMetrics::WriteMeanCoverageHistogram(const fs::path& output_path,
                                                           const io::TsvMetadata& tsv_metadata) const {
  auto out = std::ofstream(output_path);
  auto writer = csv::make_tsv_writer(out).set_auto_flush(false);
  // Write metadata
  io::WriteTsvMetadata(out, tsv_metadata);
  // Write header
  writer << std::vector{kNameMeanCoverage, kNameTargetRegionCount};
  // Write data
  for (const auto& [mean_coverage_bin, count] : region_counts_by_mean_coverage) {
    writer << std::make_tuple(mean_coverage_bin, count);
  }
}

void CoverageUniformityMetrics::WriteCoverageUniformitySummaryTsv(
    const fs::path& output_path,
    const io::TsvMetadata& tsv_metadata,
    const histogram::Histogram<u64>& filtered_coverage_histogram) const {
  const CoverageUniformitySummary summary(filtered_coverage_histogram, region_counts_by_mean_coverage);
  summary.WriteTsv(output_path, tsv_metadata);
}

namespace {
/**
 * Count the number of positions in the histogram with coverage in the specified range [lower_bound, upper_bound].
 */
u64 CountPositionsWithCoverageInRange(const histogram::Histogram<u64>& histogram, u64 lower_bound, u64 upper_bound) {
  u64 count = 0;
  for (size_t coverage = lower_bound;
       coverage < std::min(static_cast<size_t>(upper_bound + 1), histogram.counts.size());
       ++coverage) {
    count += histogram.At(coverage);
  }
  // count outliers
  for (const auto& [coverage, outlier_count] : histogram.outliers) {
    if (coverage >= lower_bound && coverage <= upper_bound) {
      count += outlier_count;
    }
  }
  return count;
}
}  // namespace

CoverageUniformitySummary::CoverageUniformitySummary(const histogram::Histogram<u64>& filtered_coverage_histogram,
                                                     const std::map<u64, u32>& mean_coverage_histogram)
    : mean_coverage(filtered_coverage_histogram.IsEmpty() ? 0.0 : histogram::ComputeMean(filtered_coverage_histogram)),
      median_coverage(filtered_coverage_histogram.IsEmpty()
                          ? 0
                          : histogram::ComputePercentile(filtered_coverage_histogram, kMedianPercentile)) {
  for (const auto& [mean_cov, count] : mean_coverage_histogram) {
    target_region_count += count;
    if (mean_cov == 0) {
      target_regions_with_no_coverage += count;
    }
  }

  // Compute fold_80_base_penalty and pct_bases_at_0_5x_to_2x_mean_coverage
  if (mean_coverage > 0 && !filtered_coverage_histogram.IsEmpty()) {
    const u64 coverage_at_20_percentile = histogram::ComputePercentile(filtered_coverage_histogram, 20);
    if (coverage_at_20_percentile > 0) {
      fold_80_base_penalty = mean_coverage / static_cast<f64>(coverage_at_20_percentile);
    }
    const u64 bases_in_range = CountPositionsWithCoverageInRange(
        filtered_coverage_histogram, static_cast<u64>(0.5 * mean_coverage), static_cast<u64>(2.0 * mean_coverage));
    const u64 total_bases = histogram::ComputeCount(filtered_coverage_histogram);
    pct_bases_at_0_5x_to_2x_mean_coverage = (static_cast<f64>(bases_in_range) / static_cast<f64>(total_bases)) * 100.0;
  }
}

void CoverageUniformitySummary::WriteTsv(const fs::path& output_path, const io::TsvMetadata& tsv_metadata) const {
  auto out = std::ofstream(output_path);
  auto writer = csv::make_tsv_writer(out).set_auto_flush(false);
  // Write metadata
  io::WriteTsvMetadata(out, tsv_metadata);
  // Write header
  writer << std::vector{kNameMetricName, kNameValue};
  // Write data
  writer << std::vector{kNameMeanCoverage, ToStringWithPrecision(mean_coverage, histogram::kMeanPrecision)};
  writer << std::make_tuple(kNameMedianCoverage, median_coverage);
  if (fold_80_base_penalty.has_value()) {
    writer << std::vector{kNameFold80BasePenalty, ToStringWithPrecision(fold_80_base_penalty.value())};
  } else {
    writer << std::vector{kNameFold80BasePenalty, std::string{kNA}};
  }
  writer << std::vector{kNamePctBasesAt05xTo2xMeanCoverage,
                        ToStringWithPrecision(pct_bases_at_0_5x_to_2x_mean_coverage)};
  writer << std::make_tuple(kNameTargetRegionCount, target_region_count);
  writer << std::make_tuple(kNameTargetRegionsWithNoCoverage, target_regions_with_no_coverage);
}

}  // namespace xoos::alignment_metrics
