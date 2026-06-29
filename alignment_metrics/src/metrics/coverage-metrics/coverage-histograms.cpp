#include "metrics/coverage-metrics/coverage-histograms.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>

#include <csv.hpp>

#include <xoos/types/int.h>

#include "coverage-distribution-summary.h"
#include "metrics/metrics-names.h"
#include "util/format-util.h"

namespace xoos::alignment_metrics {

CoverageHistograms::CoverageHistograms(u32 max_depth,
                                       const std::vector<u64>& coverage_cutoffs,
                                       const std::vector<u64>& summary_stats_percentiles,
                                       const DatasetMetadata& dataset_metadata,
                                       bool exclude_zero_bin)

    : dataset_metadata(dataset_metadata),
      coverage_cutoffs{coverage_cutoffs},
      summary_stats_percentiles{summary_stats_percentiles},
      full_coverage_histogram{max_depth + 1, 0},
      filtered_coverage_histogram{max_depth + 1, 0},
      concordant_duplex_coverage_histogram{dataset_metadata.is_duplex_dataset ? max_depth + 1 : 0, 0},
      _exclude_zero_bin{exclude_zero_bin} {
}

void CoverageHistograms::WriteTsv(const fs::path& output_path, const io::Comments& comments) const {
  histogram::Histograms<u64> histograms{std::make_tuple(kNamePositionCount, full_coverage_histogram),
                                        std::make_tuple(kNamePostFilterPositionCount, filtered_coverage_histogram)};
  if (dataset_metadata.is_duplex_dataset) {
    histograms.emplace_back(kNameConcordantDuplexCount, concordant_duplex_coverage_histogram);
  }

  histogram::WriteHistograms(
      histograms, output_path, kNameCoverage, histogram::HistogramBinOutput::kMaxLastBin, {}, comments);
}

void CoverageHistograms::WriteCoverageStatsTsv(const fs::path& output_path, const io::Comments& comments) const {
  // Make sure we always calculate the 10th, 50th, and 90th percentile
  // because we need these values for the ratio calculations.
  // NOTE: These percentiles are not necessarily written to the output file
  // unless they are explicitly requested by the user.
  std::vector<u64> required_percentiles = {10, 50, 90};
  required_percentiles.reserve(required_percentiles.size() + summary_stats_percentiles.size());

  // Add any additional percentiles specified in the summary_stats_percentiles
  for (const auto percentile : summary_stats_percentiles) {
    if (std::ranges::find(required_percentiles, percentile) == required_percentiles.end()) {
      required_percentiles.emplace_back(percentile);
    }
  }

  // Calculate summary statistics for any coverage histogram and filtered coverage histogram
  // we only include bins with non-zero counts for any coverage histogram
  histogram::HistogramSummaries summary_histograms = {
      std::make_tuple(
          kNameAnyCoverage,
          histogram::CalculateHistogramSummary(full_coverage_histogram, {}, required_percentiles, _exclude_zero_bin)),
      std::make_tuple(kNamePostFilterCoverage,
                      histogram::CalculateHistogramSummary(
                          filtered_coverage_histogram, {}, required_percentiles, _exclude_zero_bin))};
  // Calculate summary statistics for the concordant duplex coverage histogram
  if (dataset_metadata.is_duplex_dataset) {
    // we only include bins with non-zero counts for concordant duplex coverage histogram
    summary_histograms.emplace_back(
        kNameConcordantDuplexCoverage,
        histogram::CalculateHistogramSummary(
            concordant_duplex_coverage_histogram, {}, required_percentiles, _exclude_zero_bin));
  }
  // When writing to the output, we can omit the required percentiles if they are not
  // specifically requested by the user. The values are still implicitly stored in
  // the HistogramSummary struct for computing percentile ratios.
  histogram::WriteSummaryHistograms(summary_histograms, output_path, summary_stats_percentiles, {}, false, comments);

  // Write percentile ratios: ratio_90_to_10percentile and ratio_median_to_10percentile
  std::ofstream out(output_path, std::ios_base::app);
  auto writer = csv::make_tsv_writer_buffered(out);

  auto ratio_90_to_10percentile_full_histogram =
      histogram::GetRatioPercentile(std::get<1>(summary_histograms.at(0)).percentiles, 90, 10);
  auto ratio_median_to_10percentile_full_histogram =
      histogram::GetRatioPercentile(std::get<1>(summary_histograms.at(0)).percentiles, 50, 10);
  auto ratio_90_to_10percentile_filtered_histogram =
      histogram::GetRatioPercentile(std::get<1>(summary_histograms.at(1)).percentiles, 90, 10);
  auto ratio_median_to_10percentile_filtered_histogram =
      histogram::GetRatioPercentile(std::get<1>(summary_histograms.at(1)).percentiles, 50, 10);

  std::vector ratio_90_to_10_row = {fmt::format(kNameRatioXtoYPercentile, 90, 10),
                                    ratio_90_to_10percentile_full_histogram.has_value()
                                        ? std::to_string(ratio_90_to_10percentile_full_histogram.value())
                                        : kNotApplicable,
                                    ratio_90_to_10percentile_filtered_histogram.has_value()
                                        ? std::to_string(ratio_90_to_10percentile_filtered_histogram.value())
                                        : kNotApplicable};
  std::vector ratio_median_to_10_row = {fmt::format(kNameRatioXtoYPercentile, "median", 10),
                                        ratio_median_to_10percentile_full_histogram.has_value()
                                            ? std::to_string(ratio_median_to_10percentile_full_histogram.value())
                                            : kNotApplicable,
                                        ratio_median_to_10percentile_filtered_histogram.has_value()
                                            ? std::to_string(ratio_median_to_10percentile_filtered_histogram.value())
                                            : kNotApplicable};

  if (dataset_metadata.is_duplex_dataset) {
    auto ratio_90_to_10percentile_concordant_duplex_histogram =
        histogram::GetRatioPercentile(std::get<1>(summary_histograms.at(2)).percentiles, 90, 10);
    auto ratio_median_to_10percentile_concordant_duplex_histogram =
        histogram::GetRatioPercentile(std::get<1>(summary_histograms.at(2)).percentiles, 50, 10);
    ratio_90_to_10_row.emplace_back(ratio_90_to_10percentile_concordant_duplex_histogram.has_value()
                                        ? std::to_string(ratio_90_to_10percentile_concordant_duplex_histogram.value())
                                        : kNotApplicable);
    ratio_median_to_10_row.emplace_back(
        ratio_median_to_10percentile_concordant_duplex_histogram.has_value()
            ? std::to_string(ratio_median_to_10percentile_concordant_duplex_histogram.value())
            : kNotApplicable);
  }
  writer << ratio_90_to_10_row;
  writer << ratio_median_to_10_row;
}

void CoverageHistograms::WriteCoverageDistributionSummaryTsv(const fs::path& output_path,
                                                             const io::Comments& comments) const {
  const CoverageDistributionSummary summary(*this);
  summary.WriteTsv(output_path, comments);
}

CoverageHistograms& CoverageHistograms::operator+=(const CoverageHistograms& other) {
  full_coverage_histogram += other.full_coverage_histogram;
  concordant_duplex_coverage_histogram += other.concordant_duplex_coverage_histogram;
  filtered_coverage_histogram += other.filtered_coverage_histogram;
  return *this;
}

}  // namespace xoos::alignment_metrics
