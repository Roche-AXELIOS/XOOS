#include "metrics/coverage-metrics/coverage-distribution-summary.h"

#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "metrics/coverage-metrics/coverage-histograms.h"
#include "metrics/metrics-names.h"
#include "util/format-util.h"

namespace xoos::alignment_metrics {

constexpr u8 kPercentagePrecision = 2;

CoverageDistributionSummary::CoverageDistributionSummary(const CoverageHistograms& coverage_histograms)
    : dataset_metadata(coverage_histograms.dataset_metadata),
      coverage_cutoffs(coverage_histograms.coverage_cutoffs),
      total_positions(ComputeCount(coverage_histograms.full_coverage_histogram)) {
  // Also update the total positions, positions with any coverage, and positions with no coverage
  // for the full coverage histogram
  positions_with_any_coverage = CountPositionsWithCoverageAtLeast(coverage_histograms.full_coverage_histogram, u64{1});
  positions_with_no_coverage = CountPositionsWithNoCoverage(coverage_histograms.full_coverage_histogram);
  positions_with_post_filter_coverage =
      CountPositionsWithCoverageAtLeast(coverage_histograms.filtered_coverage_histogram, u64{1});
  positions_with_no_post_filter_coverage =
      CountPositionsWithNoCoverage(coverage_histograms.filtered_coverage_histogram);
  // Count number of positions with and without coverage for the concordant duplex coverage histogram
  // if the dataset is duplex
  if (coverage_histograms.dataset_metadata.is_duplex_dataset) {
    positions_with_concordant_duplex_coverage =
        CountPositionsWithCoverageAtLeast(coverage_histograms.concordant_duplex_coverage_histogram, u64{1});
    positions_with_no_concordant_duplex_coverage =
        CountPositionsWithNoCoverage(coverage_histograms.concordant_duplex_coverage_histogram);
    concordant_duplex_coverage_cutoff_map = std::map<u64, u64>{};
  }
  // For each cutoff value, count the number of positions with coverage at least that cutoff
  for (const auto& cutoff : coverage_histograms.coverage_cutoffs) {
    any_coverage_cutoff_map[cutoff] =
        CountPositionsWithCoverageAtLeast(coverage_histograms.full_coverage_histogram, cutoff);
    post_filter_coverage_cutoff_map[cutoff] =
        CountPositionsWithCoverageAtLeast(coverage_histograms.filtered_coverage_histogram, cutoff);
    if (coverage_histograms.dataset_metadata.is_duplex_dataset) {
      concordant_duplex_coverage_cutoff_map.value()[cutoff] =
          CountPositionsWithCoverageAtLeast(coverage_histograms.concordant_duplex_coverage_histogram, cutoff);
    }
  }
}

vec<std::string> CoverageDistributionSummary::GetHeaders() const {
  vec headers = {kNameMetricName, kNameAnyCoverage, kNamePostFilterCoverage};
  if (dataset_metadata.is_duplex_dataset) {
    headers.emplace_back(kNameConcordantDuplexCoverage);
  }
  return headers;
}

void CoverageDistributionSummary::WriteTsv(const fs::path& output_path, const io::Comments& comments) const {
  std::ofstream out(output_path);
  auto writer = csv::make_tsv_writer_buffered(out);
  // Write comments
  io::WriteTsvComments(writer, comments);
  // Write header
  writer << GetHeaders();
  // Write total positions, positions with coverage, and positions with no coverage
  vec total_positions_row{kNameTotalPositions, std::to_string(total_positions), std::to_string(total_positions)};
  vec positions_with_no_coverage_row{kNamePositionsWithNoCoverage,
                                     std::to_string(positions_with_no_coverage),
                                     std::to_string(positions_with_no_post_filter_coverage)};
  vec percentage_of_positions_with_no_coverage_row{
      kNamePercentageOfPositionsWithNoCoverage,
      ToPercentageWithPrecision(
          static_cast<f64>(positions_with_no_coverage), static_cast<f64>(total_positions), kPercentagePrecision),
      ToPercentageWithPrecision(static_cast<f64>(positions_with_no_post_filter_coverage),
                                static_cast<f64>(total_positions),
                                kPercentagePrecision)};
  vec positions_with_coverage_row{kNamePositionsWithCoverage,
                                  std::to_string(positions_with_any_coverage),
                                  std::to_string(positions_with_post_filter_coverage)};
  vec percentage_of_positions_with_coverage_row{
      kNamePercentageOfPositionsWithCoverage,
      ToPercentageWithPrecision(
          static_cast<f64>(positions_with_any_coverage), static_cast<f64>(total_positions), kPercentagePrecision),
      ToPercentageWithPrecision(static_cast<f64>(positions_with_post_filter_coverage),
                                static_cast<f64>(total_positions),
                                kPercentagePrecision)};
  if (dataset_metadata.is_duplex_dataset) {
    total_positions_row.emplace_back(std::to_string(total_positions));
    positions_with_no_coverage_row.emplace_back(std::to_string(positions_with_no_concordant_duplex_coverage.value()));
    percentage_of_positions_with_no_coverage_row.emplace_back(
        ToPercentageWithPrecision(static_cast<f64>(positions_with_no_concordant_duplex_coverage.value()),
                                  static_cast<f64>(total_positions),
                                  kPercentagePrecision));
    positions_with_coverage_row.emplace_back(std::to_string(positions_with_concordant_duplex_coverage.value()));
    percentage_of_positions_with_coverage_row.emplace_back(
        ToPercentageWithPrecision(static_cast<f64>(positions_with_concordant_duplex_coverage.value()),
                                  static_cast<f64>(total_positions),
                                  kPercentagePrecision));
  }
  writer << total_positions_row;
  writer << positions_with_no_coverage_row;
  writer << percentage_of_positions_with_no_coverage_row;
  writer << positions_with_coverage_row;
  writer << percentage_of_positions_with_coverage_row;

  // Write coverage cutoffs
  for (const u64 cutoff : coverage_cutoffs) {
    vec positions_with_at_least_x_coverage_row{fmt::format(kNamePositionsAtLeastXCoverage, cutoff),
                                               std::to_string(any_coverage_cutoff_map.at(cutoff)),
                                               std::to_string(post_filter_coverage_cutoff_map.at(cutoff))};
    vec percentage_of_positions_at_least_x_coverage_row{
        fmt::format(kNamePercentageOfPositionsAtLeastXCoverage, cutoff),
        ToPercentageWithPrecision(static_cast<f64>(any_coverage_cutoff_map.at(cutoff)),
                                  static_cast<f64>(total_positions),
                                  kPercentagePrecision),
        ToPercentageWithPrecision(static_cast<f64>(post_filter_coverage_cutoff_map.at(cutoff)),
                                  static_cast<f64>(total_positions),
                                  kPercentagePrecision)};
    if (concordant_duplex_coverage_cutoff_map.has_value()) {
      positions_with_at_least_x_coverage_row.emplace_back(
          std::to_string(concordant_duplex_coverage_cutoff_map.value().at(cutoff)));
      percentage_of_positions_at_least_x_coverage_row.emplace_back(
          ToPercentageWithPrecision(static_cast<f64>(concordant_duplex_coverage_cutoff_map.value().at(cutoff)),
                                    static_cast<f64>(total_positions),
                                    kPercentagePrecision));
    }
    writer << positions_with_at_least_x_coverage_row;
    writer << percentage_of_positions_at_least_x_coverage_row;
  }
}

}  // namespace xoos::alignment_metrics
