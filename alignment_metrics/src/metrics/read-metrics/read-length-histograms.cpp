#include "metrics/read-metrics/read-length-histograms.h"

#include <filesystem>
#include <tuple>
#include <vector>

#include <csv.hpp>

#include <xoos/types/int.h>

#include "metrics/metrics-names.h"

namespace xoos::alignment_metrics {

void ReadLengthHistograms::WriteTsv(const fs::path& output_path, const io::Comments& comments) const {
  histogram::Histograms histograms{std::make_tuple(kNamePostFilterAllReadCount, post_filter_all_read_length_histogram)};
  if (dataset_metadata.has_read_type_info) {
    histograms.emplace_back(kNamePostFilterPartialReadCount, post_filter_partial_read_length_histogram);
    histograms.emplace_back(kNamePostFilterFullReadCount, post_filter_full_read_length_histogram);
  }
  histogram::WriteHistograms(histograms,
                             output_path,
                             kNameReadLengthExcludingSoftClips,
                             histogram::HistogramBinOutput::kOmitFirstBinAndMaxLastBin,
                             {},
                             comments);
}

void ReadLengthHistograms::WriteSummaryTsv(const fs::path& output_path, const io::Comments& comments) const {
  // Make sure we always calculate the 10th, 50th, and 90th percentile
  // because we need these values for the ratio calculations.
  // NOTE: These percentiles are not necessarily written to the output file
  // unless they are explicitly requested by the user.
  std::vector<u64> required_percentiles = {10, 50, 90};
  // Add any additional percentiles specified in the summary_stats_percentiles
  for (const auto& percentile : summary_stats_percentiles) {
    if (std::ranges::find(required_percentiles, percentile) == required_percentiles.end()) {
      required_percentiles.emplace_back(percentile);
    }
  }
  // Calculate summary statistics for the full coverage histogram
  // Include required percentiles in addition to the summary_stats_percentiles
  histogram::HistogramSummaries summary_histograms{std::make_tuple(
      kNamePostFilterAllReadsLengthExcludingSoftClips,
      histogram::CalculateHistogramSummary(post_filter_all_read_length_histogram, {}, required_percentiles))};
  if (dataset_metadata.has_read_type_info) {
    summary_histograms.emplace_back(
        kNamePostFilterPartialReadsLengthExcludingSoftClips,
        histogram::CalculateHistogramSummary(post_filter_partial_read_length_histogram, {}, required_percentiles));
    summary_histograms.emplace_back(
        kNamePostFilterFullReadsLengthExcludingSoftClips,
        histogram::CalculateHistogramSummary(post_filter_full_read_length_histogram, {}, required_percentiles));
  }
  // When writing to the output, we can omit the required percentiles if they are not
  // specifically requested by the user. The values are still implicitly stored in
  // the HistogramSummary struct for computing percentile ratios.
  histogram::WriteSummaryHistograms(summary_histograms, output_path, summary_stats_percentiles, {}, false, comments);
}

}  // namespace xoos::alignment_metrics
