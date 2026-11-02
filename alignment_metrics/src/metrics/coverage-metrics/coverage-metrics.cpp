#include "metrics/coverage-metrics/coverage-metrics.h"

#include <optional>

#include <xoos/types/int.h>

#include "alignment-metrics-options.h"
#include "metadata/dataset-metadata.h"

namespace xoos::alignment_metrics {

CoverageMetrics::CoverageMetrics(const AlignmentMetricsOptions& options, const DatasetMetadata& dataset_metadata)
    : coverage_histograms(CoverageHistograms(options.max_coverage_bin,
                                             options.coverage_cutoff,
                                             options.summary_stats_percentiles,
                                             dataset_metadata,
                                             options.bed_input.empty())),
      hp_coverage_histograms(options.calculate_hp_metrics
                                 ? std::make_optional(CoverageHistograms(options.max_coverage_bin,
                                                                         options.coverage_cutoff,
                                                                         options.summary_stats_percentiles,
                                                                         dataset_metadata,
                                                                         options.bed_input.empty()))
                                 : std::nullopt),
      coverage_uniformity_metrics(options.enable_te_metrics ? std::make_optional(CoverageUniformityMetrics())
                                                            : std::nullopt) {
}

CoverageMetrics& CoverageMetrics::operator+=(const CoverageMetrics& other) {
  coverage_histograms.full_coverage_histogram += other.coverage_histograms.full_coverage_histogram;
  coverage_histograms.filtered_coverage_histogram += other.coverage_histograms.filtered_coverage_histogram;
  if (coverage_histograms.dataset_metadata.is_duplex_dataset &&
      other.coverage_histograms.dataset_metadata.is_duplex_dataset) {
    coverage_histograms.concordant_duplex_coverage_histogram +=
        other.coverage_histograms.concordant_duplex_coverage_histogram;
  }
  if (hp_coverage_histograms && other.hp_coverage_histograms) {
    hp_coverage_histograms->full_coverage_histogram += other.hp_coverage_histograms->full_coverage_histogram;
    hp_coverage_histograms->filtered_coverage_histogram += other.hp_coverage_histograms->filtered_coverage_histogram;
    if (hp_coverage_histograms->dataset_metadata.is_duplex_dataset &&
        other.hp_coverage_histograms->dataset_metadata.is_duplex_dataset) {
      hp_coverage_histograms->concordant_duplex_coverage_histogram +=
          other.hp_coverage_histograms->concordant_duplex_coverage_histogram;
    }
  }
  if (coverage_uniformity_metrics && other.coverage_uniformity_metrics) {
    coverage_uniformity_metrics->operator+=(*other.coverage_uniformity_metrics);
  }
  return *this;
}

void CoverageMetrics::AddEmptyHistogramData(const u64 num_bases) {
  coverage_histograms.full_coverage_histogram.AddCountToHistogram(0, num_bases);
  coverage_histograms.filtered_coverage_histogram.AddCountToHistogram(0, num_bases);
  if (coverage_histograms.dataset_metadata.is_duplex_dataset) {
    coverage_histograms.concordant_duplex_coverage_histogram.AddCountToHistogram(0, num_bases);
  }
}

/**
 * @brief Adds coverage depth counts to coverage histograms for a single genomic position.
 *
 * @param all_depth Total depth including all bases (low-quality and discordant)
 * @param concordant_duplex_depth Depth of concordant duplex bases only
 * @param filtered_depth Depth of high-quality filtered bases
 * @param is_part_of_hp Whether this position is part of a homopolymer sequence
 */
void CoverageMetrics::AddHistogramData(const u32 all_depth,
                                       const u32 concordant_duplex_depth,
                                       const u32 filtered_depth,
                                       const bool is_part_of_hp) {
  if (coverage_histograms.dataset_metadata.is_duplex_dataset) {
    coverage_histograms.concordant_duplex_coverage_histogram.AddCountToHistogram(concordant_duplex_depth, u64{1});
    if (is_part_of_hp && hp_coverage_histograms) {
      hp_coverage_histograms->concordant_duplex_coverage_histogram.AddCountToHistogram(concordant_duplex_depth, u64{1});
    }
  }
  coverage_histograms.filtered_coverage_histogram.AddCountToHistogram(filtered_depth, u64{1});
  coverage_histograms.full_coverage_histogram.AddCountToHistogram(all_depth, u64{1});
  if (is_part_of_hp && hp_coverage_histograms) {
    hp_coverage_histograms->filtered_coverage_histogram.AddCountToHistogram(filtered_depth, u64{1});
    hp_coverage_histograms->full_coverage_histogram.AddCountToHistogram(all_depth, u64{1});
  }
}

}  // namespace xoos::alignment_metrics
