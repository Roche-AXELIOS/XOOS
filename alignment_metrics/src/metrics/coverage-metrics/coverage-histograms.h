#pragma once

#include <vector>

#include <xoos/histogram/histogram-summary.h>
#include <xoos/io/metadata-util.h>
#include <xoos/types/int.h>

#include "alignment-metrics-options.h"
#include "metadata/dataset-metadata.h"

namespace xoos::alignment_metrics {

/**
 * A struct to hold coverage histograms for a dataset.
 *
 * It contains histograms for all bases and for bases that pass the configured filters.
 *
 * It contains a summary of coverage statistics with basic information like min, max, mean, median,
 * and percentiles. It also contains a summary of the coverage distribution, which includes information on
 * the number and percentage of positions have coverage at least a certain cutoff.
 */
class CoverageHistograms {
 public:
  CoverageHistograms() = default;

  /**
   * Construct and initialize a CoverageHistograms object.
   *
   * @param max_depth Maximum depth bin (i.e., number of bins - 1) for the coverage histograms
   * @param coverage_cutoffs Coverage cutoffs to be included in the coverage distribution summary
   * @param summary_stats_percentiles Percentiles to be included in the coverage statistics summary
   * @param dataset_metadata Metadata for the dataset.
   * @param exclude_zero_bin Whether to exclude the zero-coverage bin when calculating summary statistics. We include
   * positions with zero coverage when calculating mean and median coverage only if this is set to false. Should be
   * false when calculating coverage statistics for targeted (TE panel or exome) sequencing datasets and true otherwise.
   */
  explicit CoverageHistograms(u32 max_depth,
                              const std::vector<u64>& coverage_cutoffs = kDefaultCoverageCutoffs,
                              const std::vector<u64>& summary_stats_percentiles = kDefaultSummaryStatsPercentiles,
                              const DatasetMetadata& dataset_metadata = DatasetMetadata(),
                              bool exclude_zero_bin = true);

  // CoverageHistograms is a data holder whose histograms are accumulated in place (including via
  // operator+=) and then written out; these fields are intentionally public.
  // NOLINTBEGIN(misc-non-private-member-variables-in-classes,cppcoreguidelines-non-private-member-variables-in-classes)
  // Dataset metadata associated with the input BAM.
  DatasetMetadata dataset_metadata{};
  // Coverage cutoffs to be included in the summary
  std::vector<u64> coverage_cutoffs;
  // Percentiles to be included in the summary
  std::vector<u64> summary_stats_percentiles;
  // Coverage histograms for all bases (including low-quality discordant duplex bases
  // and simplex bases in the case of duplex; or all bases in the case of simplex)
  histogram::Histogram<u64> full_coverage_histogram;
  // Coverage histograms for non-YC filtered bases (only applicable to bam without yc tags)
  histogram::Histogram<u64> filtered_coverage_histogram;
  // Coverage histograms for concordant duplex bases, retained for internal pileup accounting.
  histogram::Histogram<u64> concordant_duplex_coverage_histogram;
  // NOLINTEND(misc-non-private-member-variables-in-classes,cppcoreguidelines-non-private-member-variables-in-classes)

  /**
   * Output the coverage histograms to a TSV file.
   *
   * The output contains columns for any coverage and post-filter coverage.
   */
  void WriteTsv(const fs::path& output_path, const io::TsvMetadata& tsv_metadata) const;
  /**
   * Output a summary of the coverage statistics to a TSV file.
   *
   * The summary statistics include: min, max, mean, median, and
   * a list of percentiles as specified by the user when creating the object.
   * Additionally, we output the ratio between the 90th and 10th percentiles,
   * and the ratio between the median and the 10th percentile.
   *
   * The output contains statistics for any coverage and post-filter coverage.
   */
  void WriteCoverageStatsTsv(const fs::path& output_path, const io::TsvMetadata& tsv_metadata) const;
  /**
   * Output a summary of the coverage distribution to a TSV file.
   *
   * The summary distribution contains information on how many positions have
   * coverage at least a certain cutoff.
   *
   * The output contains summary distributions for any coverage and post-filter coverage.
   */
  void WriteCoverageDistributionSummaryTsv(const fs::path& output_path, const io::TsvMetadata& tsv_metadata) const;

  CoverageHistograms& operator+=(const CoverageHistograms& other);

 private:
  bool _exclude_zero_bin = false;
};

inline u64 CountPositionsWithNoCoverage(const histogram::Histogram<u64>& histogram) {
  return histogram.At(0);
}

inline u64 CountPositionsWithCoverageAtLeast(const histogram::Histogram<u64>& histogram, u64 coverage_cutoff) {
  u64 count = 0;
  for (size_t coverage = coverage_cutoff; coverage < histogram.counts.size(); ++coverage) {
    count += histogram.At(coverage);
  }
  // count outliers
  for (const auto& [coverage, outlier_count] : histogram.outliers) {
    if (coverage >= coverage_cutoff) {
      count += outlier_count;
    }
  }
  return count;
}

}  // namespace xoos::alignment_metrics
