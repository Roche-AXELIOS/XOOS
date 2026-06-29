#pragma once

#include <vector>

#include <xoos/histogram/histogram-summary.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "alignment-metrics-options.h"
#include "metadata/dataset-metadata.h"

namespace xoos::alignment_metrics {

/**
 * A struct to hold read length histograms for a dataset.
 *
 * If the dataset contains read type information (i.e. whether a read is full length or partial),
 * it will contain the following histograms:
 * - partial_read_length_histogram: Histogram for partial reads.
 * - full_read_length_histogram: Histogram for full reads.
 * - all_read_length_histogram: Combined histogram for all reads.
 *
 * Otherwise, it will only contain all_read_length_histogram
 *
 * It contains a summary of read length statistics with basic information like min, max, mean, median,
 * and percentiles.
 */
struct ReadLengthHistograms {
  ReadLengthHistograms() = default;

  explicit ReadLengthHistograms(const u32 max_read_length,
                                const std::vector<u64>& summary_stats_percentiles = kDefaultSummaryStatsPercentiles,
                                const DatasetMetadata& dataset_metadata = DatasetMetadata())
      : dataset_metadata(dataset_metadata),
        summary_stats_percentiles{summary_stats_percentiles},
        // Initialize histograms with the specified max depth
        post_filter_partial_read_length_histogram(dataset_metadata.has_read_type_info ? max_read_length + 1 : 0, 0),
        post_filter_full_read_length_histogram(dataset_metadata.has_read_type_info ? max_read_length + 1 : 0, 0),
        post_filter_all_read_length_histogram(max_read_length + 1, 0) {
  }

  bool exclude_softclips{};
  // Dataset metadata used to determine which rows to omit
  DatasetMetadata dataset_metadata;
  // Percentiles to be included in the summary
  std::vector<u64> summary_stats_percentiles;
  histogram::Histogram<u64> post_filter_partial_read_length_histogram;
  histogram::Histogram<u64> post_filter_full_read_length_histogram;
  histogram::Histogram<u64> post_filter_all_read_length_histogram;

  /**
   * Output the read length histograms to a TSV file.
   *
   * For datasets with read type information (i.e. whether a read is full length or partial),
   * the output will contain one column for the read length distribution for all reads,
   * one column for the partial read length distribution, and one column for the full read length distribution.
   */
  void WriteTsv(const fs::path& output_path, const io::Comments& comments) const;
  /**
   * Output a summary of the read length statistics to a TSV file.
   *
   * The summary statistics include: min, max, mean, median, and
   * a list of percentiles as specified by the user when creating the object.
   */
  void WriteSummaryTsv(const fs::path& output_path, const io::Comments& comments) const;

  ReadLengthHistograms& operator+=(const ReadLengthHistograms& other) {
    post_filter_partial_read_length_histogram += other.post_filter_partial_read_length_histogram;
    post_filter_full_read_length_histogram += other.post_filter_full_read_length_histogram;
    post_filter_all_read_length_histogram += other.post_filter_all_read_length_histogram;
    return *this;
  }
};

}  // namespace xoos::alignment_metrics
