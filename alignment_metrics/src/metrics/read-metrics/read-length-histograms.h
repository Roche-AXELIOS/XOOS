#pragma once

#include <vector>

#include <xoos/histogram/histogram-summary.h>
#include <xoos/io/metadata-util.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "alignment-metrics-options.h"
#include "metadata/dataset-metadata.h"

namespace xoos::alignment_metrics {

/**
 * A struct to hold read length histograms for a dataset.
 *
 * Two variants of read length are tracked for every mapped read:
 * - including soft-clipped bases (the full query length), and
 * - excluding soft-clipped bases.
 *
 * For each variant, if the dataset contains read type information (i.e. whether a read is full length or partial),
 * it will contain the following histograms:
 * - partial_read_length_histogram: Histogram for partial reads.
 * - full_read_length_histogram: Histogram for full reads.
 * - all_read_length_histogram: Combined histogram for all reads.
 *
 * Otherwise, it will only contain the all_read_length_histogram for each variant.
 *
 * For each of the above, two sets of histograms are tracked:
 * - mapped_*: count all mapped reads (primary alignments) before read-level filters are applied.
 * - post_filter_*: count the subset of mapped reads that also passed the mapping quality and flag filters.
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
        // Initialize histograms with the specified max depth. Add 1 to the size to account for the 0-bin
        mapped_partial_read_length_histogram(dataset_metadata.has_read_type_info ? max_read_length + 1 : 0, 0),
        mapped_full_read_length_histogram(dataset_metadata.has_read_type_info ? max_read_length + 1 : 0, 0),
        mapped_all_read_length_histogram(max_read_length + 1, 0),
        mapped_partial_read_length_excluding_softclips_histogram(
            dataset_metadata.has_read_type_info ? max_read_length + 1 : 0, 0),
        mapped_full_read_length_excluding_softclips_histogram(
            dataset_metadata.has_read_type_info ? max_read_length + 1 : 0, 0),
        mapped_all_read_length_excluding_softclips_histogram(max_read_length + 1, 0),
        post_filter_partial_read_length_histogram(dataset_metadata.has_read_type_info ? max_read_length + 1 : 0, 0),
        post_filter_full_read_length_histogram(dataset_metadata.has_read_type_info ? max_read_length + 1 : 0, 0),
        post_filter_all_read_length_histogram(max_read_length + 1, 0),
        post_filter_partial_read_length_excluding_softclips_histogram(
            dataset_metadata.has_read_type_info ? max_read_length + 1 : 0, 0),
        post_filter_full_read_length_excluding_softclips_histogram(
            dataset_metadata.has_read_type_info ? max_read_length + 1 : 0, 0),
        post_filter_all_read_length_excluding_softclips_histogram(max_read_length + 1, 0) {
  }

  // Dataset metadata used to determine which rows to omit
  DatasetMetadata dataset_metadata;
  // Percentiles to be included in the summary
  std::vector<u64> summary_stats_percentiles;
  // Mapped reads (before read-level filters), read length including soft-clipped bases
  histogram::Histogram<u64> mapped_partial_read_length_histogram;
  histogram::Histogram<u64> mapped_full_read_length_histogram;
  histogram::Histogram<u64> mapped_all_read_length_histogram;
  // Mapped reads (before read-level filters), read length excluding soft-clipped bases
  histogram::Histogram<u64> mapped_partial_read_length_excluding_softclips_histogram;
  histogram::Histogram<u64> mapped_full_read_length_excluding_softclips_histogram;
  histogram::Histogram<u64> mapped_all_read_length_excluding_softclips_histogram;
  // Reads passing read-level filters, read length including soft-clipped bases
  histogram::Histogram<u64> post_filter_partial_read_length_histogram;
  histogram::Histogram<u64> post_filter_full_read_length_histogram;
  histogram::Histogram<u64> post_filter_all_read_length_histogram;
  // Reads passing read-level filters, read length excluding soft-clipped bases
  histogram::Histogram<u64> post_filter_partial_read_length_excluding_softclips_histogram;
  histogram::Histogram<u64> post_filter_full_read_length_excluding_softclips_histogram;
  histogram::Histogram<u64> post_filter_all_read_length_excluding_softclips_histogram;

  /**
   * Output the read length (including soft-clipped bases) histograms to a TSV file.
   *
   * For datasets with read type information (i.e. whether a read is full length or partial),
   * the output will contain one column for the read length distribution for all reads,
   * one column for the partial read length distribution, and one column for the full read length distribution.
   */
  void WriteTsv(const fs::path& output_path, const io::TsvMetadata& tsv_metadata) const;
  /**
   * Output the read length (excluding soft-clipped bases) histograms to a TSV file.
   *
   * See WriteTsv for the column layout.
   */
  void WriteExcludingSoftClipsTsv(const fs::path& output_path, const io::TsvMetadata& tsv_metadata) const;
  /**
   * Output a summary of the read length (including soft-clipped bases) statistics to a TSV file.
   *
   * The summary statistics include: min, max, mean, median, and
   * a list of percentiles as specified by the user when creating the object.
   */
  void WriteSummaryTsv(const fs::path& output_path, const io::TsvMetadata& tsv_metadata) const;
  /**
   * Output a summary of the read length (excluding soft-clipped bases) statistics to a TSV file.
   *
   * See WriteSummaryTsv for the row layout.
   */
  void WriteExcludingSoftClipsSummaryTsv(const fs::path& output_path, const io::TsvMetadata& tsv_metadata) const;

  ReadLengthHistograms& operator+=(const ReadLengthHistograms& other) {
    mapped_partial_read_length_histogram += other.mapped_partial_read_length_histogram;
    mapped_full_read_length_histogram += other.mapped_full_read_length_histogram;
    mapped_all_read_length_histogram += other.mapped_all_read_length_histogram;
    mapped_partial_read_length_excluding_softclips_histogram +=
        other.mapped_partial_read_length_excluding_softclips_histogram;
    mapped_full_read_length_excluding_softclips_histogram +=
        other.mapped_full_read_length_excluding_softclips_histogram;
    mapped_all_read_length_excluding_softclips_histogram += other.mapped_all_read_length_excluding_softclips_histogram;
    post_filter_partial_read_length_histogram += other.post_filter_partial_read_length_histogram;
    post_filter_full_read_length_histogram += other.post_filter_full_read_length_histogram;
    post_filter_all_read_length_histogram += other.post_filter_all_read_length_histogram;
    post_filter_partial_read_length_excluding_softclips_histogram +=
        other.post_filter_partial_read_length_excluding_softclips_histogram;
    post_filter_full_read_length_excluding_softclips_histogram +=
        other.post_filter_full_read_length_excluding_softclips_histogram;
    post_filter_all_read_length_excluding_softclips_histogram +=
        other.post_filter_all_read_length_excluding_softclips_histogram;
    return *this;
  }
};

}  // namespace xoos::alignment_metrics
