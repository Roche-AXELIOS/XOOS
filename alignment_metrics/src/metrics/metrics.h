#pragma once

#include <htslib/sam.h>

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/metadata-util.h>
#include <xoos/types/int.h>

#include "alignment-metrics-options.h"
#include "core/alignment.h"
#include "io/alignment-reader.h"
#include "metadata/alignment-metadata.h"
#include "metadata/dataset-metadata.h"
#include "metrics/accuracy-metrics/base-level-accuracy-summary.h"
#include "metrics/accuracy-metrics/errors-by-cluster-size.h"
#include "metrics/accuracy-metrics/errors-by-read-type.h"
#include "metrics/accuracy-metrics/errors-by-substitution-type.h"
#include "metrics/accuracy-metrics/hp-stat.h"
#include "metrics/accuracy-metrics/qscore-stats.h"
#include "metrics/coverage-metrics/coverage-histograms.h"
#include "metrics/coverage-metrics/coverage-uniformity-metrics.h"
#include "metrics/read-metrics/read-counts-by-cluster-size.h"
#include "metrics/read-metrics/read-length-histograms.h"
#include "metrics/read-metrics/read-metrics-summary.h"
#include "metrics/read-metrics/read-metrics-te.h"

namespace xoos::alignment_metrics {

/**
 * @brief Container for all sequencing accuracy-related metrics.
 *
 * Aggregates error statistics, quality score analysis, and substitution patterns across the dataset.
 * Includes base-level accuracy summaries, error stratification by cluster size and read type,
 * quality score calibration data, and homopolymer-specific error profiles. Optional metrics are
 * populated based on dataset characteristics (e.g., cluster info, HP analysis enabled).
 */
class AccuracyMetrics {
 public:
  AccuracyMetrics() = default;
  explicit AccuracyMetrics(const AlignmentMetricsOptions& options, const DatasetMetadata& dataset_metadata);

  BaseLevelAccuracySummary base_level_accuracy_summary;
  std::optional<ErrorsByClusterSize> errors_by_cluster_size;
  ErrorsBySubstitutionType errors_by_substitution_type;
  QscoreStats qscore_stats;
  std::optional<ErrorsByReadType> errors_by_read_type;
  std::optional<HpStats> hp_stats;

  AccuracyMetrics& operator+=(const AccuracyMetrics& other);
};

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
   * If the number of bases is not provided, it increments the 0-bin by 1.
   */
  void AddEmptyHistogramData(u32 num_bases = 1);
  /**
   * Update the coverage histograms with the given depth values
   */
  void AddHistogramData(u32 all_depth, u32 concordant_duplex_depth, u32 filtered_depth, bool is_part_of_hp = false);
  CoverageMetrics& operator+=(const CoverageMetrics& other);
};

/**
 * @brief Container for all read-level and alignment statistics.
 *
 * Aggregates read counts, alignment flags, mapping statistics, and read length distributions to provide comprehensive
 * summaries of sequencing run quality and read characteristics.
 */
class ReadMetrics {
  friend class Metrics;  // Allow Metrics to access private members for updating TE and unmapped metrics
 public:
  ReadMetrics() = default;
  explicit ReadMetrics(const AlignmentMetricsOptions& options, const DatasetMetadata& dataset_metadata);

  ReadMetricsSummary read_metrics_summary;
  ReadLengthHistograms read_length_histograms;
  std::optional<ReadCountsByClusterSize> read_counts_by_cluster_size;
  std::optional<TargetEnrichmentReadMetrics> te_read_metrics;

  /**
   * Update read metrics using the given alignment. We assume that the calling code has already checked that the
   * alignment is not double counted
   *
   * This updates various read counts, read length distributions, and base counts except `aligned_bases` and
   * `soft_clipped_bases` which are counted when processing the CIGAR ops of the alignment.
   *
   * @param alignment Alignment to process; used to extract BAM flags and read length
   * @param alignment_metadata Metadata for the alignment containing information such as full/partial status
   * @param passed_filter Whether the read passed the minimum mapping quality and BAM flag filters
   */
  void AddRead(const Alignment& alignment, const AlignmentMetadata& alignment_metadata, bool passed_filter);
  /**
   * Updates the unmapped read metrics in the final metrics object
   * by querying the BAM index for the number of unplaced unmapped reads (no coordinate).
   * Placed unmapped reads (reads with the unmapped flag 0x4 but has a coordinate) are
   * counted separately when querying the BAM file by regions.
   *
   * The metric object is updated in place. This function should only be called once
   * for the final metrics object (after merging per-thread metric objects) to avoid
   * double counting.
   */
  void UpdateUnmappedReadMetrics(const AlignmentReader& alignment_reader, bool passed_filter, u16 trim);
  ReadMetrics& operator+=(const ReadMetrics& other);

 private:
  bool _unmapped_metrics_updated = false;
};

/**
 * @brief Top-level container for all alignment metrics categories.
 *
 * Aggregates accuracy, coverage, and read metrics into a single object for unified processing.
 * Each metric category is optional and populated based on user-specified options and dataset
 * characteristics.
 */
class Metrics {
 public:
  std::optional<AccuracyMetrics> accuracy_metrics;
  std::optional<CoverageMetrics> coverage_metrics;
  std::optional<ReadMetrics> read_metrics;

  explicit Metrics(const AlignmentMetricsOptions& options, const DatasetMetadata& dataset_metadata);

  /**
   * Write the metrics contained in this object to the specified output directory with the given comments
   * as the TSV header.
   */
  void WriteMetrics(const fs::path& output_dir, const io::Comments& comments) const;

  Metrics& operator+=(const Metrics& other);
};
}  // namespace xoos::alignment_metrics
