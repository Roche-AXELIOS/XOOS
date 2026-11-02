#pragma once

#include <optional>
#include <string>

#include <xoos/io/metadata-util.h>
#include <xoos/types/fs.h>

#include "alignment-metrics-options.h"
#include "metadata/dataset-metadata.h"
#include "metrics/accuracy-metrics/accuracy-metrics.h"
#include "metrics/coverage-metrics/coverage-metrics.h"
#include "metrics/coverage-metrics/gc-bias-metrics.h"
#include "metrics/read-metrics/read-metrics.h"

namespace xoos::alignment_metrics {

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
  std::optional<GcBiasMetrics> gc_bias_metrics;

  explicit Metrics(const AlignmentMetricsOptions& options, const DatasetMetadata& dataset_metadata);

  /**
   * Writes all populated metric categories to TSV files in organized subdirectories
   * (accuracy/, coverage/, read_metrics/) under the given output directory.
   * When `output_prefix` is non-empty, it is prepended to every output filename with an
   * underscore separator (e.g. prefix "SAMPLE" produces "SAMPLE_coverage_stats.tsv").
   *
   * @param output_dir Root directory for output files
   * @param tsv_metadata Command-line and sample metadata written as TSV comment headers
   * @param output_prefix Filename prefix; empty string means no prefix
   */
  void WriteMetrics(const fs::path& output_dir,
                    const io::TsvMetadata& tsv_metadata,
                    const std::string& output_prefix) const;

  Metrics& operator+=(const Metrics& other);
};
}  // namespace xoos::alignment_metrics
