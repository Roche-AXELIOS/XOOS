#pragma once

#include <optional>

#include "alignment-metrics-options.h"
#include "metadata/dataset-metadata.h"
#include "metrics/accuracy-metrics/base-level-accuracy-summary.h"
#include "metrics/accuracy-metrics/errors-by-cluster-size.h"
#include "metrics/accuracy-metrics/errors-by-read-type.h"
#include "metrics/accuracy-metrics/errors-by-substitution-type.h"
#include "metrics/accuracy-metrics/hp-stat.h"
#include "metrics/accuracy-metrics/qscore-stats.h"

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

}  // namespace xoos::alignment_metrics
