#pragma once

#include <memory>
#include <vector>

#include "alignment-metrics-options.h"
#include "metadata/dataset-metadata.h"
#include "metrics/metrics.h"

namespace xoos::alignment_metrics {

// Creates the output directory for the metrics files.
void CreateOutputDir(const AlignmentMetricsOptions& options);

// Class responsible for running all metrics calculations and writing the metric results to files.
class AlignmentMetrics {
 public:
  explicit AlignmentMetrics(const AlignmentMetricsOptions& opts);

  // Runs the metrics calculations and writes the results to files.
  void CalculateMetrics();

  // Merges the metrics from each worker thread into a single set of metrics on this local thread.
  void MergeMetrics(const std::vector<std::shared_ptr<Metrics>>& per_thread_metrics);

  // Input options
  const AlignmentMetricsOptions& options;

  // Dataset metadata that is used to determine what metrics to calculate
  // since some metrics are only applicable to certain datasets
  DatasetMetadata dataset_metadata;

  // Finals metrics aggreagted over all regions that is used to compute the metrics
  Metrics final_metrics;
};
}  // namespace xoos::alignment_metrics
