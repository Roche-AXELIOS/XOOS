#pragma once

#include <memory>
#include <optional>
#include <vector>

#include "alignment-metrics-options.h"
#include "metadata/dataset-metadata.h"
#include "metrics/stratified-metrics.h"

namespace xoos::alignment_metrics {

// Creates the output directory for the metrics files.
void CreateOutputDir(const AlignmentMetricsOptions& options);

// Class responsible for running all metrics calculations and writing the metric results to files.
class AlignmentMetrics {
 public:
  explicit AlignmentMetrics(const AlignmentMetricsOptions& opts);

  // Runs the metrics calculations and writes the results to files.
  void CalculateMetrics();

  // Merges the per-label metrics from each worker thread into a single set of metrics on this local thread.
  void MergeMetrics(const std::vector<std::shared_ptr<StratifiedMetrics>>& per_thread_metrics);

  // Non-owning reference to options owned by the caller for this object's lifetime; stored by
  // reference by design, making the runner intentionally non-assignable.
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-const-or-ref-data-members)
  const AlignmentMetricsOptions& options;

  // Dataset metadata that is used to determine what metrics to calculate
  // since some metrics are only applicable to certain datasets
  DatasetMetadata dataset_metadata;

  // Final per-label metrics aggregated over all regions. Constructed in CalculateMetrics once the
  // region labels set is known (from the RegionLookupTable), hence optional.
  std::optional<StratifiedMetrics> final_metrics;
};
}  // namespace xoos::alignment_metrics
