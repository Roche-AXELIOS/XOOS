#include "run-alignment-metrics.h"

#include "alignment-metrics-options.h"
#include "alignment-metrics.h"

namespace xoos::alignment_metrics {

// Runner for alignment metrics for the given options.
void RunAlignmentMetrics(const AlignmentMetricsOptions& opts) {
  AlignmentMetrics alignment_metrics(opts);
  alignment_metrics.CalculateMetrics();
}
}  // namespace xoos::alignment_metrics
