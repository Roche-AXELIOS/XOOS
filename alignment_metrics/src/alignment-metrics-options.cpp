#include "alignment-metrics-options.h"

namespace xoos::alignment_metrics {

bool NeedHpMasking(const AlignmentMetricsOptions& options) {
  return !options.disable_hp_quality_modification &&
         (options.metric_types.has_accuracy_metrics || options.metric_types.has_coverage_metrics);
}

bool NeedBaseTypeDecoding(const AlignmentMetricsOptions& options) {
  return !options.disable_base_type_decoding &&
         (options.metric_types.has_accuracy_metrics || options.metric_types.has_coverage_metrics);
}

}  // namespace xoos::alignment_metrics
