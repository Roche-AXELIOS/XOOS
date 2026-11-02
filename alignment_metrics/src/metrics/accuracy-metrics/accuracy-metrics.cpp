#include "metrics/accuracy-metrics/accuracy-metrics.h"

#include <optional>

#include "alignment-metrics-options.h"
#include "metadata/dataset-metadata.h"

namespace xoos::alignment_metrics {

AccuracyMetrics::AccuracyMetrics(const AlignmentMetricsOptions& options, const DatasetMetadata& dataset_metadata)
    : errors_by_cluster_size(dataset_metadata.has_cluster_info
                                 ? std::make_optional(ErrorsByClusterSize(options.max_cluster_size_bin))
                                 : std::nullopt),
      errors_by_substitution_type(ErrorsBySubstitutionType(dataset_metadata)),
      qscore_stats(QscoreStats(dataset_metadata)),
      errors_by_read_type((dataset_metadata.has_read_type_info || dataset_metadata.has_strand_info)
                              ? std::make_optional(ErrorsByReadType(dataset_metadata))
                              : std::nullopt),
      hp_stats(options.calculate_hp_metrics ? std::make_optional(HpStats(options.max_hp_length)) : std::nullopt) {
}

AccuracyMetrics& AccuracyMetrics::operator+=(const AccuracyMetrics& other) {
  base_level_accuracy_summary += other.base_level_accuracy_summary;
  if (errors_by_cluster_size && other.errors_by_cluster_size) {
    errors_by_cluster_size->operator+=(*other.errors_by_cluster_size);
  }
  errors_by_substitution_type += other.errors_by_substitution_type;
  qscore_stats += other.qscore_stats;
  if (errors_by_read_type && other.errors_by_read_type) {
    errors_by_read_type->operator+=(*other.errors_by_read_type);
  }
  if (hp_stats && other.hp_stats) {
    hp_stats->operator+=(*other.hp_stats);
  }
  return *this;
}

}  // namespace xoos::alignment_metrics
