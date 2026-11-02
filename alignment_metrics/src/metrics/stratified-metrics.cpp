#include "metrics/stratified-metrics.h"

#include <filesystem>
#include <utility>

#include <xoos/error/error.h>

namespace xoos::alignment_metrics {

StratifiedMetrics::StratifiedMetrics(const AlignmentMetricsOptions& options,
                                     const DatasetMetadata& dataset_metadata,
                                     const RegionLabelsSet& labels)
    : _labels_empty(labels.Empty()), _label_names(labels.Names()) {
  // Always hold at least one container: in non-stratified mode a single top-level container, else
  // one container per label in registration (bit-index) order.
  const size_t container_count = _labels_empty ? 1 : labels.Size();
  _metrics.reserve(container_count);
  for (size_t i = 0; i < container_count; ++i) {
    _metrics.emplace_back(options, dataset_metadata);
  }
}

StratifiedMetrics::StratifiedMetrics(Metrics metrics) : _labels_empty(true) {
  _metrics.emplace_back(std::move(metrics));
}

size_t StratifiedMetrics::Size() const {
  return _metrics.size();
}

Metrics& StratifiedMetrics::At(const size_t index) {
  return _metrics.at(index);
}

const Metrics& StratifiedMetrics::At(const size_t index) const {
  return _metrics.at(index);
}

StratifiedMetrics& StratifiedMetrics::operator+=(const StratifiedMetrics& other) {
  if (_metrics.size() != other._metrics.size()) {
    throw error::Error("Cannot merge StratifiedMetrics with {} containers into one with {} containers",
                       other._metrics.size(),
                       _metrics.size());
  }
  for (size_t i = 0; i < _metrics.size(); ++i) {
    _metrics.at(i) += other._metrics.at(i);
  }
  return *this;
}

void StratifiedMetrics::WriteMetrics(const fs::path& output_dir,
                                     const io::TsvMetadata& tsv_metadata,
                                     const std::string& output_prefix) const {
  if (_labels_empty) {
    _metrics.front().WriteMetrics(output_dir, tsv_metadata, output_prefix);
    return;
  }
  for (size_t i = 0; i < _metrics.size(); ++i) {
    _metrics.at(i).WriteMetrics(output_dir / _label_names.at(i), tsv_metadata, output_prefix);
  }
}

}  // namespace xoos::alignment_metrics
