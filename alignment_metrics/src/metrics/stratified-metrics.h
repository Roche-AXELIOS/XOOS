#pragma once

#include <concepts>
#include <string>

#include <xoos/io/metadata-util.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "alignment-metrics-options.h"
#include "core/region-labels-set.h"
#include "metadata/dataset-metadata.h"
#include "metrics/metrics.h"

namespace xoos::alignment_metrics {

/**
 * @brief Holds one Metrics container per region label, stratifying output by label.
 *
 * In non-stratified mode (the region labels set is empty) this holds a single Metrics container and
 * writes to the top-level output directory, preserving the original output layout. In stratified
 * mode it holds one container per label (in label bit-index order) and writes each to a
 * per-label subdirectory `<output_dir>/<label>/{accuracy,coverage,read_metrics}`.
 *
 * A genomic position or read may belong to several labels at once; callers route each unit of work
 * to every container whose label bit is set (see `ForEachContainer`).
 */
class StratifiedMetrics {
 public:
  StratifiedMetrics(const AlignmentMetricsOptions& options,
                    const DatasetMetadata& dataset_metadata,
                    const RegionLabelsSet& labels);

  /**
   * Construct a non-stratified StratifiedMetrics wrapping a single existing Metrics container. Mainly a
   * convenience for tests that build and populate a Metrics directly.
   */
  explicit StratifiedMetrics(Metrics metrics);

  /**
   * Number of Metrics containers held (label count, or 1 in non-stratified mode).
   */
  size_t Size() const;

  /**
   * Access the container for the given label bit index.
   */
  Metrics& At(size_t index);
  const Metrics& At(size_t index) const;

  /**
   * Invoke `fn(Metrics&)` for every container whose label bit is set in `label_mask`. In
   * non-stratified mode the single container is always visited regardless of `label_mask`.
   */
  template <std::invocable<Metrics&> Fn>
  void ForEachContainer(const u16 label_mask, const Fn& fn) {
    if (_labels_empty) {
      fn(_metrics.front());
      return;
    }
    for (size_t index = 0; index < _metrics.size(); ++index) {
      if ((label_mask & (1U << index)) != 0) {
        fn(_metrics.at(index));
      }
    }
  }

  /**
   * Invoke `fn(Metrics&)` for every container. Used for work that must be attributed to all labels
   * (e.g. unplaced unmapped reads).
   */
  template <std::invocable<Metrics&> Fn>
  void ForAllContainers(const Fn& fn) {
    for (auto& metrics : _metrics) {
      fn(metrics);
    }
  }

  /**
   * Merge another StratifiedMetrics element-wise into this one (per-thread aggregation).
   */
  StratifiedMetrics& operator+=(const StratifiedMetrics& other);

  /**
   * Writes all containers to disk. In non-stratified mode writes directly to `output_dir`;
   * in stratified mode writes each container to `output_dir/<label_name>/`.
   */
  void WriteMetrics(const fs::path& output_dir,
                    const io::TsvMetadata& tsv_metadata,
                    const std::string& output_prefix) const;

 private:
  // True when built from an empty region labels set (non-stratified single-container mode).
  bool _labels_empty;
  // Label names in bit-index order; empty in non-stratified mode.
  vec<std::string> _label_names;
  // One container per label (or exactly one in non-stratified mode).
  vec<Metrics> _metrics;
};

}  // namespace xoos::alignment_metrics
