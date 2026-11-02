#include "coverage/bam-to-corrected-coverage.h"

#include <fstream>
#include <utility>

#include <taskflow/taskflow.hpp>

#include <xoos/log/logging.h>

#include "coverage/calculate-coverage-taskflow-graph.h"
#include "coverage/coverage-check.h"
#include "gc-correct/gc-correct-algorithm.h"

namespace xoos::cnc {
tf::Graph BamToCorrectedCoverageTaskFlowGraph::BuildGraph() {
  tf::Graph graph;
  tf::FlowBuilder builder(graph);
  tf::Task coverage_task = builder.emplace([this](tf::Subflow& subflow) {
    CalculateCoverageTaskflowGraph calculate_coverage_tf_graph(
        _bam_file, _augmented_intervals, _exclude_flags, _ignore_dn);
    subflow.composed_of(calculate_coverage_tf_graph);
    subflow.join();
    _uncorrected_coverage = calculate_coverage_tf_graph.GetResult();
    if (_coverage_out_fname.has_value()) {
      Logging::Info("Writing coverage to {}", _coverage_out_fname->string());
      std::ofstream ofs(_coverage_out_fname.value());
      _uncorrected_coverage.Write(ofs, _command_line_info);
    }
    // Check median coverage before proceeding to GC correction
    _coverage_check_result = CheckMedianCoverage(_uncorrected_coverage);
    if (!_coverage_check_result->is_sufficient) {
      Logging::Error(
          "Median coverage ({:.2f}x) is below minimum threshold ({:.0f}x). "
          "Skipping copy number calling.",
          _coverage_check_result->median_coverage,
          kMinMedianCoverageThreshold);
    }
  });
  const tf::Task gc_correct_task = builder.emplace([this]() { RunGcCorrect(); });
  coverage_task.precede(gc_correct_task);
  return graph;
}

void BamToCorrectedCoverageTaskFlowGraph::RunGcCorrect() {
  if (IsLowCoverage()) {
    Logging::Info("Skipping GC correction due to low coverage");
    return;
  }
  const GCCorrectResult gc_correct_res = GCCorrect(_uncorrected_coverage, _augmented_intervals, _gc_correct_first_span);
  _res.region = std::move(_uncorrected_coverage.region);
  _res.count = gc_correct_res.count;
  _res.total_coverage = gc_correct_res.total_coverage;
  _res.average_coverage = gc_correct_res.average_coverage;
  _res.on_target = std::move(_uncorrected_coverage.on_target);
  _res.mean_mapping_quality = std::move(_uncorrected_coverage.mean_mapping_quality);
  _res.has_rescued_secondaries = _uncorrected_coverage.has_rescued_secondaries.load();
  _res.UpdateRegionMap();
  if (_corrected_coverage_out_fname.has_value()) {
    Logging::Info("Writing corrected coverage to {}", _corrected_coverage_out_fname->string());
    std::ofstream ofs(_corrected_coverage_out_fname.value());
    _res.Write(ofs, _command_line_info);
  }
}
}  // namespace xoos::cnc
