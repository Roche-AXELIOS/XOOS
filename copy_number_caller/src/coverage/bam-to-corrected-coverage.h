#pragma once
#include <optional>
#include <string>
#include <utility>

#include <taskflow/core/graph.hpp>

#include <xoos/io/metadata-util.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>

#include "coverage/coverage-check.h"
#include "coverage/coverage.h"
#include "intervals/augmented-intervals.h"

namespace xoos::cnc {
class BamToCorrectedCoverageTaskFlowGraph {
 public:
  BamToCorrectedCoverageTaskFlowGraph(fs::path bam_file,
                                      AugmentedIntervals augmented_intervals,
                                      const std::optional<fs::path>& coverage_out_fname,
                                      const std::optional<fs::path>& corrected_coverage_out_fname,
                                      std::string exclude_flags,
                                      const bool ignore_dn,
                                      const f64 gc_correct_first_span,
                                      io::CommandLineInfo command_line_info)
      : _bam_file(std::move(bam_file)),
        _augmented_intervals(std::move(augmented_intervals)),
        _coverage_out_fname(coverage_out_fname),
        _corrected_coverage_out_fname(corrected_coverage_out_fname),
        _exclude_flags(std::move(exclude_flags)),
        _ignore_dn(ignore_dn),
        _gc_correct_first_span(gc_correct_first_span),
        _command_line_info(std::move(command_line_info)),
        _graph(BuildGraph()) {
  }

  // NOLINTNEXTLINE
  tf::Graph& graph() {
    return _graph;
  }

  CoverageRecords& GetResult() {
    return _res;
  }

  bool IsLowCoverage() const {
    return cnc::IsLowCoverage(_coverage_check_result);
  }

  const std::optional<CoverageCheckResult>& GetCoverageCheckResult() const {
    return _coverage_check_result;
  }

 private:
  tf::Graph BuildGraph();
  // Runs GC correction on the coverage computed by the coverage task and populates _res. Extracted out of the
  // BuildGraph() lambda so the lambda itself stays within the taskflow node's expected size.
  void RunGcCorrect();
  fs::path _bam_file;
  AugmentedIntervals _augmented_intervals;
  CoverageRecords _uncorrected_coverage;
  CoverageRecords _res;
  std::optional<fs::path> _coverage_out_fname;
  std::optional<fs::path> _corrected_coverage_out_fname;
  std::string _exclude_flags;
  bool _ignore_dn;
  f64 _gc_correct_first_span;
  io::CommandLineInfo _command_line_info;
  tf::Graph _graph;
  std::optional<CoverageCheckResult> _coverage_check_result;
};
}  // namespace xoos::cnc
