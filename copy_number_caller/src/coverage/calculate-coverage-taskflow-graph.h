#pragma once
#include <string>
#include <utility>

#include <taskflow/core/graph.hpp>

#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "coverage/coverage.h"
#include "intervals/augmented-intervals.h"

namespace xoos::cnc {

// Default SAM flag mask excluding unmapped, secondary, QC-fail, and duplicate reads.
constexpr u32 kDefaultExcludeFlags = 1796;

class CalculateCoverageTaskflowGraph {
 public:
  CalculateCoverageTaskflowGraph(fs::path bam_file,
                                 AugmentedIntervals augmented_intervals,
                                 std::string exclude_flags,
                                 const bool ignore_dn)
      : _bam_file(std::move(bam_file)),
        _augmented_intervals(std::move(augmented_intervals)),
        _exclude_flags(std::move(exclude_flags)),
        _ignore_DN(ignore_dn),
        _graph(BuildGraph()) {
  }

  // NOLINTNEXTLINE
  tf::Graph& graph() {
    return _graph;
  }

  const CoverageRecords& GetResult() const {
    return _res;
  }

  CoverageRecords& GetResult() {
    return _res;
  }

 private:
  tf::Graph BuildGraph();
  fs::path _bam_file;
  AugmentedIntervals _augmented_intervals;
  std::string _exclude_flags;
  u32 _exclude_flags_int = kDefaultExcludeFlags;
  bool _ignore_DN;
  s32 _step = 0;
  CoverageRecords _res;
  tf::Graph _graph;
};
}  // namespace xoos::cnc
