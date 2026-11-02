#pragma once

#include <taskflow/core/graph.hpp>

#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "observations/observations.h"
#include "segmentation/genomic-segments.h"
#include "segmentation/interval-trees.h"
#include "segmentation/segmentation-options.h"

namespace xoos::cnc {

class SegmentBySegmentsTaskflowGraph {
 public:
  SegmentBySegmentsTaskflowGraph() = delete;

  /**
   * @brief create a Taskflow Graph for performing CBS segmentation, parallelizing on seed segments
   * @param seed_segments a list of known seed segments (a CBS procedure per seed segment). If empty, then this pipeline
   * will return an empty result
   * @param obvs input observations
   * @param max_p maximum p-value for calling a segment significant (default: 0.01)
   * @param n_permutations  number of permutations for significance testing
   * @param min_obs_per_segment minimum number of observations for a segment
   * NOTE: in a Taskflow Graph, we will run into segfaults if any shared data has a lifetime that does not span the
   * whole of the taskflow pipeline. Therefore, we have to we capture pass-by-reference values that this taskflow graph
   * needs from the parent taskflow as a "dependency injection via constructor"
   * To get the "result" from this taskflow graph, call GetResult()
   */
  SegmentBySegmentsTaskflowGraph(const vec<GenomicSegment>& seed_segments,
                                 const Observations& obvs,
                                 SegmentationOptions options,
                                 bool keep_seed_segments,
                                 bool no_single_obvs_subsegments)
      : _seed_segments(seed_segments),
        _obvs(obvs),
        _options(options),
        _keep_seed_segments(keep_seed_segments),
        _no_single_obvs_subsegments(no_single_obvs_subsegments),
        _graph(BuildGraph()) {
  }

  SegmentBySegmentsTaskflowGraph(const SegmentBySegmentsTaskflowGraph& rhs) = delete;
  SegmentBySegmentsTaskflowGraph(SegmentBySegmentsTaskflowGraph&& rhs) noexcept = default;
  SegmentBySegmentsTaskflowGraph& operator=(const SegmentBySegmentsTaskflowGraph& rhs) = delete;
  SegmentBySegmentsTaskflowGraph& operator=(SegmentBySegmentsTaskflowGraph&& rhs) = delete;
  ~SegmentBySegmentsTaskflowGraph() = default;

  // NOLINTNEXTLINE
  tf::Graph& graph() {
    return _graph;
  }

  const vec<GenomicSegment>& GetResult() {
    return _res;
  }

 private:
  vec<GenomicSegment> SegmentByOneSegment(const GenomicSegment& seed_seg);
  tf::Graph BuildGraph();
  vec<vec<GenomicSegment>> _buffer;
  // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members) - deferred binding: _seed_segments may alias the
  // parent graph's _logr_segments, which is populated at runtime; a value member would copy it while still empty.
  const vec<GenomicSegment>& _seed_segments;
  const Observations& _obvs;
  // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)
  IntervalTrees _obvs_tree;
  vec<GenomicSegment> _res;
  s32 _n_seed_segments = 0;
  SegmentationOptions _options;
  bool _keep_seed_segments = false;
  bool _no_single_obvs_subsegments = true;
  tf::Graph _graph;
};

}  // namespace xoos::cnc
