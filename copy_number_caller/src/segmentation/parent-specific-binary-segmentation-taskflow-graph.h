#pragma once

#include <string>

#include <taskflow/core/graph.hpp>

#include <xoos/sex_predict/sex.h>
#include <xoos/types/vec.h>

#include "cli/sample-metadata-options.h"
#include "observations/observations.h"
#include "segmentation/genomic-segments.h"
#include "segmentation/segment-by-segments-taskflow-graph.h"
#include "segmentation/segmentation-mode.h"

namespace xoos::cnc {

class ParentSpecificBinarySegmentationTaskflowGraph {
 public:
  /**
   * @brief Taskflow graph that will run parallel implementation of the Parent-Specific-Binary-Segmentation algorithm
   * described in Olshen et. al 2011 (https://doi.org/10.1093/bioinformatics/btr329). Briefly, this algorithm performs
   * Circular Binary Segmentation (CBS) first on the tumor-normal logr values, then feeds the resulting segments into a
   * second CBS step based on B-allele fractions at germline heterozygous variants
   * @param logrs Observations representing logrs
   * @param dh_vals Observations representing dh_vals
   * @param seed_segments (optional) std::vector of GenomicSegments representing seed segments. This function will
   * parallelize on these seed segments. If nullopt, then parallelization will happen over chromosomes
   * @param max_p maximum p-value for calling a segment significant (default: 0.01)
   * @param n_permutations  number of permutations for significance testing
   * @param min_obs_per_segment minimum number of observations for a segment
   * NOTE: in a Taskflow Graph, we will run into segfaults if any shared data has a lifetime that does not span the
   * whole of the taskflow pipeline. Therefore, we have to we capture pass-by-reference values that this taskflow graph
   * needs from the parent taskflow as a "dependency injection via constructor"
   * To get the "result" from this taskflow graph, call GetResult()
   */
  ParentSpecificBinarySegmentationTaskflowGraph(const Observations& logrs,
                                                const Observations& dh_vals,
                                                const vec<GenomicSegment>& seed_segments,
                                                const SegmentationOptions& segmentation_options,
                                                const SampleMetadataOptions& sample_metadata_options)
      : _logrs(logrs),
        _dh_vals(dh_vals),
        _seed_segments(LoadSeedSegments(seed_segments,
                                        sample_metadata_options.sex.value_or(sex_predict::Sex::kUnknown),
                                        segmentation_options.segmentation_mode)),
        _options(segmentation_options),
        _seg_by_seg_graph(_seed_segments, _logrs, _options, false, true),
        _baf_seg_graph(_logr_segments, _dh_vals, _options, true, false),
        _graph(BuildGraph()) {
  }

  ParentSpecificBinarySegmentationTaskflowGraph(ParentSpecificBinarySegmentationTaskflowGraph&&) noexcept = default;
  ParentSpecificBinarySegmentationTaskflowGraph(const ParentSpecificBinarySegmentationTaskflowGraph&) noexcept = delete;
  ParentSpecificBinarySegmentationTaskflowGraph& operator=(ParentSpecificBinarySegmentationTaskflowGraph&&) noexcept =
      delete;
  ParentSpecificBinarySegmentationTaskflowGraph& operator=(
      const ParentSpecificBinarySegmentationTaskflowGraph&) noexcept = delete;
  ~ParentSpecificBinarySegmentationTaskflowGraph() = default;

  // NOLINTNEXTLINE
  tf::Graph& graph() {
    return _graph;
  }

  vec<GenomicSegment>& GetLogRSegments() {
    return _logr_segments;
  }

  // BAF segments before hierarchical pruning. Empty when no BAF source was provided.
  vec<GenomicSegment>& GetBafSegments() {
    return _dh_segments;
  }

  vec<GenomicSegment>& GetResult() {
    return _res;
  }

 private:
  tf::Graph BuildGraph();
  static vec<GenomicSegment> LoadSeedSegments(const vec<GenomicSegment>& seed_segments,
                                              sex_predict::Sex sex,
                                              SegmentationMode mode);
  vec<GenomicSegment> _logr_segments;
  vec<GenomicSegment> _dh_segments;
  // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members) - taskflow lifetime injection: these bind to
  // caller-owned inputs that must outlive the running graph, and are passed by reference into the BAF sub-graph below.
  const Observations& _logrs;
  const Observations& _dh_vals;
  const vec<GenomicSegment> _seed_segments;
  // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)
  vec<std::string> _contigs;
  vec<GenomicSegment> _res;
  SegmentationOptions _options;
  SegmentBySegmentsTaskflowGraph _seg_by_seg_graph;
  SegmentBySegmentsTaskflowGraph _baf_seg_graph;
  tf::Graph _graph;
};

}  // namespace xoos::cnc
