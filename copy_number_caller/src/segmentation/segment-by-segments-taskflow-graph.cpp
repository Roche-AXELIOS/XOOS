#include "segmentation/segment-by-segments-taskflow-graph.h"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <stdexcept>

#include <taskflow/algorithm/for_each.hpp>
#include <taskflow/core/flow_builder.hpp>
#include <taskflow/core/graph.hpp>
#include <taskflow/core/task.hpp>

#include <xoos/log/logging.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "segmentation/cbs-wrapper.h"
#include "segmentation/genomic-segments.h"

namespace xoos::cnc {
/**
 * @brief build the Taskflow Graph for SegmentBySegmentsTaskflowGraph.
 * init -> [[CBS segmentation per seed segment ]] -> merge subsegments
 * @return
 */
tf::Graph SegmentBySegmentsTaskflowGraph::BuildGraph() {
  _options.no_single_obvs_subsegments = _no_single_obvs_subsegments;
  tf::Graph graph;
  tf::FlowBuilder builder(graph);
  tf::Task init = builder.emplace([this]() {
    if (_obvs.starts.size() == 0 || _obvs.obvs.size() == 0) {
      Logging::Error("observations are empty!");
      throw std::runtime_error("observations are empty!");
    }
    if (_seed_segments.empty()) {
      Logging::Error("seed segments are empty!");
      throw std::runtime_error("seed segments are empty!");
    }
    Logging::Info("Segmenting on {} seed segments and {} observations", _seed_segments.size(), _obvs.obvs.size());
    _obvs_tree = IntervalTrees(_obvs);
    _buffer = vec<vec<GenomicSegment>>(_seed_segments.size());
    _n_seed_segments = static_cast<s32>(_seed_segments.size());
  });
  // we have to use std::ref here to make sure _n_seed_segments isn't 0 (which it might be when this function is called)
  tf::Task pf = builder.for_each_index(0, std::ref(_n_seed_segments), 1, [this](const size_t seed_seg_idx) {
    const auto& seed_seg = _seed_segments[seed_seg_idx];
    const vec<GenomicSegment> genomic_segments = SegmentByOneSegment(seed_seg);
    _buffer[seed_seg_idx] = genomic_segments;
  });
  const tf::Task merge = builder.emplace([this]() {
    for (const auto& subvec : _buffer) {
      for (const auto& gseg : subvec) {
        _res.push_back(gseg);
      }
    }
    Logging::Info("Found {} segments", _res.size());
  });
  init.precede(pf);
  pf.precede(merge);
  return graph;
}

/**
 * @brief Perform CBS segmentation on a "seed" segment -> the new segments will be subsequences of this seed segment
 * @param seed_seg - GenomicSegment from which to begin CBS
 * @param obvs input observations
 * @param obvs_tree IntervalTrees index of obvs
 * @param max_p maximum p-value for calling a segment significant (default: 0.01)
 * @param n_permutations  number of permutations for significance testing
 * @param min_obs_per_segment minimum number of observations for a segment
 * @return list of GenomicSegments
 */
vec<GenomicSegment> SegmentBySegmentsTaskflowGraph::SegmentByOneSegment(const GenomicSegment& seed_seg) {
  vec<GenomicSegment> genomic_segments;
  vec<size_t> obvs_idxs = _obvs_tree.LookUp(seed_seg.contig, seed_seg.start, seed_seg.end);
  // return original segment if there are no observations overlapping it
  if (obvs_idxs.empty()) {
    Logging::Info(
        "No observations found for seed segment: {}:{}-{}", seed_seg.contig, seed_seg.start + 1, seed_seg.end);
    if (_keep_seed_segments) {
      genomic_segments.push_back(seed_seg);
    }
  } else {
    const size_t start = *(std::ranges::min_element(obvs_idxs));
    const size_t end = *(std::ranges::max_element(obvs_idxs)) + 1;
    const Observations seg_obvs = _obvs.FilterByRange(start, end);
    Logging::Debug("finding segments in Segment {}:{}-{}", seed_seg.contig, seed_seg.start + 1, seed_seg.end);
    auto segments = CircularBinarySegmentationWrapper(seg_obvs.obvs, _options);
    for (auto& s : segments) {
      s.start += start;
      s.end += start;
    }
    genomic_segments = ConvertSegmentsToGenomicSegments(_obvs, segments);
    // adjust the start and ends of the subsegments to match the seed segment. this is done in case the first
    // observation or the last observation also overlap the previous or next segment, resp.
    const size_t gidx = genomic_segments.size() - segments.size();
    genomic_segments[gidx].start = std::max(genomic_segments[gidx].start, seed_seg.start);
    genomic_segments[genomic_segments.size() - 1].end =
        std::min(genomic_segments[genomic_segments.size() - 1].end, seed_seg.end);
    for (auto& gseg : genomic_segments) {
      gseg.in_allosome = seed_seg.in_allosome;
      gseg.in_pseudo_autosomal_region = seed_seg.in_pseudo_autosomal_region;
    }
  }
  return genomic_segments;
}

}  // namespace xoos::cnc
