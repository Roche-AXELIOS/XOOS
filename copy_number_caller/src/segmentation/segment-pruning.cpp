#include "segmentation/segment-pruning.h"

#include <fastcluster.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <optional>

#include <xoos/log/logging.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "segmentation/fastcluster-helper.h"
#include "segmentation/interval-trees.h"
#include "segmentation/segmentation-options.h"

namespace xoos::cnc {
/**
 * @brief breaks a set of logr genomicsegmnets according to the partitioning defined by dh_segments. for a
 set of dh segments overlapping a logr segment, the dh segments are extended such that the endpoints of each segment
 correspond to the LogR observation that surrounds the midpoint between the two BAF segments.
 * @param logr_segments
 * @param logrs
 * @param dh_segments
 * @param dhs
 * @return a list of genomicsegment the same size as dh_segments, with the boundaries extended to fully cover the
 genomic space of logr_segments

 */
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
vec<GenomicSegment> BreakLogrSegmentsByDHSegments(const vec<GenomicSegment>& logr_segments,
                                                  const Observations& logrs,
                                                  const vec<GenomicSegment>& dh_segments,
                                                  const Observations& /*dhs*/) {
  vec<GenomicSegment> new_segments = dh_segments;
  const IntervalTrees dh_trees(dh_segments);
  for (const auto& logr_seg : logr_segments) {
    // find all the dh segments that overlap this logr segment
    vec<size_t> dh_seg_idxs = dh_trees.LookUp(logr_seg.contig, logr_seg.start, logr_seg.end);
    // special handling when there are no or just one BAF segment overlapping logr_seg
    if (dh_seg_idxs.size() == 1) {
      new_segments[dh_seg_idxs[0]].start = logr_seg.start;
      new_segments[dh_seg_idxs[0]].end = logr_seg.end;
      continue;
    }
    if (dh_seg_idxs.empty()) {
      continue;
    }
    new_segments[dh_seg_idxs[0]].start = logr_seg.start;
    size_t left_idx = 0;
    size_t right_idx = 1;
    size_t left_seg_idx = dh_seg_idxs[left_idx];
    size_t right_seg_idx = dh_seg_idxs[right_idx];
    while (right_idx < dh_seg_idxs.size()) {
      // genomic position between the end of the left DH segment and the start of the right DH segment
      size_t midpoint = (new_segments[left_seg_idx].end + new_segments[right_seg_idx].start) / 2;
      // find the logr value in logr_seg that precedes midpoint
      for (size_t logr_idx = logr_seg.arr_start.value() + 1; logr_idx < logr_seg.arr_end.value(); ++logr_idx) {
        if (logrs.ends[static_cast<Eigen::Index>(logr_idx)] >= midpoint) {
          // If it's the case that the left segment falls entirely within the
          // observation, keep merging it with the right segment until the right
          // end is greater than the end of the observation
          if (new_segments[left_seg_idx].start >= logrs.starts[static_cast<Eigen::Index>(logr_idx)]) {
            // Recompute the boundary between the (possibly already-grown) left segment and the
            // current right segment each iteration and test it *before* merging. Using the
            // freshly computed midpoint in the loop guard prevents a single tiny DH segment
            // from triggering a runaway merge that swallows large downstream segments and
            // destroys a legitimate boundary (e.g. a CN-LOH -> balanced transition).
            while (right_idx < dh_seg_idxs.size()) {
              midpoint = (new_segments[left_seg_idx].end + new_segments[right_seg_idx].start) / 2;
              if (midpoint >= logrs.ends[static_cast<Eigen::Index>(logr_idx)]) {
                break;
              }
              new_segments[left_seg_idx].end = new_segments[right_seg_idx].end;
              new_segments[left_seg_idx].arr_end = new_segments[right_seg_idx].arr_end;
              new_segments[right_seg_idx].start = 0;
              new_segments[right_seg_idx].end = 0;
              new_segments[right_seg_idx].arr_start = 0;
              new_segments[right_seg_idx].arr_end = 0;
              ++right_idx;
              if (right_idx < dh_seg_idxs.size()) {
                right_seg_idx = dh_seg_idxs[right_idx];
              }
            }
            if (right_idx == dh_seg_idxs.size()) {
              break;
            }
            // Merging may have advanced the boundary well past the current observation, so
            // re-resolve the logr observation that now contains the boundary rather than
            // assuming it is the immediately following one.
            while (logr_idx < logr_seg.arr_end.value() && logrs.ends[static_cast<Eigen::Index>(logr_idx)] < midpoint) {
              ++logr_idx;
            }
            // if it's the last logr idx, then we extend the left segment to the end, invalidate the right segment, and
            // stop
            if (logr_idx >= logr_seg.arr_end.value()) {
              new_segments[left_seg_idx].end = logrs.ends[static_cast<Eigen::Index>(logr_seg.arr_end.value() - 1)];
              new_segments[left_seg_idx].arr_end = new_segments[right_seg_idx].arr_end;
              new_segments[right_seg_idx].start = 0;
              new_segments[right_seg_idx].end = 0;
              // Clear array indices too so the invalidated segment is fully reset,
              // consistent with the inner merge loop above and avoiding stale
              // array indices being observed before the final erase_if.
              new_segments[right_seg_idx].arr_start = 0;
              new_segments[right_seg_idx].arr_end = 0;
              // Signal exhaustion so the outer while loop exits rather than
              // advancing left_seg_idx to this now-invalidated segment and
              // restoring its end in the next iteration.
              right_idx = dh_seg_idxs.size();
              break;
            }
          }
          // end of left segment should be end of PREVIOUS logr observation
          new_segments[left_seg_idx].end = logrs.ends[static_cast<Eigen::Index>(logr_idx - 1)];
          // start of right segment should be start of this logr observation
          new_segments[right_seg_idx].start = logrs.starts[static_cast<Eigen::Index>(logr_idx)];
          break;
        }
      }
      // If the for loop broke because DH segments were exhausted (inner-while or last-logr-idx
      // case), right_idx has been set to dh_seg_idxs.size(). Exit before accessing tail code
      // which would advance left_seg_idx onto the just-invalidated segment.
      if (right_idx >= dh_seg_idxs.size()) {
        break;
      }
      left_idx = right_idx;
      left_seg_idx = dh_seg_idxs[left_idx];
      right_idx = left_idx + 1;
      if (right_idx >= dh_seg_idxs.size()) {
        break;
      }
      right_seg_idx = dh_seg_idxs[right_idx];
    }
  }
  // Drop degenerate segments. This catches both the (0,0) tombstones planted by the merge
  // branch above and zero/negative-length segments whose two boundaries snapped onto the same
  // contiguous log-ratio bin edge (e.g. chrX:119469007-119469007). A legitimate leading
  // segment starts at genomic position 0 but has end > 0, so start < end and it is naturally
  // exempt.
  std::erase_if(new_segments, [](const GenomicSegment& gs) { return gs.start >= gs.end; });
  // re-scan new segments and make sure that logr-segments are _fully covered_. If not, extend the segments to the logr
  // segments boundaries
  const IntervalTrees new_segment_trees(new_segments);
  for (const auto& logr_seg : logr_segments) {
    vec<size_t> new_seg_idxs = new_segment_trees.LookUp(logr_seg.contig, logr_seg.start, logr_seg.end);
    // inspect first and last new segments. Adjust their breakpoints if necessary
    if (!new_seg_idxs.empty()) {
      new_segments[new_seg_idxs.front()].start = std::min(new_segments[new_seg_idxs.front()].start, logr_seg.start);
      new_segments[new_seg_idxs.back()].end = std::max(new_segments[new_seg_idxs.back()].end, logr_seg.end);
    }
  }
  return new_segments;
}

/**
 * @brief Pruning segments by clustering segments with similar mean logr and weighted SNP AF if provided
 * @param logr_segments
 * @param joint_clustering bool, whether or not to perform join clustering
 * @param h height parameter in Hierarchical clustering
 * @param baf_weight multiplier applied to the mirrored-BAF distance between two segments, but only when both are
 *        copy-neutral (|mean_logr| < kSegmentationDefaultNeutralLogrThreshold). Values >1 increase the separation of
 *        copy-neutral allelic-imbalance segments (e.g. cnLOH) from balanced neighbours, making them less likely to be
 *        merged, while leaving genuine gain/loss segments unweighted. Only has an effect when joint_clustering is true.
 * @return a list of pruned genomicsegment
 */
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
vec<GenomicSegment> PruneSegmentsByClustering(const vec<GenomicSegment>& logr_segments,
                                              const bool joint_clustering,
                                              const std::optional<f64> h,
                                              const f64 baf_weight) {
  auto time1 = std::chrono::high_resolution_clock::now();
  constexpr s32 kMinSnps = 5;
  vec<GenomicSegment> new_segments;

  size_t n_segments = logr_segments.size();
  if (n_segments <= 1) {
    Logging::Info("Less than two segments found, skip Pruning...", n_segments);
    new_segments = logr_segments;
    return new_segments;
  }
  Logging::Info("Pruning {} segments with hierarchical clustering", n_segments);

  vec<bool> ignore_seg(n_segments, false);
  vec<Point> points;

  for (size_t k = 0; k < n_segments; ++k) {
    ignore_seg[k] = !logr_segments[k].num_obs.has_value() || logr_segments[k].num_obs == 0 ||
                    !logr_segments[k].mean_logr.has_value() || std::isnan(logr_segments[k].mean_logr.value());
    if (joint_clustering) {
      ignore_seg[k] = ignore_seg[k] && (!logr_segments[k].num_snps.has_value() || logr_segments[k].num_snps < kMinSnps);
    }
  }

  // Some segments have no logr, ignore these segments for clustering, but maintain an index to retrieve them
  vec<s32> original_idx(n_segments, -1);
  s32 k = 0;
  f64 mean_af = 0.0;
  for (size_t i = 0; i < n_segments; ++i) {
    if (logr_segments[i].mean_logr.has_value() && !std::isnan(logr_segments[i].mean_logr.value())) {
      if (joint_clustering && logr_segments[i].mean_dh.has_value() && !std::isnan(logr_segments[i].mean_dh.value())) {
        constexpr f64 kDhToMirroredAf = 2.0;
        // mirrored AF = 1/2 DH, stored raw. The baf_weight is applied pairwise inside GetDistance (gated on both
        // segments being logR-neutral), so it must NOT be baked into the coordinate here.
        mean_af = logr_segments[i].mean_dh.value() / kDhToMirroredAf;
      }
      points.emplace_back(logr_segments[i].mean_logr.value(), mean_af);
      original_idx[i] = k++;
    }
  }

  constexpr f64 kDefaultBandwidth = 0.25;
  f64 h_new = kDefaultBandwidth;
  if (h.has_value()) {
    h_new = h.value();
  }
  Logging::Debug("Setting clustering paramter H = {}", h_new);

  vec<s32> clusters = Hclust(
      points, hclust_fast_methods::HCLUST_METHOD_COMPLETE, h_new, baf_weight, kSegmentationDefaultNeutralLogrThreshold);

  Logging::Debug("Merging segments based on clustering");

  // Setting segment mean logr as weighted mean of segments in the same cluster, weighted by num of snps
  vec<f64> logr_mean(clusters.size(), 0.0);
  vec<s32> total_num_snp(clusters.size(), 0);
  vec<f64> logr_mean_sum(clusters.size(), 0.0);
  for (size_t i = 0; i < n_segments; ++i) {
    if (ignore_seg[i] || original_idx[i] < 0) {
      continue;
    }
    const s32 cluster_id = clusters[original_idx[i]];
    total_num_snp[cluster_id] += static_cast<s32>(logr_segments[i].num_obs.value());
    logr_mean_sum[cluster_id] +=
        logr_segments[i].mean_logr.value() * static_cast<s32>(logr_segments[i].num_obs.value());
  }

  for (size_t i = 0; i < total_num_snp.size(); ++i) {
    if (total_num_snp[i] > 0) {
      logr_mean[i] = logr_mean_sum[i] / total_num_snp[i];
    }
  }

  // Mering consecutive segments with the same cluster id
  GenomicSegment gs = logr_segments[0];
  for (size_t i = 1; i < n_segments; ++i) {
    if (ignore_seg[i - 1] || ignore_seg[i] || logr_segments[i].contig != logr_segments[i - 1].contig ||
        original_idx[i] < 0 || original_idx[i - 1] < 0 || clusters[original_idx[i]] != clusters[original_idx[i - 1]]) {
      new_segments.emplace_back(gs);
      gs = logr_segments[i];
    } else {
      gs.end = logr_segments[i].end;
      gs.arr_end = logr_segments[i].arr_end;
      gs.num_obs = gs.num_obs.value() + logr_segments[i].num_obs.value();
      gs.mean_logr = logr_mean[clusters[original_idx[i]]];
      if (gs.num_snps.has_value() && logr_segments[i].num_snps.has_value()) {
        gs.num_snps = gs.num_snps.value() + logr_segments[i].num_snps.value();
      }
      // mean_dh is invalidated by merging; clear it so callers re-populate from raw observations
      gs.mean_dh = std::nullopt;
    }
  }
  new_segments.emplace_back(gs);
  const auto time2 = std::chrono::high_resolution_clock::now();
  auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count();
  Logging::Info("Found {} segments after Pruning in {}ms", new_segments.size(), total_time);

  return new_segments;
}

vec<GenomicSegment> PruneSegmentsByClusteringPerSeedSegment(const vec<GenomicSegment>& seed_segments,
                                                            const vec<GenomicSegment>& segments_to_prune,
                                                            const bool joint_clustering,
                                                            const std::optional<f64> h,
                                                            const f64 baf_weight) {
  vec<GenomicSegment> res;
  const IntervalTrees seg_trees(segments_to_prune);
  for (const auto& seed_seg : seed_segments) {
    const vec<size_t> seg_idxs = seg_trees.LookUp(seed_seg.contig, seed_seg.start, seed_seg.end);
    if (!seg_idxs.empty()) {
      vec<GenomicSegment> segs;
      segs.reserve(seg_idxs.size());
      for (const auto& idx : seg_idxs) {
        segs.emplace_back(segments_to_prune[idx]);
      }
      vec<GenomicSegment> pruned_segs = PruneSegmentsByClustering(segs, joint_clustering, h, baf_weight);
      res.insert(res.end(), pruned_segs.begin(), pruned_segs.end());
    } else {
      Logging::Debug(
          "No segments found for seed segment: {}:{}-{}. Skipping pruning for this one and removing this seed segment "
          "from the output...",
          seed_seg.contig,
          seed_seg.start,
          seed_seg.end);
    }
  }
  return res;
}
}  // namespace xoos::cnc
