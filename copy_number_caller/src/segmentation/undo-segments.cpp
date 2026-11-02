#include "segmentation/undo-segments.h"

#include <algorithm>
#include <cassert>

#include <xoos/log/logging.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "interval-trees.h"
#include "util/eigen-math.h"

namespace xoos::cnc {

constexpr f64 kMedianAbsoluteDeviationScaleFactor = 1.4826;
constexpr size_t kUndoSegmentsMinNumMedians = 2;
constexpr f64 kSqrtBase = 2.0;
constexpr f64 kSqrtExponent = 0.5;

/**
 * @brief Calculate the Median Absolute Deviation (MAD) of a vector of observations
 */
namespace {
f64 MedianAbsoluteDeviation(const Eigen::VectorXd& obvs) {
  // median absolute deviation (MAD) is defined as the median of the absolute deviations from the median
  const f64 median = eigen::Median(obvs);
  const Eigen::VectorXd abs_dev = (obvs.array() - median).abs().matrix();
  return kMedianAbsoluteDeviationScaleFactor * eigen::Median(abs_dev);
}
}  // namespace

/**
 * @brief Calculate the trimmed standard deviation of a vector of observations
 */
f64 GetTrimmedStandardDeviation(const Eigen::VectorXd& arr) {
  // Estimate the standard deviation from the median absolute deviation of successive differences,
  // scaled by 1/sqrt(2) to correct for the differencing.
  const auto n = static_cast<u64>(arr.size());
  const Eigen::VectorXd diff = arr.segment(1, n - 1) - arr.segment(0, n - 1);
  f64 trimmed_sd = MedianAbsoluteDeviation(diff);
  trimmed_sd = trimmed_sd / std::pow(kSqrtBase, kSqrtExponent);
  return trimmed_sd;
}

/**
 * @brief holds segment with its associated observations, for each of iteration
 */
struct SegmentWithObservations {
  GenomicSegment segment;
  Observations obvs;
};

Observations MergeObservations(const Observations& obvs1, const Observations& obvs2) {
  Observations merged_obvs;
  merged_obvs.contigs.reserve(obvs1.contigs.size() + obvs2.contigs.size());
  merged_obvs.contigs.insert(merged_obvs.contigs.begin(), obvs1.contigs.begin(), obvs1.contigs.end());
  merged_obvs.contigs.insert(merged_obvs.contigs.end(), obvs2.contigs.begin(), obvs2.contigs.end());
  merged_obvs.starts = eigen::JoinCols(obvs1.starts, obvs2.starts);
  merged_obvs.ends = eigen::JoinCols(obvs1.ends, obvs2.ends);
  merged_obvs.obvs = eigen::JoinCols(obvs1.obvs, obvs2.obvs);
  return merged_obvs;
}

/**
 * @brief Undo genomic segments by merging adjacent segments if their median observation difference is less than a
 * factor (max_trimmed_sd_factor) of the trimmed standard deviation (trimmed_sd)
 * @param segments input segments with associated observations
 * @param trimmed_sd trimmed standard deviation of the observations (use GetTrimmedStandardDeviation to calculate)
 * @param max_trimmed_sd_factor maximum factor of trimmed standard deviation to consider for merging segments
 * @return merged segments
 */
namespace {
vec<SegmentWithObservations> UndoGenomicSegments(const vec<SegmentWithObservations>& segments,
                                                 f64 trimmed_sd,
                                                 f64 max_trimmed_sd_factor) {
  if (segments.empty()) {
    return {};
  }
  const f64 max_threshold = trimmed_sd * max_trimmed_sd_factor;
  // merging will be done in a series of iterations. In each iteration, we will transfer segments from
  // `prev_merged_segments` to `merged_segments`, merging segments in `prev_merged_segments` if their median observation
  // difference is less than the threshold, and then compare the size of `merged_segments` with that of
  // `prev_merged_segments` to determine whether to continue iterating (i.e. whether any merging has been done in this
  // iteration)
  vec<SegmentWithObservations> prev_merged_segments(segments);
  vec<SegmentWithObservations> merged_segments;
  merged_segments.reserve(prev_merged_segments.size());
  vec<f64> medians;
  medians.reserve(segments.size());
  vec<u8> breakpoints_to_merge(segments.size());
  while (merged_segments.size() != prev_merged_segments.size()) {
    merged_segments.clear();
    medians.clear();
    // calculate the medians of all the segments
    for (const auto& seg : prev_merged_segments) {
      medians.emplace_back(eigen::Median(seg.obvs.obvs));
    }
    // this vector will tell us which segments to merge together.  0 means don't merge, 1 means merge
    breakpoints_to_merge.clear();
    breakpoints_to_merge.resize(medians.size());
    // breakpoints[i] == 1 implies we need to merge segment[i] with segment[i-1]. breakpoints[i] == 0 means do
    // not merge segment[i] with segment[i-1].
    std::ranges::fill(breakpoints_to_merge, 0);
    // for each element in medians, compare to the previous and next element
    for (size_t i = 1; i < medians.size() - 1; ++i) {
      const f64 prev_diff = std::abs(medians[i] - medians[i - 1]);
      const f64 next_diff = std::abs(medians[i] - medians[i + 1]);
      // make sure the distance to this neighbor is less than the difference to the next neighbor. Cannot merge both
      // neighbors in same pass
      if (prev_diff <= next_diff && prev_diff <= max_threshold) {
        // if it is 1, leave it as 1. If it is 0, change to 1
        breakpoints_to_merge[i] += 1 - breakpoints_to_merge[i];
      } else if (next_diff < prev_diff && next_diff <= max_threshold) {
        breakpoints_to_merge[i + 1] += 1 - breakpoints_to_merge[i + 1];
      } else {
        // if neither neighbor is less than the threshold, then we don't merge this segment with either neighbor, so do
        // nothing and leave the value as is
      }
    }
    if (medians.size() >= kUndoSegmentsMinNumMedians) {
      // handle the last breakpoint
      const f64 prev_diff = std::abs(medians[medians.size() - 1] - medians[medians.size() - 2]);
      if (prev_diff <= max_threshold) {
        // if it is 1, leave it as 1. If it is 0, change to 1
        breakpoints_to_merge[medians.size() - 1] += 1 - breakpoints_to_merge[medians.size() - 1];
      }
    }
    // begin merge process. We iteratively merge segments with their left
    // first segment automatically added to merged_segments because it can't be merged to the left.
    merged_segments.emplace_back(prev_merged_segments[0]);
    // walk through the breakpoints and medians, and simultaneously merge segments from prev_merged_segments and add to
    // the back of merged_segments if breakpoint is 1, otherwise just add the segment to merged_segments without merging
    for (size_t i = 1; i < medians.size(); ++i) {
      if (breakpoints_to_merge[i] > 0) {
        merged_segments.back().segment.end = prev_merged_segments[i].segment.end;
        // The merged segment now spans through the absorbed segment, so its observation index range
        // must be extended too. Adjacent segments are contiguous in the global observation array, so
        // the merged segment's arr_end is the absorbed segment's arr_end. Failing to update this
        // leaves arr_end stale (pointing at the end of only the first sub-segment), which later
        // truncates downstream consumers that scan observations using arr_start/arr_end (e.g.
        // BreakLogrSegmentsByDHSegments), dropping coverage and destroying real segment boundaries.
        merged_segments.back().segment.arr_end = prev_merged_segments[i].segment.arr_end;
        // also transfer all the observations of the to-be-merged segment
        merged_segments.back().obvs = MergeObservations(merged_segments.back().obvs, prev_merged_segments[i].obvs);
      } else {
        merged_segments.emplace_back(prev_merged_segments[i]);
      }
    }
    // in the next round, we will merge the segments even further (if possible), so swap merged_segments with
    // prev_merged_segments. merged_segments will be cleared at the beginning of the next loop
    std::swap(prev_merged_segments, merged_segments);
  }
  return prev_merged_segments;
}
}  // namespace

/**
 * @brief wrapper around UndoSegments to work with GenomicSegment and Observations
 */
vec<GenomicSegment> UndoSegments(const vec<GenomicSegment>& segments,
                                 const Observations& obvs,
                                 const f64 trimmed_sd,
                                 const f64 max_trimmed_sd_factor) {
  // convert segments and obvs into SegmentWithObservations
  vec<SegmentWithObservations> segs_with_obvs;
  const IntervalTrees obvs_trees(obvs);
  for (const auto& seg : segments) {
    const Observations seg_obvs = obvs.FilterByIdxs(obvs_trees.LookUp(seg.contig, seg.start, seg.end));
    segs_with_obvs.emplace_back(SegmentWithObservations{.segment = seg, .obvs = seg_obvs});
  }
  // call UndoGenomicSegments
  const vec<SegmentWithObservations> undone_segs =
      UndoGenomicSegments(segs_with_obvs, trimmed_sd, max_trimmed_sd_factor);
  // extract GenomicSegments from undone_segs
  vec<GenomicSegment> ret;
  ret.reserve(undone_segs.size());
  for (const auto& seg_with_obvs : undone_segs) {
    ret.emplace_back(seg_with_obvs.segment);
  }
  return ret;
}

/**
 * @brief Undo segments in a single, non-iterative merge pass at a fixed max_trimmed_sd_factor. Unlike
 * IterativeUndoGenomicSegments, the factor is not escalated and no target segment count is applied. Note that except
 * for the contig, start, and end fields of each GenomicSegment member, all other fields are invalidated.
 * @param seed_segments seed segments
 * @param segments input segments
 * @param obvs observations vector
 * @param max_trimmed_sd_factor max trimmed sd factor used for the single merge pass
 * @return merged segments
 */
vec<GenomicSegment> NonIterativeUndoGenomicSegments(const vec<GenomicSegment>& seed_segments,
                                                    const vec<GenomicSegment>& segments,
                                                    const Observations& obvs,
                                                    f64 max_trimmed_sd_factor) {
  if (segments.empty()) {
    return segments;
  }
  Logging::Info("Non-iteratively undoing segments at max_trimmed_sd_factor = {} ({} input segments)",
                max_trimmed_sd_factor,
                segments.size());
  // associate segments with seed segments, and observations with segments
  const IntervalTrees segment_trees(segments);
  const IntervalTrees obvs_trees(obvs);
  vec<vec<SegmentWithObservations>> segments_per_seed(seed_segments.size());
  for (size_t i = 0; i < seed_segments.size(); ++i) {
    auto segment_idxs = segment_trees.LookUp(seed_segments[i].contig, seed_segments[i].start, seed_segments[i].end);
    for (auto segment_idx : segment_idxs) {
      const GenomicSegment& seg = segments[segment_idx];
      segments_per_seed[i].emplace_back(seg, obvs.FilterByIdxs(obvs_trees.LookUp(seg.contig, seg.start, seg.end)));
    }
  }
  // undo segments under each seed once at the fixed factor
  const f64 trimmed_sd = GetTrimmedStandardDeviation(obvs.obvs);
  vec<SegmentWithObservations> all_undone_segments;
  for (const auto& seg_vec : segments_per_seed) {
    if (!seg_vec.empty()) {
      vec<SegmentWithObservations> undone_segments = UndoGenomicSegments(seg_vec, trimmed_sd, max_trimmed_sd_factor);
      all_undone_segments.insert(all_undone_segments.end(), undone_segments.begin(), undone_segments.end());
    }
  }
  Logging::Info("Non-iterative undoing completed: {}->{} segments", segments.size(), all_undone_segments.size());
  vec<GenomicSegment> ret;
  ret.reserve(all_undone_segments.size());
  for (const auto& seg_with_obvs : all_undone_segments) {
    ret.emplace_back(seg_with_obvs.segment);
  }
  return ret;
}

/**
 * @brief Iteratively undo segments until the number of segments is less than or equal to max_num_segments. At least one
 * merge pass is always performed at start_max_trimmed_sd_factor before the max_num_segments target is evaluated, so
 * over-segmented input is still merged even when the input count is already at/below the target. Each subsequent round
 * adds increment_max_trimmed_sd_factor to the factor (a constant per-round ramp) and stops once the target is reached
 * or no further segments can be undone. Note that except for the contig, start, end, arr_start, and arr_end fields of
 * each GenomicSegment member, all other fields in these members are invalidated
 * @param seed_segments seed segments
 * @param segments input segments
 * @param obvs observations vector (the trimmed standard deviation is computed internally from these observations)
 * @param start_max_trimmed_sd_factor starting max trimmed sd factor, used for the first merge pass
 * @param increment_max_trimmed_sd_factor amount added to the max trimmed sd factor at each subsequent round
 * @param max_num_segments maximum number of segments to achieve
 * @return merged segments
 */
vec<GenomicSegment> IterativeUndoGenomicSegments(const vec<GenomicSegment>& seed_segments,
                                                 const vec<GenomicSegment>& segments,
                                                 const Observations& obvs,
                                                 const f64 start_max_trimmed_sd_factor,
                                                 const f64 increment_max_trimmed_sd_factor,
                                                 const size_t max_num_segments) {
  if (segments.empty()) {
    Logging::Info("No segments to undo");
    return segments;
  }
  Logging::Info("Iteratively undoing segments. Targeting {} -> {} segments", segments.size(), max_num_segments);

  // associate segments with seed segments, and observations with segments
  const IntervalTrees segment_trees(segments);
  const IntervalTrees obvs_trees(obvs);
  vec<vec<SegmentWithObservations>> segments_per_seed(seed_segments.size());
  for (size_t i = 0; i < seed_segments.size(); ++i) {
    const auto segment_idxs =
        segment_trees.LookUp(seed_segments[i].contig, seed_segments[i].start, seed_segments[i].end);
    for (const auto segment_idx : segment_idxs) {
      const GenomicSegment& seg = segments[segment_idx];
      segments_per_seed[i].emplace_back(seg, obvs.FilterByIdxs(obvs_trees.LookUp(seg.contig, seg.start, seg.end)));
    }
  }
  // prepare for undoing segments
  vec<SegmentWithObservations> all_undone_segments;
  for (const auto& vec : segments_per_seed) {
    all_undone_segments.insert(all_undone_segments.end(), vec.begin(), vec.end());
  }
  f64 max_trimmed_sd_factor = start_max_trimmed_sd_factor;
  size_t iter = 0;
  const f64 trimmed_sd = GetTrimmedStandardDeviation(obvs.obvs);
  // Always perform at least one merge pass at the initial factor (iter == 0, so no increment is applied) before
  // evaluating the target segment count. Each subsequent round adds a single increment_max_trimmed_sd_factor to the
  // factor (a constant per-round ramp) and continues only while the target has not been reached and progress is still
  // being made (fewer segments than the previous round).
  while (true) {
    const size_t prev_num_segments = all_undone_segments.size();
    all_undone_segments.clear();
    if (iter > 0) {
      max_trimmed_sd_factor += increment_max_trimmed_sd_factor;
    }
    for (const auto& seg_vec : segments_per_seed) {
      if (!seg_vec.empty()) {
        vec<SegmentWithObservations> undone_segments = UndoGenomicSegments(seg_vec, trimmed_sd, max_trimmed_sd_factor);
        all_undone_segments.insert(all_undone_segments.end(), undone_segments.begin(), undone_segments.end());
      }
    }
    Logging::Info("Iteration {}: max_trimmed_sd_factor = {}, number of segments = {}",
                  iter,
                  max_trimmed_sd_factor,
                  all_undone_segments.size());
    ++iter;
    // Stop once the target is reached (segments <= max_num_segments) or no further progress was made this round.
    if (all_undone_segments.size() <= max_num_segments || all_undone_segments.size() >= prev_num_segments) {
      break;
    }
  }
  Logging::Info("Iterative undoing completed in {} iterations: {}->{} segments",
                iter,
                segments.size(),
                all_undone_segments.size());
  vec<GenomicSegment> ret;
  ret.reserve(max_num_segments);
  for (const auto& seg_with_obvs : all_undone_segments) {
    ret.emplace_back(seg_with_obvs.segment);
  }
  return ret;
}
}  // namespace xoos::cnc
