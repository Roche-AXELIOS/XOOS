#pragma once

#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "segmentation/max-t.h"
#include "segmentation/segmentation-mode.h"

namespace xoos::cnc {

// Mode for the log-ratio segment undoing procedure.
// kDisabled: no segment undoing is performed; segments pass through untouched.
// kIterative: perform at least one merge pass at undo_logr_segments_sd_factor, then repeatedly relax the standard
//   deviation factor (by increment_undo_logr_segments_sd_factor) and re-run the merge pass until the segment count is
//   <= max_num_segments or no further merging occurs. The first pass always runs even when the input segment count is
//   already <= max_num_segments.
// kNonIterative: run a single merge pass at undo_logr_segments_sd_factor; max_num_segments and the increment
//   are ignored.
enum class SegmentUndoMode : u8 {
  kDisabled,
  kIterative,
  kNonIterative
};

const f64 kSegmentationDefaultMaxP = 0.01;
const size_t kSegmentationDefaultMinObsPerSegment = 2;
const size_t kSegmentationDefaultNPermutations = 10'000;
const size_t kSegmentationDefaultMaxNumSegments = 500;
const CbsMaxTMethod kSegmentationDefaultCbsMethod = CbsMaxTMethod::kBruteForce;
const f64 kSegmentationDefaultPruningClusteringParameter = 0.25;
const f64 kSegmentationDefaultBafDistanceWeight = 2.0;
// logR-neutral gate for BAF weighting: the baf_distance_weight is applied to the BAF-axis distance
// only when both segments in a pair have |mean_logr| below this threshold. This protects copy-neutral
// LOH (neutral logR, allelic imbalance) without distorting genuine gain/loss segments.
const f64 kSegmentationDefaultNeutralLogrThreshold = 0.05;
const f64 kSegmentationDefaultMinTForAutomaticSegmentation = 7;
const SegmentUndoMode kSegmentationDefaultUndoLogrMode = SegmentUndoMode::kDisabled;
const f64 kSegmentationDefaultUndoLogrSegmentsSdFactor = 3.0;
const f64 kSegmentationDefaultIncrementUndoLogrSegmentsSdFactor = 1.0;
const bool kSegmentationDefaultUndoDhSegments = false;
const f64 kSegmentationDefaultUndoDhSegmentsSdFactor = 3.0;
const f64 kSegmentationDefaultIncrementUndoDhSegmentsSdFactor = 1.0;

struct SegmentationOptions {
  f64 max_p = kSegmentationDefaultMaxP;
  size_t min_obs_per_segment = kSegmentationDefaultMinObsPerSegment;
  size_t n_permutations = kSegmentationDefaultNPermutations;
  bool disable_hierarchical_pruning = false;
  bool disable_merging = false;
  bool no_single_obvs_subsegments = true;
  bool truncate_outliers = false;
  CbsMaxTMethod cbs_method = kSegmentationDefaultCbsMethod;
  SegmentationMode segmentation_mode{SegmentationMode::kUnknown};
  f64 pruning_clustering_parameter = kSegmentationDefaultPruningClusteringParameter;
  f64 baf_distance_weight = kSegmentationDefaultBafDistanceWeight;
  f64 min_t_for_automatic_segmentation = kSegmentationDefaultMinTForAutomaticSegmentation;
  SegmentUndoMode undo_logr_mode = kSegmentationDefaultUndoLogrMode;
  f64 undo_logr_segments_sd_factor = kSegmentationDefaultUndoLogrSegmentsSdFactor;
  f64 increment_undo_logr_segments_sd_factor = kSegmentationDefaultIncrementUndoLogrSegmentsSdFactor;
  bool undo_dh_segments = kSegmentationDefaultUndoDhSegments;
  f64 undo_dh_segments_sd_factor = kSegmentationDefaultUndoDhSegmentsSdFactor;
  f64 increment_undo_dh_segments_sd_factor = kSegmentationDefaultIncrementUndoDhSegmentsSdFactor;
  size_t max_num_segments = kSegmentationDefaultMaxNumSegments;
};
}  // namespace xoos::cnc
