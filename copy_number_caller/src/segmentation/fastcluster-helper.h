#pragma once
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "fastcluster.h"

namespace xoos::cnc {

/**
 * @brief a point in segment-clustering space: logR on one axis, mirrored BAF (raw dh/2) on the other.
 * The BAF value is stored unweighted; any weighting is applied pairwise in GetDistance so it can be
 * gated on the logR-neutral status of both points (see GetDistance).
 */
struct Point {
  f64 logr;
  f64 baf;

  Point(const f64 logr_val, const f64 baf_val) : logr(logr_val), baf(baf_val) {
  }
};

/**
 * @brief Euclidean distance between two points with logR-neutral-gated BAF weighting.
 * The BAF-axis difference is scaled by baf_weight only when both points are copy-neutral
 * (|logr| < neutral_logr_threshold). This protects copy-neutral LOH pairs (neutral logR, differing
 * allelic imbalance) from being merged, while leaving genuine gain/loss segments unweighted.
 * Symmetric in p1/p2: the gate predicate and the weighted difference are both order-independent.
 * @param baf_weight multiplier for the BAF difference when both points are neutral (1.0 = no-op).
 * @param neutral_logr_threshold |logr| cutoff below which a point is considered copy-neutral.
 */
static inline f64 GetDistance(const Point& p1,
                              const Point& p2,
                              const f64 baf_weight,
                              const f64 neutral_logr_threshold) {
  const f64 dx = p1.logr - p2.logr;
  f64 dy = p1.baf - p2.baf;
  const bool both_neutral =
      (std::abs(p1.logr) < neutral_logr_threshold) && (std::abs(p2.logr) < neutral_logr_threshold);
  if (both_neutral) {
    dy *= baf_weight;
  }
  return std::sqrt((dx * dx) + (dy * dy));
}

/**
 * @brief Wrapper function of the hclust_fast
 * @param baf_weight multiplier applied to the BAF-axis distance for logR-neutral pairs (see GetDistance).
 * @param neutral_logr_threshold |logr| cutoff defining copy-neutral points for the BAF-weight gate.
 * @return Vector of cluster id from 0, 1, ..., ncluster - 1
 */
inline vec<s32> Hclust(const vec<Point>& points,
                       const s32 hclust_method,
                       const f64 h,
                       const f64 baf_weight,
                       const f64 neutral_logr_threshold) {
  auto npoints = static_cast<s32>(points.size());

  vec<f64> distmat;
  vec<s32> merge;
  vec<f64> height;
  vec<s32> labels;

  distmat.resize(npoints * (npoints - 1) / 2);
  merge.resize(2 * static_cast<std::size_t>(npoints - 1));
  height.resize(npoints - 1);
  labels.resize(npoints);

  // Compute the upper triangle of the distance matrix excluding diagnol values
  s32 k = 0;
  for (s32 i = 0; i < npoints; ++i) {
    for (s32 j = i + 1; j < npoints; ++j) {
      distmat[k++] = GetDistance(points[i], points[j], baf_weight, neutral_logr_threshold);
    }
  }

  if (npoints > 1) {
    hclust_fast(npoints, distmat.data(), hclust_method, merge.data(), height.data());
    cutree_cdist(npoints, merge.data(), height.data(), h, labels.data());
    return labels;
  }
  if (npoints == 1) {
    return {0};
  }
  throw std::runtime_error("hclust: npoints is 0!");
}

}  // namespace xoos::cnc
