#include "observations/filter-observations.h"

#include <Eigen/Dense>
#include <cstddef>

#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "observations/observations.h"
#include "segmentation/interval-trees.h"
#include "util/eigen-math.h"

namespace xoos::cnc {

/**
 * @brief Get observations from `right` that overlap with `left`
 * @param left Observations object to filter from
 * @param right Observations object to filter with
 * @return Observations object containing observations from right that overlap with left. Sets IsBAFSegObvs
 * automatically
 */
Observations GetOverlappingObservations(const Observations& left, const Observations& right) {
  if (right.contigs.empty()) {
    return {};
  }
  Observations ret;
  vec<u64> starts;
  vec<u64> ends;
  vec<f64> obvs;
  const vec<s32> dps;
  const IntervalTrees left_trees(left);
  for (size_t i = 0; i < right.contigs.size(); ++i) {
    const auto ei = static_cast<Eigen::Index>(i);
    const vec<size_t> idxs = left_trees.LookUp(right.contigs[i], right.starts[ei], right.ends[ei]);
    // we don't want to keep the observation if there are no overlapping regions
    if (idxs.empty()) {
      continue;
    }
    ret.contigs.push_back(right.contigs[i]);
    starts.push_back(right.starts[ei]);
    ends.push_back(right.ends[ei]);
    obvs.push_back(right.obvs[ei]);
    if (!right.dps.empty()) {
      ret.dps.push_back(right.dps[i]);
    }
  }
  ret.starts = eigen::FromStdVecU64(starts);
  ret.ends = eigen::FromStdVecU64(ends);
  ret.obvs = eigen::FromStdVec(obvs);
  return ret;
}

/**
 * @brief Get observations from `right` that overlap with `left`
 * @param left list of region strings to filter from
 * @param right Observations object to filter with
 * @return Observations object containing observations from right that overlap with left. Does not set IsBAFSegObvs
 * automatically
 */
Observations GetOverlappingObservations(const vec<std::string>& left, const Observations& right) {
  if (right.contigs.empty()) {
    return {};
  }
  Observations ret;
  vec<u64> starts;
  vec<u64> ends;
  vec<f64> obvs;
  const vec<s32> dps;
  const IntervalTrees left_trees(left);
  for (size_t i = 0; i < right.contigs.size(); ++i) {
    const auto ei = static_cast<Eigen::Index>(i);
    const vec<size_t> idxs = left_trees.LookUp(right.contigs[i], right.starts[ei], right.ends[ei]);
    // we don't want to keep the observation if there are no overlapping regions
    if (idxs.empty()) {
      continue;
    }  // keep the observation from `right`
    ret.contigs.push_back(right.contigs[i]);
    starts.push_back(right.starts[ei]);
    ends.push_back(right.ends[ei]);
    obvs.push_back(right.obvs[ei]);
    if (!right.dps.empty()) {
      ret.dps.push_back(right.dps[i]);
    }
  }
  ret.starts = eigen::FromStdVecU64(starts);
  ret.ends = eigen::FromStdVecU64(ends);
  ret.obvs = eigen::FromStdVec(obvs);
  return ret;
}
}  // namespace xoos::cnc
