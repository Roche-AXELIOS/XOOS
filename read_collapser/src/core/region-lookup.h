#pragma once

#include <map>

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "core/region.h"

namespace xoos::read_collapser {

/**
 * Represents a zero-based half-open interval [start, end).
 */
struct Interval {
  s64 start{};
  s64 end{};

  bool Contains(s64 position) const;
  bool Overlaps(const Interval& other) const;
};

/**
 * An identifier that uniquely identifies a sub-region within a super region
 * specified by the super region ID and the sub-region index.
 */
struct SubRegionId {
  size_t super_region_id;
  size_t sub_region_index;
};

/**
 * A lookup table for efficiently finding the closest overlapping genomic regions.
 *
 * In read_collapser, we process reads by regions. However, a read can overlap with multiple regions.
 * To avoid processing the same read multiple times, we use a heuristic where a read is assigned
 * to the region whose start position is closest to the start position of the read.
 *
 * This data structure organizes genomic intervals by chromosome and provides fast lookup
 * capabilities to find the region whose start position is closest to a given query interval.
 */
class RegionLookupTable {
 private:
  /**
   * Sub-regions grouped by chromosome.
   */
  std::map<s32, vec<Interval>, std::less<>> _regions_by_tid;
  /**
   * Super region is a group of regions that are processed together.
   * It is used to group regions that are close to each other in the genome.
   *
   * This lookup vector is used to store the indices of the first sub-region in each super region.
   * Using this index, we can jump directly to any sub-region within a super region without iterating through all
   * regions in `regions_by_chr`.
   */
  vec<size_t> _super_region_offsets;

 public:
  /**
   * Constructs a RegionLookupTable from a vector of super regions.
   */
  explicit RegionLookupTable(const vec<SuperRegion>& super_regions);

  /**
   * Consider the range defined by `interval_to_search` on the contig identified by `tid`. Given a sub-region identified
   * by the super region ID and sub-region index that overlaps with this range, this function finds the sub-region whose
   * start position is the closest to the start position of the range among all sub-regions indexed in this
   * RegionLookupTable. `sub_region_id` is used as a hint to narrow down the search space.
   *
   * NOTE: This function assumes that the regions are sorted by their start position and that
   * the [position, position + max_distance) overlaps with the sub-region identified by the super region ID and sub
   * region index.
   */
  std::optional<Interval> FindClosestOverlappingRegion(s32 tid,
                                                       const Interval& interval_to_search,
                                                       const SubRegionId& sub_region_id) const;
};

}  // namespace xoos::read_collapser
