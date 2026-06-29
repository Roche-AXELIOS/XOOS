#include "core/region-lookup.h"

#include <optional>
#include <utility>

#include <xoos/error/error.h>
#include <xoos/io/bed-region.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/malloc-ptr.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>

namespace xoos::read_collapser {

bool Interval::Contains(const s64 position) const {
  return position >= start && position < end;
}

bool Interval::Overlaps(const Interval& other) const {
  return start < other.end && end > other.start;
}

RegionLookupTable::RegionLookupTable(const vec<SuperRegion>& super_regions) {
  for (const auto& super_region : super_regions) {
    _super_region_offsets.push_back(_regions_by_tid[super_region.front().tid].size());
    for (const auto& sub_region : super_region) {
      _regions_by_tid[sub_region.tid].emplace_back(Interval{sub_region.start, sub_region.end});
    }
  }
}

std::optional<Interval> RegionLookupTable::FindClosestOverlappingRegion(const s32 tid,
                                                                        const Interval& interval_to_search,
                                                                        const SubRegionId& sub_region_id) const {
  // Validate inputs
  const size_t super_region_id = sub_region_id.super_region_id;
  const size_t sub_region_idx = sub_region_id.sub_region_index;
  const auto it = _regions_by_tid.find(tid);
  // No regions for this chromosome or super region ID is out of bounds
  if (it == _regions_by_tid.end() || super_region_id >= _super_region_offsets.size()) {
    return std::nullopt;
  }
  const auto& intervals = it->second;
  size_t curr_offset = _super_region_offsets.at(super_region_id) + sub_region_idx;
  // Start offset is out of bounds
  if (curr_offset >= intervals.size()) {
    return std::nullopt;
  }

  // Initialize search state
  const s64 position = interval_to_search.start;
  s64 closest_distance = std::abs(position - intervals.at(curr_offset).start);
  std::optional<Interval> closest_interval = intervals.at(curr_offset).Overlaps(interval_to_search)
                                                 ? std::make_optional(intervals.at(curr_offset))
                                                 : std::nullopt;
  const bool search_forward = std::cmp_greater_equal(position, intervals.at(curr_offset).start);
  const s64 direction = search_forward ? 1 : -1;

  // Search for the closest overlapping interval in the determined direction
  // Stop if we reach the beginning or end of the intervals
  while ((direction > 0 && curr_offset < intervals.size() - 1) || (direction < 0 && curr_offset > 0)) {
    curr_offset += direction;
    const s64 distance = std::abs(position - intervals.at(curr_offset).start);

    // If the distance is already greater than the closest distance found, we can stop searching
    if (distance > closest_distance) {
      break;
    }

    // Update the closest interval if we find a better one
    if (distance < closest_distance && intervals.at(curr_offset).Overlaps(interval_to_search)) {
      closest_distance = distance;
      closest_interval = intervals.at(curr_offset);
    } else if (distance == closest_distance && intervals.at(curr_offset).Contains(position)) {
      // If the distance is the same, then we only update if the interval contains the position we are searching for.
      // This serves as a tie-breaker in case we have multiple intervals at the same distance
      // to the position we are searching for
      closest_interval = intervals.at(curr_offset);
    }
  }
  return closest_interval;
}

}  // namespace xoos::read_collapser
