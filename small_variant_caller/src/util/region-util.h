#pragma once

#include <optional>
#include <string>

#include <xoos/io/bed-region.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>

namespace xoos::svc {

/**
 * Synopsis:
 * This file contains utility functions and data structures for handling genomic regions,
 * intervals, and target regions. It provides functionality for checking overlaps, merging
 * intervals, and reading/writing BED files.
 */

using BedRegion = io::BedRegion;

// An interval on a chromosome, defined by a start and end position
struct Interval {
  u64 start{};  // 0-based start coordinate
  u64 end{};    // end exclusive
  auto operator<=>(const Interval&) const = default;
};

// A partitioned region to compute features or filter variants in parallel
struct TargetRegion {
  std::string chrom;
  u64 start{};  // 0-based inclusive start position on the reference
  u64 end{};    // 0-based exclusive end position on the reference
};

using ChromIntervalsMap = StrMap<vec<Interval>>;
bool IntervalOverlapOrContained(const Interval& a, const Interval& b);
bool IntervalOverlap(u64 start1, u64 end1, u64 start2, u64 end2);
bool IntervalOverlap(const Interval& i, u64 start, u64 end);
bool IntervalOverlap(const Interval& a, const Interval& b);
bool IntervalOverlap(const vec<Interval>& sorted_intervals, u64 start, u64 end);
vec<Interval> MergeIntervals(vec<Interval>& sorted_intervals, u64 collapse_dist = 0);
std::optional<StrMap<vec<Interval>>> GetChromIntervalMap(const std::optional<fs::path>& region_bed);
void WriteBed(const ChromIntervalsMap& chrom_map, const fs::path& out_bed_path);
std::optional<vec<BedRegion>> GetBedRegions(const std::optional<fs::path>& region_bed);

/**
 * @brief Convert a vector of BedRegion objects to a map of chromosome to intervals
 * @param bed_regions Optional vector of BedRegion objects to convert
 * @return Map of chromosome name to vector of sorted and merged intervals
 */
StrMap<vec<Interval>> ToIntervals(const std::optional<vec<BedRegion>>& bed_regions);

}  // namespace xoos::svc
