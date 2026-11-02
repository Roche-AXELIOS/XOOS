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

/// BED region representation provided by the shared I/O library.
using BedRegion = io::BedRegion;

/**
 * @brief Interval on a chromosome using 0-based, half-open coordinates `[start, end)`.
 */
struct Interval {
  /// 0-based inclusive start position.
  u64 start{};
  /// 0-based exclusive end position.
  u64 end{};
  /// Compare intervals lexicographically by start position and then end position.
  auto operator<=>(const Interval&) const = default;
};

/**
 * @brief Partitioned genomic region used for parallel feature computation or variant filtering.
 */
struct TargetRegion {
  /// Chromosome name.
  std::string chrom;
  /// 0-based inclusive start position on the reference.
  u64 start{};
  /// 0-based exclusive end position on the reference.
  u64 end{};
};

/// Chromosome-keyed collection of sorted genomic intervals.
using ChromIntervalsMap = StrMap<vec<Interval>>;

/**
 * @brief Return whether two 0-based, half-open regions overlap.
 * @param start1 First region start, inclusive.
 * @param end1 First region end, exclusive.
 * @param start2 Second region start, inclusive.
 * @param end2 Second region end, exclusive.
 * @return `true` when the regions overlap.
 */
bool IntervalOverlap(u64 start1, u64 end1, u64 start2, u64 end2);

/**
 * @brief Return whether an interval overlaps a 0-based, half-open region.
 * @param i Interval to test.
 * @param start Region start, inclusive.
 * @param end Region end, exclusive.
 * @return `true` when the interval and region overlap.
 */
bool IntervalOverlap(const Interval& i, u64 start, u64 end);

/**
 * @brief Return whether two intervals overlap.
 * @param a First interval.
 * @param b Second interval.
 * @return `true` when the intervals overlap.
 */
bool IntervalOverlap(const Interval& a, const Interval& b);

/**
 * @brief Return whether any sorted interval overlaps a 0-based, half-open region.
 * @param sorted_intervals Intervals sorted by start position.
 * @param start Region start, inclusive.
 * @param end Region end, exclusive.
 * @return `true` when at least one interval overlaps the region.
 */
bool IntervalOverlap(const vec<Interval>& sorted_intervals, u64 start, u64 end);

/**
 * @brief Return whether a 0-based position is in one of the chromosome's half-open intervals.
 * @param intervals_by_chromosome Sorted, non-overlapping intervals keyed by chromosome.
 * @param chromosome Chromosome containing the position.
 * @param position 0-based position to query.
 * @return `true` when an interval contains the position.
 */
bool IsPositionInIntervals(const ChromIntervalsMap& intervals_by_chromosome,
                           const std::string& chromosome,
                           u64 position);

/**
 * @brief Merge intervals that overlap or are within a specified distance.
 * @param sorted_intervals Intervals sorted by start position.
 * @param collapse_dist Maximum distance between intervals that should be merged.
 * @return Sorted, non-overlapping merged intervals.
 */
vec<Interval> MergeIntervals(vec<Interval>& sorted_intervals, u64 collapse_dist = 0);

/**
 * @brief Parse a BED file into sorted, merged intervals keyed by chromosome.
 * @param region_bed Optional path to the BED file.
 * @return Parsed intervals, or `std::nullopt` when no path is supplied.
 * @throws error::Error When the supplied BED file contains no regions.
 */
std::optional<ChromIntervalsMap> GetChromIntervalMap(const std::optional<fs::path>& region_bed);

/**
 * @brief Write chromosome intervals as BED records.
 * @param chrom_map Intervals keyed by chromosome.
 * @param out_bed_path Output BED path.
 */
void WriteBed(const ChromIntervalsMap& chrom_map, const fs::path& out_bed_path);

/**
 * @brief Parse an optional BED file into BED regions.
 * @param region_bed Optional path to the BED file.
 * @return Parsed BED regions, or `std::nullopt` when no path is supplied.
 */
std::optional<vec<BedRegion>> GetBedRegions(const std::optional<fs::path>& region_bed);

/**
 * @brief Convert BED regions into sorted, merged intervals keyed by chromosome.
 * @param bed_regions Optional BED regions to convert.
 * @return Converted intervals, or an empty map when no regions are supplied.
 */
StrMap<vec<Interval>> ToIntervals(const std::optional<vec<BedRegion>>& bed_regions);

}  // namespace xoos::svc
