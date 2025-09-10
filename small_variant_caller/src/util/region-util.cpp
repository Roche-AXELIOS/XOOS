#include "region-util.h"

#include <algorithm>
#include <fstream>

#include <csv.hpp>

#include <xoos/types/vec.h>
#include <xoos/util/string-functions.h>

#include "xoos/error/error.h"

namespace xoos::svc {

/**
 * @brief Return whether two given regions overlap.
 * @param start1 Region 1 start (inclusive)
 * @param end1 Region 1 end (exclusive)
 * @param start2 Region 2 start (inclusive)
 * @param end2 Region 2 end (exclusive)
 * @return Regions overlap or not
 */
bool IntervalOverlap(const u64 start1, const u64 end1, const u64 start2, const u64 end2) {
  return (start1 <= start2 && start2 < end1) || (start1 < end2 && end2 <= end1) ||
         (start2 <= start1 && start1 < end2) || (start2 < end1 && end1 <= end2);
}

/**
 * @brief Return whether Interval overlaps with a given region.
 * @param i Interval
 * @param start Region start (inclusive)
 * @param end Region end (exclusive)
 * @return Region overlap or not
 */
bool IntervalOverlap(const Interval& i, const u64 start, const u64 end) {
  return IntervalOverlap(i.start, i.end, start, end);
}

/**
 * @brief Return whether two Interval overlap.
 * @param a Interval a
 * @param b Interval b
 * @return Interval overlap or not
 */
bool IntervalOverlap(const Interval& a, const Interval& b) {
  return IntervalOverlap(a.start, a.end, b.start, b.end);
}

/**
 * @brief Return whether intervals overlap with given region.
 * @param sorted_intervals Sorted vector of Interval
 * @param start Region start (inclusive)
 * @param end Region end (exclusive)
 * @return Region overlap or not
 */
bool IntervalOverlap(const vec<Interval>& sorted_intervals, const u64 start, const u64 end) {
  for (const auto& i : sorted_intervals) {
    if (i.start > end) {
      return false;
    }
    if (IntervalOverlap(i, start, end)) {
      return true;
    }
  }
  return false;
}

/**
 * @brief Merge intervals if they overlap or are within callapsible distance of each other.
 * @param sorted_intervals Sorted vector of Interval
 * @param collapse_dist Distance to collapse nearby intervals
 * @return Sorted vector of non-overlapping Interval
 */
vec<Interval> MergeIntervals(vec<Interval>& sorted_intervals, const u64 collapse_dist) {
  vec<Interval> results;
  if (!sorted_intervals.empty()) {
    results.emplace_back(sorted_intervals[0]);
    for (size_t i = 1; i < sorted_intervals.size(); ++i) {
      const Interval& curr = sorted_intervals[i];
      Interval& prev = results.back();
      if (IntervalOverlap(prev, curr) || prev.end + collapse_dist >= curr.start) {
        prev.end = std::max(prev.end, curr.end);
        // we do not need to update `prev.start` because intervals are sorted
      } else {
        results.emplace_back(curr);
      }
    }
  }
  return results;
}

/**
 * @brief Parse BED file to create a map of chromosome to vector of Interval.
 * @param region_bed Path to BED file
 * @return Map of chromosome to vector of Interval
 */
std::optional<ChromIntervalsMap> GetChromIntervalMap(const std::optional<fs::path>& region_bed) {
  if (!region_bed) {
    return std::nullopt;
  }

  // parse regions from BED file
  std::ifstream bed_file_stream(region_bed.value());
  std::string line;
  ChromIntervalsMap results;
  while (std::getline(bed_file_stream, line)) {
    if (line.starts_with("browser") || line.starts_with("track") || line.starts_with("#")) {
      continue;
    }
    string::Trim(line);
    if (line.empty()) {
      continue;
    }
    const BedRegion r = io::ParseRegion(line);
    const u64 start = r.start >= 0 ? static_cast<u64>(r.start) : 0;
    const u64 end = r.end >= 0 ? static_cast<u64>(r.end) : 0;
    results[r.chromosome].emplace_back(Interval{start, end});
  }
  if (results.empty()) {
    throw error::Error("BED file contains no regions");
  }
  // sort and merge overlapping intervals for each chromosome
  for (auto& [chrom, intervals] : results) {
    std::sort(intervals.begin(), intervals.end());
    results[chrom] = MergeIntervals(intervals);
  }
  return results;
}

/**
 * @brief Write a BED file from a map of chromosome to vector of intervals.
 * @param chrom_map Map of chromosome to vector of intervals
 * @param out_bed_path Output BED file path
 */
void WriteBed(const ChromIntervalsMap& chrom_map, const fs::path& out_bed_path) {
  std::ofstream ofs(out_bed_path, std::ios::out);
  auto writer = csv::make_tsv_writer(ofs);
  for (const auto& [chrom, intervals] : chrom_map) {
    for (const auto& i : intervals) {
      writer << vec<std::string>({chrom, std::to_string(i.start), std::to_string(i.end)});
    }
  }
}

/**
 * @brief Extract a vector of BedRegion fom a region string or a path to BED file.
 * @param region_bed Path to BED file
 * @return Vector of BedRegion
 */
std::optional<vec<BedRegion>> GetBedRegions(const std::optional<fs::path>& region_bed) {
  if (!region_bed) {
    return std::nullopt;
  }

  return io::ParseBedFile(region_bed.value());
}

StrMap<vec<Interval>> ToIntervals(const std::optional<vec<BedRegion>>& bed_regions) {
  if (!bed_regions) {
    return {};
  }

  // Use a vector of Intervals per chromosome vs BedRegion structs for storing BED regions. Makes lookups easier for
  // getting which intervals overlap a given position on a chromosome.
  StrMap<vec<Interval>> chrom_intervals;

  // Convert the vector of structs into a map so we can easily query
  for (const auto& region : *bed_regions) {
    chrom_intervals[region.chromosome].emplace_back(region.start, region.end);
  }

  for (auto& [chrom, intervals] : chrom_intervals) {
    std::ranges::sort(intervals);
    chrom_intervals[chrom] = MergeIntervals(intervals);
  }

  return chrom_intervals;
}

}  // namespace xoos::svc
