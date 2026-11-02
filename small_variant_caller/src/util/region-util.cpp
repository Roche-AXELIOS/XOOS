#include "region-util.h"

#include <algorithm>

#include "locked-tsv-writer.h"
#include "xoos/error/error.h"

namespace xoos::svc {

bool IntervalOverlap(const u64 start1, const u64 end1, const u64 start2, const u64 end2) {
  return (start1 <= start2 && start2 < end1) || (start1 < end2 && end2 <= end1) ||
         (start2 <= start1 && start1 < end2) || (start2 < end1 && end1 <= end2);
}

bool IntervalOverlap(const Interval& i, const u64 start, const u64 end) {
  return IntervalOverlap(i.start, i.end, start, end);
}

bool IntervalOverlap(const Interval& a, const Interval& b) {
  return IntervalOverlap(a.start, a.end, b.start, b.end);
}

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

bool IsPositionInIntervals(const ChromIntervalsMap& intervals_by_chromosome,
                           const std::string& chromosome,
                           const u64 position) {
  const auto chromosome_it = intervals_by_chromosome.find(chromosome);
  if (chromosome_it == intervals_by_chromosome.end()) {
    return false;
  }

  for (const auto& interval : chromosome_it->second) {
    if (position < interval.start) {
      return false;
    }
    if (position < interval.end) {
      return true;
    }
  }
  return false;
}

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

std::optional<ChromIntervalsMap> GetChromIntervalMap(const std::optional<fs::path>& region_bed) {
  if (!region_bed) {
    return std::nullopt;
  }

  const vec<BedRegion> bed_regions = io::ParseBedFile(region_bed->string());
  ChromIntervalsMap results;
  for (const BedRegion& region : bed_regions) {
    results[region.chromosome].emplace_back(
        Interval{.start = static_cast<u64>(region.start), .end = static_cast<u64>(region.end)});
  }
  if (results.empty()) {
    throw error::Error("BED file contains no regions");
  }
  // sort and merge overlapping intervals for each chromosome
  for (auto& [chrom, intervals] : results) {
    std::ranges::sort(intervals);
    results[chrom] = MergeIntervals(intervals);
  }
  return results;
}

void WriteBed(const ChromIntervalsMap& chrom_map, const fs::path& out_bed_path) {
  LockedTsvWriter writer(out_bed_path);
  for (const auto& [chrom, intervals] : chrom_map) {
    for (const auto& i : intervals) {
      writer.AppendRow({chrom, std::to_string(i.start), std::to_string(i.end)});
    }
  }
}

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
