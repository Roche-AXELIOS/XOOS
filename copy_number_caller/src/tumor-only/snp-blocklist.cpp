#include "tumor-only/snp-blocklist.h"

#include <algorithm>

#include <xoos/io/bed-region.h>
#include <xoos/log/logging.h>

namespace xoos::cnc {

SnpBlocklist::SnpBlocklist(const fs::path& bed_file) {
  const auto bed_regions = io::ParseBedFile(bed_file.string());
  for (const auto& region : bed_regions) {
    // Skip empty or malformed intervals defensively; a BED region should be half-open with start < end.
    if (region.end > region.start) {
      _intervals[region.chromosome].emplace_back(region.start, region.end);
    }
  }
  size_t total = 0;
  for (auto& [chrom, intervals] : _intervals) {
    // Sort by start, then coalesce overlapping/adjacent intervals so Contains() can binary-search.
    std::ranges::sort(intervals);
    vec<std::pair<s64, s64>> merged;
    for (const auto& interval : intervals) {
      if (!merged.empty() && interval.first <= merged.back().second) {
        merged.back().second = std::max(merged.back().second, interval.second);
      } else {
        merged.push_back(interval);
      }
    }
    intervals = std::move(merged);
    total += intervals.size();
  }
  Logging::Info(
      "Loaded SNP blocklist with {} regions across {} contigs from {}", total, _intervals.size(), bed_file.string());
}

bool SnpBlocklist::Contains(const std::string& chromosome, const s64 position) const {
  const auto chrom_it = _intervals.find(chromosome);
  if (chrom_it == _intervals.end()) {
    return false;
  }
  const auto& intervals = chrom_it->second;
  // Find the first interval whose start is strictly greater than position; the candidate
  // containing interval is the one immediately before it (largest start <= position).
  const auto it = std::ranges::upper_bound(intervals, position, std::ranges::less{}, &std::pair<s64, s64>::first);
  if (it == intervals.begin()) {
    return false;
  }
  const auto& [start, end] = *std::prev(it);
  // Half-open: start inclusive, end exclusive.
  return position >= start && position < end;
}

bool SnpBlocklist::Empty() const {
  return _intervals.empty();
}

size_t SnpBlocklist::Size() const {
  size_t total = 0;
  for (const auto& [chrom, intervals] : _intervals) {
    total += intervals.size();
  }
  return total;
}

}  // namespace xoos::cnc
