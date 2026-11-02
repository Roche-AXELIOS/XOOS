#include "util/region-util.h"

#include <algorithm>
#include <array>
#include <limits>
#include <optional>
#include <string>
#include <utility>

#include <csv.hpp>

#include <xoos/error/error.h>
#include <xoos/io/fasta-reader.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/log/logging.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "core/interval.h"
#include "io/bed-reader.h"

namespace xoos::alignment_metrics {

RegionsByChromosome MergeOverlappingRegions(RegionsByChromosome regions) {
  RegionsByChromosome merged_regions;
  for (auto& [chromosome, intervals] : regions) {
    if (!intervals.empty()) {
      if (!std::ranges::is_sorted(intervals)) {
        std::ranges::sort(intervals);
      }
      merged_regions[chromosome].emplace_back(intervals.at(0));
      for (size_t i = 1; i < intervals.size(); ++i) {
        const auto& curr = intervals.at(i);
        auto& prev = merged_regions[chromosome].back();
        if (curr.start < prev.end) {
          prev.end = std::max(prev.end, curr.end);
        } else {
          merged_regions[chromosome].emplace_back(curr);
        }
      }
    }
  }
  return merged_regions;
}

namespace {

// Extend a partition end position rightward so it does not fall in the middle of a homopolymer run, returning the
// adjusted end. Reads reference bases from `fasta_reader` on `chromosome` starting at `end`.
s64 ExtendEndToHomopolymerBoundary(io::FastaReader& fasta_reader, const std::string& chromosome, const s64 end) {
  s32 buffer_length = 1;
  if (end + buffer_length > std::numeric_limits<s32>::max()) {
    throw error::Error("End position is greater than max int");
  }
  auto ref_seq = fasta_reader.GetSequence(chromosome, static_cast<s32>(end - 1), static_cast<s32>(end + buffer_length));
  if (ref_seq.length() >= 2) {
    const char last_base = ref_seq[0];
    while (last_base == ref_seq[ref_seq.length() - 1] && last_base != 'N') {
      if (end + buffer_length + 1 > std::numeric_limits<s32>::max()) {
        throw error::Error("End position is greater than max int");
      }
      const std::string extension = fasta_reader.GetSequence(
          chromosome, static_cast<s32>(end + buffer_length), static_cast<s32>(end + buffer_length + 1));
      if (extension.empty()) {
        break;
      }
      ++buffer_length;
      ref_seq += extension;
    }
  }
  return end + buffer_length - 1;
}

// Partition a single interval into region_size-sized regions, appending them to `out`. When `hp_aware` and a
// reader is available, non-final region ends are extended to the end of any homopolymer they land in.
void PartitionInterval(const Interval& interval,
                       const std::string& chromosome,
                       const u32 region_size,
                       const bool hp_aware,
                       std::optional<io::FastaReader>& fasta_reader,
                       vec<Interval>& out) {
  // Determine the full length of the target region and the number of regions it will be partitioned into.
  // Interval is 0-based half-open [start, end), so length = end - start.
  const s64 target_len = interval.end - interval.start;
  const s64 target_region_count = (target_len + region_size - 1) / region_size;

  // Generate a region for each partitioned region, the last region may be smaller than region_size
  // if the target region is not evenly divisible by region_size
  s64 beg = interval.start;
  s64 end = interval.start + region_size;
  for (auto j = 0; j < target_region_count; ++j) {
    // We have reached the end of the region
    if (beg >= interval.end) {
      break;
    }
    // Check if the end position is in the middle of a homopolymer
    if (j < target_region_count - 1 && hp_aware && fasta_reader.has_value()) {
      // If not the last region, then extend the end position to the end of the homopolymer
      end = ExtendEndToHomopolymerBoundary(fasta_reader.value(), chromosome, end);
    }
    end = std::min(end, interval.end);
    // Every partitioned chunk inherits the label bit vector of the interval it was cut from.
    out.emplace_back(beg, end, interval.labels);

    beg = end;
    end = std::min(beg + region_size, interval.end);
  }
}

}  // namespace

RegionsByChromosome PartitionRegions(const RegionsByChromosome& regions,
                                     const u32 region_size,
                                     const std::optional<fs::path>& reference_path,
                                     const bool hp_aware) {
  RegionsByChromosome partitioned_regions;
  std::optional<io::FastaReader> fasta_reader;
  if (reference_path.has_value() && exists(reference_path.value())) {
    fasta_reader.emplace(reference_path.value());
  }
  for (const auto& [chromosome, intervals] : regions) {
    for (const auto interval : intervals) {
      PartitionInterval(interval, chromosome, region_size, hp_aware, fasta_reader, partitioned_regions[chromosome]);
    }
  }
  return partitioned_regions;
}

}  // namespace xoos::alignment_metrics
