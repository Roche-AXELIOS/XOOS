#pragma once

#include <string>
#include <utility>
#include <vector>

#include <xoos/types/int.h>

namespace xoos::io {

/**
 *  A region read as directly from a bed file.
 */
struct BedRegion {
  // Chromosome name
  std::string chromosome;
  /**
   * Start position (0-based, inclusive)
   *
   * NOTE: Some modules may be incorrectly using 1-based start positions with this struct
   */
  s64 start;
  /**
   * End position (0-based, exclusive)
   *
   * NOTE: Some modules may be incorrectly using 1-based end positions with this struct
   */
  s64 end;

  BedRegion(std::string in_chromosome, const s64 in_start, const s64 in_end)
      : chromosome(std::move(in_chromosome)), start(in_start), end(in_end) {
  }

  bool Overlap(const s64 pos) const {
    return pos >= start && pos < end;
  }

  // Default comparison operator for lexicographical comparison
  // Equivalent to using `std::tie(chromosome, start, end)`
  auto operator<=>(const BedRegion&) const = default;
  bool operator==(const BedRegion&) const = default;
};

inline bool RegionsOverlap(const BedRegion& r1, const BedRegion& r2) {
  return (r1.chromosome == r2.chromosome) && (r1.start <= r2.end) && (r1.end >= r2.start);
}

inline BedRegion PadRegion(const BedRegion& region, const s64 padding) {
  return {region.chromosome, (region.start >= padding) ? region.start - padding : 0l, region.end + padding};
}

std::vector<BedRegion> PadBedRegions(std::vector<BedRegion>& regions, s64 padding);

// Parse a BED region string into a tuple of chrom, start, end
BedRegion ParseRegion(const std::string& region);

std::vector<BedRegion> ParseBedFile(const std::string& bed_file_name);

std::vector<BedRegion> PartitionRegions(const std::vector<BedRegion>& input_regions, s64 region_size);

std::vector<BedRegion> MergeOverlappingRegions(std::vector<BedRegion>&& input_regions);

}  // namespace xoos::io
