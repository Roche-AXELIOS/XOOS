#include "xoos/io/bed-region.h"

#include <filesystem>
#include <regex>

#include <csv.hpp>

#include <xoos/types/fs.h>

namespace xoos::io {

std::vector<BedRegion> PadBedRegions(std::vector<BedRegion>& regions, s64 padding) {
  if (padding == 0 || regions.empty()) {
    return regions;
  }

  std::vector<BedRegion> padded_regions;
  auto out_region = PadRegion(regions[0], padding);
  for (size_t i = 1; i < regions.size(); ++i) {
    if ((regions[i].start < regions[i - 1].start) && (regions[i].chromosome == regions[i - 1].chromosome)) {
      // NOSONAR
      throw std::runtime_error("PadBedRegions: Input bed regions are not sorted.");
    }
    auto next_region = PadRegion(regions[i], padding);
    if (RegionsOverlap(out_region, next_region)) {
      out_region.end = std::max(out_region.end, next_region.end);
    } else {
      padded_regions.push_back(out_region);
      out_region = next_region;
    }
  }
  padded_regions.push_back(out_region);
  return padded_regions;
}

/**
 * Parse a BED region string into a tuple of chrom, start, end
 * Allowable formats:
 * chr1:100-200
 * chr1\t100\t200
 **/
BedRegion ParseRegion(const std::string& region) {
  static const std::regex kRegionExpression(
      R"(^([0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*)[\t,:]{1}([0-9]+)[\t,-]{1}([0-9]+)[\t]*.*)");
  std::smatch match;
  if (std::regex_match(region, match, kRegionExpression)) {
    return {match[1].str(), static_cast<s64>(std::stoul(match[2].str())), static_cast<s64>(std::stoul(match[3].str()))};
  }
  throw std::runtime_error(std::string("Unable to parse region: " + region));
}

std::vector<BedRegion> ParseBedFile(const std::string& bed_file_name) {
  // since a BED file is a specialized TSV we can configure
  // a standard CSV parser to read a BED file.
  const csv::CSVFormat format = csv::CSVFormat().delimiter('\t').no_header().trim({' '});
  // if the file is empty, the reader will throw an error so we have to return an empty bed file first
  if (fs::file_size(bed_file_name) == 0) {
    return {};
  }

  csv::CSVReader reader(bed_file_name, format);

  std::vector<BedRegion> regions;
  for (const auto& row : reader) {
    auto chrom = row[0].get<std::string>();
    // skip two special cases which are used to configure UCSC genome browser
    if (chrom.starts_with("browser") || chrom.starts_with("track") || chrom.starts_with("#")) {
      continue;
    }
    if (row.size() < 3) {
      throw std::runtime_error("BED file must have at least 3 columns");
    }
    if (row[1].is_str() || row[2].is_str()) {
      throw std::runtime_error("BED file must have numeric start and end columns");
    }
    if (row[1].get<s64>() < 0 || row[2].get<s64>() < 0) {
      throw std::runtime_error("BED file must have non-negative start and end columns");
    }
    if (row[1].get<s64>() >= row[2].get<s64>()) {
      throw std::runtime_error("BED file must have start < end");
    }
    regions.emplace_back(chrom, row[1].get<s64>(), row[2].get<s64>());
  }
  return regions;
}

/**
 * Partition the regions into smaller regions of size region_size, this is done to
 * achieve better parallelism when computing the coverage histogram. The regions are fixed size
 * except for the last region which may be limited by reaching the end of the original region.
 * @param input_regions
 * @param region_size
 * @return
 */
std::vector<BedRegion> PartitionRegions(const std::vector<BedRegion>& input_regions, s64 region_size) {
  auto partitioned_regions = std::vector<BedRegion>();
  for (const auto& region : input_regions) {
    // Determine the full length of the target region and the number of regions it will be partitioned into
    const s64 target_len = region.end - region.start + 1;
    const s64 target_region_count = (target_len + region_size - 1) / region_size;

    // Generate a region for each partitioned region, the last region may be smaller than region_size
    // if the target region is not evenly divisible by region_size
    for (auto j = 0; j < target_region_count; ++j) {
      const s64 beg = region.start + j * region_size;
      const s64 end = std::min(region.end, region.start + (j + 1) * region_size);
      partitioned_regions.emplace_back(region.chromosome, beg, end);
    }
  }
  return partitioned_regions;
}

/**
 * Merge overlapping (not including adjacent) bed regions.
 *
 * If the input regions are not sorted, the input regions will be sorted.
 *
 * @param input_regions Sorted list of bed regions.
 */
std::vector<BedRegion> MergeOverlappingRegions(std::vector<BedRegion>&& input_regions) {
  if (!std::ranges::is_sorted(input_regions)) {
    std::ranges::sort(input_regions);
  }
  std::vector<BedRegion> merged_regions;
  if (!input_regions.empty()) {
    merged_regions.push_back(input_regions[0]);
    for (size_t i = 1; i < input_regions.size(); ++i) {
      const auto& curr = input_regions[i];
      auto& prev = merged_regions.back();
      if (prev.chromosome == curr.chromosome && curr.start < prev.end) {
        prev.end = std::max(prev.end, curr.end);
      } else {
        merged_regions.push_back(curr);
      }
    }
  }
  return merged_regions;
}

}  // namespace xoos::io
