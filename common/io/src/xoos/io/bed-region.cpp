#include "xoos/io/bed-region.h"

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <regex>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <htslib/bgzf.h>

#include <csv.hpp>

#include <xoos/error/error.h>
#include <xoos/types/int.h>

#include "xoos/io/htslib-util/htslib-ptr.h"

namespace xoos::io {

namespace fs = std::filesystem;

namespace {

// Extension marking a compressed BED file. Compression is detected by extension alone, so a
// compressed file that is not named `.gz` is still read as plain text.
constexpr std::string_view kGzipExtension{".gz"};

// Number of decompressed bytes requested per read from a compressed BED file.
constexpr size_t kDecompressChunkSize = 64UL * 1024UL;

// Column positions defined by the BED format. Columns from `kFirstExtraColumn` onward are optional
// and are captured verbatim in `BedRegion::data`.
constexpr size_t kChromosomeColumn = 0;
constexpr size_t kStartColumn = 1;
constexpr size_t kEndColumn = 2;
constexpr size_t kFirstExtraColumn = 3;

// Minimum column counts and the errors reported when a row has fewer, for each parsing mode.
constexpr size_t kMinRegionColumns = 3;
constexpr size_t kMinExtraColumns = 4;

/// Regex match group indices for ParseRegion: group 0 is the full match.
constexpr size_t kRegexChromGroup = 1;
constexpr size_t kRegexStartGroup = 2;
constexpr size_t kRegexEndGroup = 3;

}  // namespace

bool IsGzipPath(const fs::path& path) {
  return path.extension() == kGzipExtension;
}

/**
 * Decompress a gzip- or BGZF-compressed BED file into an in-memory stream.
 *
 * The whole file is decompressed up front because the CSV parser seeks within its source, which a
 * decompression stream cannot do cheaply. Chunks are written straight into the returned stream so
 * that no second copy of the decompressed text is held. A BED file holds one interval per line, so
 * the decompressed text is comparable in size to the interval set the caller builds from it.
 */
std::stringstream DecompressToStream(const fs::path& path) {
  const BgzfPtr file{bgzf_open(path.c_str(), "r")};  // NOSONAR(S5950) — bgzf_open returns an owning raw pointer
  if (file == nullptr) {
    throw error::Error("Could not open compressed BED file '{}'", path);
  }
  std::stringstream content;
  std::string chunk(kDecompressChunkSize, '\0');
  while (true) {
    // bgzf_read returns a signed byte count, negative on failure.
    const s64 bytes_read =
        bgzf_read(file.get(), chunk.data(), chunk.size());  // NOSONAR(S5356) — htslib C API requires void*
    if (bytes_read < 0) {
      throw error::Error("Could not decompress BED file '{}'", path);
    }
    if (bytes_read == 0) {
      break;
    }
    content.write(chunk.data(), bytes_read);
  }
  return content;
}

// A BED file is a specialized TSV, so a standard CSV parser configured for tabs can read it.
csv::CSVFormat MakeBedFormat() {
  return csv::CSVFormat().delimiter('\t').no_header().trim({' '});
}

/**
 * Validate a non-empty BED row and convert it to a BedRegion. Returns nullopt for lines that carry
 * no region, namely the UCSC genome browser configuration lines and comments. `expect_extra_columns`
 * raises the minimum column count from 3 to 4 so callers that need the fourth column onward fail
 * fast on rows that lack it.
 */
std::optional<BedRegion> ParseBedRow(const csv::CSVRow& row, bool expect_extra_columns) {
  std::string chrom = row[kChromosomeColumn].get<std::string>();
  // skip two special cases which are used to configure UCSC genome browser
  if (chrom.starts_with("browser") || chrom.starts_with("track") || chrom.starts_with("#")) {
    return std::nullopt;
  }
  if (row.size() < (expect_extra_columns ? kMinExtraColumns : kMinRegionColumns)) {
    throw error::Error(
        expect_extra_columns
            ? "BED file must have at least 4 columns in order to additionally parse data from the fourth column onward"
            : "BED file must have at least 3 columns");
  }
  // Reject non-numeric fields, including empty ones: an empty column parses as CSV_NULL (not
  // CSV_STRING), so guarding on is_str() alone would fall through to get<s64>() and surface the
  // csv parser's raw "Not a number." message instead of this BED-specific one.
  if (!row[kStartColumn].is_num() || !row[kEndColumn].is_num()) {
    throw error::Error("BED file must have numeric start and end columns");
  }
  const s64 start = row[kStartColumn].get<s64>();
  const s64 end = row[kEndColumn].get<s64>();
  if (start < 0 || end < 0) {
    throw error::Error("BED file must have non-negative start and end columns");
  }
  if (start >= end) {
    throw error::Error("BED file must have start < end");
  }
  std::vector<std::string> data;
  for (size_t i = kFirstExtraColumn; i < row.size(); ++i) {
    data.push_back(row[i].get<std::string>());
  }
  return BedRegion{std::move(chrom), start, end, std::move(data)};
}

std::vector<BedRegion> PadBedRegions(std::vector<BedRegion>& regions, const s64 padding) {
  if (padding == 0 || regions.empty()) {
    return regions;
  }

  std::vector<BedRegion> padded_regions;
  auto out_region = PadRegion(regions[0], padding);
  for (size_t i = 1; i < regions.size(); ++i) {
    if ((regions[i].start < regions[i - 1].start) && (regions[i].chromosome == regions[i - 1].chromosome)) {
      // NOSONAR
      throw error::Error("PadBedRegions: Input bed regions are not sorted.");
    }
    const auto next_region = PadRegion(regions[i], padding);
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
    return {match[kRegexChromGroup].str(),
            static_cast<s64>(std::stoul(match[kRegexStartGroup].str())),
            static_cast<s64>(std::stoul(match[kRegexEndGroup].str()))};
  }
  throw error::Error("Unable to parse region: {}", region);
}

std::vector<BedRegion> ParseBedFile(const std::string& bed_file_name) {
  std::vector<BedRegion> regions;
  ForEachBedRegionInFile(bed_file_name, [&regions](BedRegion region) {
    // Drop any extra columns to preserve the original 3-column BedRegion shape.
    regions.emplace_back(std::move(region.chromosome), region.start, region.end);
  });
  return regions;
}

/**
 * Same logic as ParseBedFile, except also record additional columns from the fourth column onward in the bed file and
 * store them in the `data` field of BedRegion.
 */
std::vector<BedRegion> ParseBedFileWithExtraColumns(const std::string& bed_file_name) {
  std::vector<BedRegion> regions;
  ForEachBedRegionInFile(bed_file_name, [&regions](BedRegion region) { regions.push_back(std::move(region)); }, true);
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
std::vector<BedRegion> PartitionRegions(const std::vector<BedRegion>& input_regions, const s64 region_size) {
  auto partitioned_regions = std::vector<BedRegion>();
  for (const auto& region : input_regions) {
    // Determine the full length of the target region and the number of regions it will be partitioned into
    const s64 target_len = region.end - region.start;
    const s64 target_region_count = (target_len + region_size - 1) / region_size;

    // Generate a region for each partitioned region, the last region may be smaller than region_size
    // if the target region is not evenly divisible by region_size
    for (auto j = 0; j < target_region_count; ++j) {
      const s64 beg = region.start + (j * region_size);
      const s64 end = std::min(region.end, region.start + ((j + 1) * region_size));
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
  // Claim ownership of the caller's vector and use it as the working buffer.
  std::vector<BedRegion> regions = std::move(input_regions);
  if (!std::ranges::is_sorted(regions)) {
    std::ranges::sort(regions);
  }
  std::vector<BedRegion> merged_regions;
  if (!regions.empty()) {
    merged_regions.push_back(regions[0]);
    for (size_t i = 1; i < regions.size(); ++i) {
      const auto& curr = regions[i];
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
