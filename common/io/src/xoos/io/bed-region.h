#pragma once

#include <concepts>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <csv.hpp>

#include <xoos/types/fs.h>
#include <xoos/types/int.h>

namespace xoos::io {

/**
 * A genomic region using 0-based, half-open coordinates [start, end).
 *
 * This is the canonical coordinate type for all internal XOOS code.
 * BED files natively use this convention, so no conversion is needed when reading BED.
 *
 * Conversion rules for external APIs:
 *   - sam_itr_queryi():    pass (start, end) directly (also 0-based half-open)
 *   - sam_itr_querys():    format a xoos::genome::Region, which renders as 1-based closed
 *   - faidx_fetch_seq():   pass (start, end - 1) to convert to 0-based closed
 *   - VCF output:          add 1 to position for 1-based VCF POS column
 *
 * To parse or render samtools-style region strings, use xoos::genome::Region and
 * xoos::genome::ParseSamtoolsLikeRegion in <xoos/genome/region.h>.
 */
struct BedRegion {
  // Chromosome name
  std::string chromosome;
  // Start position (0-based, inclusive)
  s64 start;
  // End position (0-based, exclusive)
  s64 end;
  // Extra data from the fourth column onward in the BED file
  std::vector<std::string> data;

  BedRegion(std::string in_chromosome, const s64 in_start, const s64 in_end)
      : chromosome(std::move(in_chromosome)), start(in_start), end(in_end) {
  }

  BedRegion(std::string in_chromosome, const s64 in_start, const s64 in_end, std::vector<std::string> in_data)
      : chromosome(std::move(in_chromosome)), start(in_start), end(in_end), data(std::move(in_data)) {
  }

  // Check if a 0-based position falls within this half-open interval
  bool Overlap(const s64 pos) const {
    return pos >= start && pos < end;
  }

  // Length of the interval in bases
  s64 Length() const {
    return end - start;
  }

  // Default comparison operator for lexicographical comparison
  // Equivalent to using `std::tie(chromosome, start, end)`
  // operator== is synthesized by the compiler from operator<=>
  auto operator<=>(const BedRegion&) const = default;
};

// Half-open interval overlap: two regions [s1, e1) and [s2, e2) overlap iff s1 < e2 && s2 < e1
inline bool RegionsOverlap(const BedRegion& r1, const BedRegion& r2) {
  return (r1.chromosome == r2.chromosome) && (r1.start < r2.end) && (r2.start < r1.end);
}

// Expand a 0-based half-open region by `padding` bases on each side, clamping start at 0.
inline BedRegion PadRegion(const BedRegion& region, const s64 padding) {
  return {region.chromosome, (region.start >= padding) ? region.start - padding : 0L, region.end + padding};
}

std::vector<BedRegion> PadBedRegions(std::vector<BedRegion>& regions, s64 padding);

// Parse a BED region string into a BedRegion.
// Accepts tab-delimited "chr\tstart\tend" (BED format, 0-based half-open) or
// "chr:start-end" where start and end are treated as 0-based half-open.
// For samtools-style 1-based closed strings, use xoos::genome::ParseSamtoolsLikeRegion instead.
BedRegion ParseRegion(const std::string& region);

// Whether `path`'s file name ends in the extension that marks a compressed BED file.
bool IsGzipPath(const fs::path& path);

// Decompress a gzip- or BGZF-compressed BED file into an in-memory stream for parsing.
std::stringstream DecompressToStream(const fs::path& path);

// A BED file is a specialized TSV, so a standard CSV parser configured for tabs can read it.
csv::CSVFormat MakeBedFormat();

// Validate a BED row and convert it to a BedRegion, or nullopt for browser/track/comment lines.
// `expect_extra_columns` raises the minimum column count from 3 to 4.
std::optional<BedRegion> ParseBedRow(const csv::CSVRow& row, bool expect_extra_columns);

// Parse a BED file, ignoring any columns beyond the third.
// A `.gz` file name extension selects gzip/BGZF decompression.
std::vector<BedRegion> ParseBedFile(const std::string& bed_file_name);

// use this if you need to parse additional columns from the bed file (fourth column onward), these will be stored in
// the `data` field of BedRegion
// will fail if the bed file has less than 4 columns
// A `.gz` file name extension selects gzip/BGZF decompression.
std::vector<BedRegion> ParseBedFileWithExtraColumns(const std::string& bed_file_name);

std::vector<BedRegion> PartitionRegions(const std::vector<BedRegion>& input_regions, s64 region_size);

std::vector<BedRegion> MergeOverlappingRegions(std::vector<BedRegion>&& input_regions);

// Invoke `callback` for every region yielded by `reader`
template <std::invocable<BedRegion> Callback>
void ForEachBedRegionInReader(csv::CSVReader& reader, bool expect_extra_columns, const Callback& callback) {
  bool saw_empty_row = false;
  bool saw_region = false;
  for (const csv::CSVRow& row : reader) {
    if (row.empty()) {
      saw_empty_row = true;
      continue;
    }
    std::optional<BedRegion> region = ParseBedRow(row, expect_extra_columns);
    if (!region.has_value()) {
      continue;
    }
    callback(std::move(*region));
    saw_region = true;
  }
  if (saw_empty_row && !saw_region) {
    throw std::runtime_error("BED file is empty or does not contain any valid regions");
  }
}

/**
 * Parse a BED file and invoke `callback` for each valid region, in file order.
 *
 * Rows are validated (numeric non-negative coordinates, start < end, at least 3 columns — or 4 when
 * `expect_extra_columns` is set); columns beyond the third are captured in `BedRegion::data`. UCSC
 * browser/track/comment (`#`) lines and blank lines are skipped, but a blank-only file is rejected
 * as malformed. A `.gz` file name selects gzip/BGZF decompression.
 *
 * @param bed_file_name Path to the BED file.
 * @param callback Invoked once per valid region.
 * @param expect_extra_columns Require a fourth column, failing rows that lack one.
 * @throws std::runtime_error if the file is malformed or cannot be decompressed.
 */
template <std::invocable<BedRegion> Callback>
void ForEachBedRegionInFile(const fs::path& bed_file_name, const Callback& callback, bool expect_extra_columns) {
  const csv::CSVFormat format = MakeBedFormat();
  if (IsGzipPath(bed_file_name)) {
    std::stringstream stream = DecompressToStream(bed_file_name);
    if (stream.tellp() == std::stringstream::pos_type(0)) {
      return;
    }
    csv::CSVReader reader(stream, format);
    ForEachBedRegionInReader(reader, expect_extra_columns, callback);
    return;
  }
  if (fs::file_size(bed_file_name) == 0) {
    return;
  }
  csv::CSVReader reader(bed_file_name.string(), format);
  ForEachBedRegionInReader(reader, expect_extra_columns, callback);
}

// Parse a BED file, invoking `callback` for each valid region and ignoring any columns beyond the third.
template <std::invocable<BedRegion> Callback>
void ForEachBedRegionInFile(const fs::path& bed_file_name, const Callback& callback) {
  ForEachBedRegionInFile(bed_file_name, callback, false);
}

}  // namespace xoos::io
