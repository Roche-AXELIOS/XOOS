#pragma once

#include <concepts>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

#include <htslib/hts.h>
#include <htslib/kseq.h>

#include <xoos/error/error.h>
#include <xoos/io/bed-region.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/io/htslib-util/kstring.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

namespace xoos::io {

/**
 * Parse one BED line without csv-parser's chunk buffering. The parser accepts whitespace-delimited BED fields and
 * returns no value for browser, track, and comment directives. Tab-delimited extra columns are preserved.
 * @param line The BED line to parse.
 * @param expect_extra_columns Whether to expect extra columns beyond the standard BED fields (chromosome, start, end).
 * @return An optional `BedRegion` if the line contains a valid BED region; otherwise, `std::nullopt`.
 */
std::optional<BedRegion> ParseBedRowFast(std::string_view line, bool expect_extra_columns);

template <typename Fn>
concept BedRegionCallback = std::invocable<Fn&, BedRegion&>;

/**
 * Parse a BED file and invoke `callback` for each valid region in file order.
 * Parsing stops early if the callback returns false.
 * @param bed_file_name The path to the BED file to parse.
 * @param callback The callback to invoke for each valid BED region.
 * The callback should accept a `BedRegion&` and may return `bool` to indicate whether parsing should continue.
 * @param expect_extra_columns Whether to expect extra columns beyond the standard BED fields (chromosome, start, end)
 */
template <BedRegionCallback Callback>
// We use forwarding references for `callback` to allow the use of
// temporary lambda expressions as well as other callable objects.
// std::forward is not used because the callable passed as an argument is
// invoked as lvalue repeatedly within the loop.
// A forwarded/moved callable cannot be safely invoked more than once.
// NOLINTBEGIN(cppcoreguidelines-missing-std-forward)
void ForEachBedRegionInFileFast(const fs::path& bed_file_name,  // NOSONAR cpp:S5536 - called by user code
                                Callback&& callback,            // NOSONAR cpp:S5425 - not forwarded, documented above
                                const bool expect_extra_columns) {
  const HtsFilePtr fp = HtsOpen(bed_file_name, "r");
  auto fetch_bed_line = [&fp, &bed_file_name](Kstring& line) {
    const s32 result = hts_getline(fp.get(), KS_SEP_LINE, line.Get());
    if (result == -1) {
      return false;
    }
    if (result < -1) {
      throw error::Error("Could not decompress BED file '{}'", bed_file_name);
    }
    return true;
  };
  Kstring line;
  bool saw_region = false;
  bool saw_empty_row = false;
  while (fetch_bed_line(line)) {
    auto bed_line = std::string_view(line.Get()->s, line.Get()->l);
    if (bed_line.find_first_not_of(" \t\r\n\v\f") == std::string_view::npos) {
      saw_empty_row = true;
      continue;
    }
    std::optional<BedRegion> parsed_bed_region = ParseBedRowFast(bed_line, expect_extra_columns);
    if (parsed_bed_region.has_value()) {
      saw_region = true;
      if constexpr (std::same_as<std::invoke_result_t<Callback&, BedRegion&>, bool>) {
        if (!callback(*parsed_bed_region)) {
          break;
        }
      } else {
        callback(*parsed_bed_region);
      }
    }
  }
  if (saw_empty_row && !saw_region) {
    throw error::Error("BED file is empty or does not contain any valid regions");
  }
}

// NOLINTEND(cppcoreguidelines-missing-std-forward)

/**
 * Parse a BED file and invoke `callback` for each valid region in file order.
 * Parsing stops early if the callback returns false.
 * @param bed_file_name The path to the BED file to parse.
 * @param callback The callback to invoke for each valid BED region.
 * The callback should accept a `BedRegion&` and may return `bool` to indicate whether parsing should continue.
 */
template <BedRegionCallback Callback>
void ForEachBedRegionInFileFast(const fs::path& bed_file_name,  // NOSONAR cpp:S5536 - called by user code
                                Callback&& callback) {
  ForEachBedRegionInFileFast(bed_file_name, std::forward<Callback>(callback), false);
}

}  // namespace xoos::io
