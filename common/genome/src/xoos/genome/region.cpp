#include "xoos/genome/region.h"

#include <algorithm>
#include <iterator>
#include <limits>
#include <string>
#include <string_view>

#include <xoos/error/error.h>
#include <xoos/util/parse-int.h>

namespace xoos::genome {

namespace {

// A coordinate suffix looks like "123", "123-456", or with thousands separators
// "1,000-2,000". It must contain at least one digit and consist solely of
// digits, commas and at most one dash separating two coordinate parts.
// This mirrors how samtools decides whether the text after the last colon is a
// coordinate range or part of a sequence name (e.g. an HLA allele "HLA:A").
bool LooksLikeCoordinates(const std::string_view text) {
  if (text.empty()) {
    return false;
  }
  bool seen_digit = false;
  s32 dash_count = 0;
  for (const char c : text) {
    if (c >= '0' && c <= '9') {
      seen_digit = true;
    } else if (c == ',') {
      // thousands separator, ignored during parsing
    } else if (c == '-') {
      ++dash_count;
    } else {
      return false;
    }
  }
  return seen_digit && dash_count <= 1;
}

// Remove thousands-separator commas from a coordinate token.
std::string StripCommas(std::string_view token) {
  std::string result;
  result.reserve(token.size());
  std::ranges::copy_if(token, std::back_inserter(result), [](const char c) { return c != ','; });
  return result;
}

// Length of the contig-name portion of a region string. The contig is
// everything before the *last* colon, but only when the text after that colon
// looks like a coordinate range; otherwise the whole string is the contig name.
// Both ParseSamtoolsLikeRegion and ContigName rely on this so they always agree
// on where the name ends (important for colon-containing names like "HLA:A").
std::string_view::size_type ContigNameLength(const std::string_view region) {
  const auto last_colon = region.rfind(':');
  if (last_colon == std::string_view::npos) {
    return region.size();
  }
  if (!LooksLikeCoordinates(region.substr(last_colon + 1))) {
    return region.size();
  }
  return last_colon;
}

}  // namespace

Region ParseSamtoolsLikeRegion(const std::string& region) {
  return ParseSamtoolsLikeRegion(region, RegionForm::kAny);
}

Region ParseSamtoolsLikeRegion(const std::string& region, const RegionForm form) {
  if (region.empty()) {
    throw error::Error("region string is empty");
  }

  // Whole-chromosome span used when no coordinates are supplied.
  constexpr s64 kMaxEnd = std::numeric_limits<s64>::max();

  const bool allow_whole_chromosome = form == RegionForm::kAny;

  // Resolve the contig/coordinate boundary the same way ContigName does, so the
  // two functions always agree on where the name ends.
  const std::string_view region_view{region};
  const auto name_len = ContigNameLength(region_view);
  if (name_len == region_view.size()) {
    // No coordinate suffix -> whole chromosome (e.g. "chr" or "HLA:A").
    if (!allow_whole_chromosome) {
      throw error::Error("region '{}' must include coordinates (chr:start-end)", region);
    }
    return {.chromosome = region, .start = 0, .end = kMaxEnd};
  }

  const std::string_view suffix = region_view.substr(name_len + 1);
  const std::string chromosome = region.substr(0, name_len);
  if (chromosome.empty()) {
    throw error::Error("region '{}' has an empty chromosome name", region);
  }

  // Split the coordinate suffix into start and (optional) end on the single dash.
  const auto dash = suffix.find('-');
  const auto start_token = StripCommas(dash == std::string_view::npos ? suffix : suffix.substr(0, dash));
  if (start_token.empty()) {
    throw error::Error("region '{}' has an empty start coordinate", region);
  }

  // Largest coordinate that still fits in s64 once converted. ParseU64 yields a
  // u64, so a value above INT64_MAX would wrap to a negative s64 and produce an
  // invalid Region; reject it explicitly.
  constexpr auto kMaxCoord = static_cast<u64>(std::numeric_limits<s64>::max());

  const u64 start_1based = util::ParseU64(start_token);
  if (start_1based == 0) {
    throw error::Error("region must be 1-based, but {} has start of 0", region);
  }
  if (start_1based > kMaxCoord) {
    throw error::Error("region '{}' start coordinate ({}) exceeds the maximum of {}", region, start_1based, kMaxCoord);
  }
  const s64 start_0based = static_cast<s64>(start_1based) - 1;

  if (dash == std::string_view::npos) {
    // "chr:start" -> single 1-based position -> half-open [start-1, start).
    if (form == RegionForm::kRange) {
      throw error::Error("region '{}' must be a range (chr:start-end), not a single position", region);
    }
    return {.chromosome = chromosome, .start = start_0based, .end = static_cast<s64>(start_1based)};
  }

  const auto end_token = StripCommas(suffix.substr(dash + 1));
  if (end_token.empty()) {
    throw error::Error("region '{}' has an empty end coordinate", region);
  }
  // 1-based closed end equals the 0-based half-open exclusive end.
  const u64 end_1based = util::ParseU64(end_token);
  if (end_1based > kMaxCoord) {
    throw error::Error("region '{}' end coordinate ({}) exceeds the maximum of {}", region, end_1based, kMaxCoord);
  }
  const auto end_exclusive = static_cast<s64>(end_1based);

  if (end_exclusive < static_cast<s64>(start_1based)) {
    throw error::Error("region '{}' has end ({}) before start ({})", region, end_exclusive, start_1based);
  }

  return {.chromosome = chromosome, .start = start_0based, .end = end_exclusive};
}

std::string ContigName(const std::string_view region) {
  return std::string{region.substr(0, ContigNameLength(region))};
}

bool IsEqualContig(const std::string_view region, const std::string_view contig) {
  return region.substr(0, ContigNameLength(region)) == contig;
}

}  // namespace xoos::genome
