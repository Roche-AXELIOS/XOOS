#pragma once

#include <cmath>
#include <concepts>
#include <iomanip>
#include <sstream>
#include <string>
#include <tuple>

#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/util/missing-value.h>

namespace xoos::cnc {

/// Parse a 1-based closed region string "chr:start-end" and return (chromosome, start, end)
/// as 0-based half-open coordinates.
std::tuple<std::string, size_t, size_t> ParseRegionString(const std::string& r);

/// Spelling used for missing values in .seg and coverage BED files written before util::kNA was adopted.
inline constexpr std::string_view kLegacyMissingValue = "nan";

/**
 * @brief Check whether a field read from a CNC output file represents a missing value.
 *
 * Accepts the standard util::kNA token plus kLegacyMissingValue, so a .seg or coverage BED produced by an earlier
 * version still loads. This is read-side only: writers must emit util::kNA via FormatOrNA or FloatAsString. The
 * legacy spelling is not recognized by util::IsNA, which deliberately supports the standard token alone.
 *
 * @param value field to test
 * @return true if the field is the standard token or the legacy spelling
 */
bool IsMissingValue(std::string_view value);

constexpr s32 kDefaultFloatPrecision = 6;

/// Format a float for TSV/BED output. A NaN means the value could not be computed, so it is written as the standard
/// missing-value token rather than the platform-dependent "nan"/"-nan" spelling.
std::string FloatAsString(std::floating_point auto val, s32 precision = kDefaultFloatPrecision) {
  if (std::isnan(val)) {
    return std::string{util::kNA};
  }
  std::stringstream ss;
  ss << std::fixed << std::setprecision(precision) << static_cast<f64>(val);
  return ss.str();
}

std::string FloatAsIntString(std::floating_point auto val) {
  return FloatAsString(val, 0);
}

}  // namespace xoos::cnc
