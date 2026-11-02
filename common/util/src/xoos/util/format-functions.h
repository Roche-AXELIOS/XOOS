#pragma once

#include <concepts>
#include <string>

#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/util/missing-value.h>

namespace xoos::util {

// Default number of decimal places used when formatting floating-point values for TSV output.
constexpr u8 kDefaultPrecision = 2;

// Formats a f64 to a string with the given precision to provide output consistency.
std::string ToStringWithPrecision(f64 value, u8 precision = kDefaultPrecision);

// Formats a ratio as a numeric percentage. Output columns and row names carry the unit, so no percent sign is added.
template <typename T>
  requires std::convertible_to<T, f64>
std::string ToPercentageWithPrecision(const T numerator, const T denominator, const u8 precision = kDefaultPrecision) {
  // TODO: add a more tolerant check for zero denominator
  if (denominator == 0) {
    return std::string{kNA};
  }
  return ToStringWithPrecision(static_cast<f64>(numerator) / static_cast<f64>(denominator) * 100.0, precision);
}

std::string ToPhredOrNA(u64 numerator, u64 denominator, u8 precision);

}  // namespace xoos::util
