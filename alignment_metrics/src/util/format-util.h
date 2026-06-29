#pragma once

#include <concepts>
#include <string>

#include <xoos/types/float.h>
#include <xoos/types/fs.h>

#include "xoos/types/int.h"

namespace xoos::alignment_metrics {

// A constant string to represent "not applicable" values in TSV output.
static const std::string kNotApplicable = "NA";

// Formats a f64 to a string with the given precision to provide output consistency.
std::string ToStringWithPrecision(f64 value, u8 precision = 6);

template <typename T>
  requires std::convertible_to<T, f64>
std::string ToPercentageWithPrecision(const T numerator, const T denominator, const u8 precision = 6) {
  // TODO: add a more tolerant check for zero denominator
  if (denominator == 0) {
    return kNotApplicable;
  }
  return ToStringWithPrecision(static_cast<f64>(numerator) / static_cast<f64>(denominator) * 100.0, precision) + '%';
}

}  // namespace xoos::alignment_metrics
