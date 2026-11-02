#include "metrics.h"

#include <cmath>
#include <string>

#include <fmt/format.h>

#include <xoos/types/float.h>
#include <xoos/types/int.h>

#include "metrics-constants.h"

namespace xoos::demux {

f32 CalculatePercentage(const u64 numerator, const u64 denominator) {
  return denominator == 0 ? NAN : static_cast<f32>(static_cast<f64>(numerator) / static_cast<f64>(denominator) * 100.0);
}

std::string FormatPercentage(const u64 numerator, const u64 denominator, const s32 precision) {
  const f32 percentage = CalculatePercentage(numerator, denominator);
  if (std::isnan(percentage)) {
    return std::string(kNA);
  }
  return fmt::format("{:.{}f}", percentage, precision);
}

std::string FormatPercentage(const u64 numerator, const u64 denominator) {
  return FormatPercentage(numerator, denominator, kPercentagePrecision);
}

}  // namespace xoos::demux
