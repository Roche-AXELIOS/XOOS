#include "histogram-summary.h"

namespace xoos::histogram {

std::optional<double> GetRatioPercentile(const PercentileMap& percentiles, const u64 numerator, const u64 denominator) {
  // check if the numerator and denominator are valid percentiles
  if (!percentiles.contains(numerator) || !percentiles.contains(denominator) ||
      !percentiles.at(numerator).has_value() || !percentiles.at(denominator).has_value()) {
    return std::nullopt;
  }
  // if the denominator is zero, return std::nullopt to avoid division by zero
  if (percentiles.at(denominator).value() == 0) {
    return std::nullopt;
  }
  return static_cast<double>(percentiles.at(numerator).value()) /
         static_cast<double>(percentiles.at(denominator).value());
}

}  // namespace xoos::histogram
