#include "format-functions.h"

#include <iomanip>
#include <sstream>

#include <xoos/util/math.h>

namespace xoos::util {

std::string ToStringWithPrecision(const f64 value, const u8 precision) {
  std::ostringstream out;
  out << std::fixed << std::setprecision(precision) << value;
  return out.str();
}

std::string ToPhredOrNA(const u64 numerator, const u64 denominator, const u8 precision) {
  if (denominator == 0) {
    return std::string{kNA};
  }
  return ToStringWithPrecision(math::ErrorRateToPhred(static_cast<f64>(numerator) / static_cast<f64>(denominator)),
                               precision);
}

}  // namespace xoos::util
