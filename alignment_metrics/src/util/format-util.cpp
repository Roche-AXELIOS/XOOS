
#include "util/format-util.h"

namespace xoos::alignment_metrics {

std::string ToStringWithPrecision(const f64 value, const u8 precision) {
  std::ostringstream out;
  out << std::fixed << std::setprecision(precision) << value;
  return out.str();
}

}  // namespace xoos::alignment_metrics
