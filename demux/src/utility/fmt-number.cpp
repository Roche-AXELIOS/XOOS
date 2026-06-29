#include "utility/fmt-number.h"

#include <fmt/format.h>

#include <vector>

namespace xoos::demux {

std::string FormatWithMetricSuffix(double value) {
  const auto metric_suffixes = std::vector<std::string>{"k", "M", "G", "T", "P", "E"};

  auto current_value = value;
  auto current_suffix = std::string("");
  for (const auto& suffix : metric_suffixes) {
    if (current_value < 1000.0) {
      break;
    }
    current_value /= 1000.0;
    current_suffix = suffix;
  }
  return fmt::format("{:.1f}{}", current_value, current_suffix);
}

}  // namespace xoos::demux
