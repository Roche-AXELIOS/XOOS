#include "likelihood/likelihood-flags.h"

#include <algorithm>
#include <iterator>
#include <ranges>
#include <stdexcept>
#include <string>

#include <xoos/log/logging.h>
#include <xoos/types/vec.h>

namespace xoos::cnc {
vec<std::string> StringToFlags(const std::string& str) {
  vec<std::string> ret;
  for (const auto& word : std::views::split(str, std::string_view{";"})) {
    std::string flag;
    std::ranges::copy(word, std::back_inserter(flag));
    if (std::ranges::find(kLikelihoodFlagsAll, flag) == kLikelihoodFlagsAll.end()) {
      Logging::Error("invalid flag: {}!", flag);
      throw std::runtime_error("flag parsing error");
    }
    if (flag != kLikelihoodFlagNone) {
      ret.emplace_back(flag);
    }
  }
  return ret;
}

std::string FlagsToString(const vec<std::string>& flags) {
  if (flags.empty()) {
    return kLikelihoodFlagNone;
  }
  std::string flags_str{};
  for (const auto& flag : flags) {
    if (flag != kLikelihoodFlagNone) {
      if (!flags_str.empty()) {
        flags_str += ";";
      }
      flags_str += flag;
    }
  }
  if (flags_str.empty()) {
    return kLikelihoodFlagNone;
  }
  return flags_str;
}

}  // namespace xoos::cnc
