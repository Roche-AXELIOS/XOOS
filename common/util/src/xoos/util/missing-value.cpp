#include "missing-value.h"

#include <string_view>

namespace xoos::util {

bool IsNA(const std::string_view value) {
  return value == kNA;
}

}  // namespace xoos::util
