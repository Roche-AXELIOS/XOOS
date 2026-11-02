#include "xoos/types/str-container.h"

#include <cstddef>
#include <functional>
#include <string_view>

namespace xoos {

size_t StringHash::operator()(const std::string_view sv) const noexcept {
  return std::hash<std::string_view>{}(sv);
}

}  // namespace xoos
