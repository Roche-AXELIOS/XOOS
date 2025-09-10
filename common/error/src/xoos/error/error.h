#pragma once

#include <stdexcept>
#include <string>
#include <utility>

#include <fmt/format.h>

#include "fmt-path.h"  // IWYU pragma: keep

namespace xoos::error {

template <typename... Args>
std::runtime_error Error(fmt::format_string<Args...> fmt, Args&&... args) {
  return std::runtime_error(fmt::format(fmt, std::forward<Args>(args)...));
}

std::runtime_error Error(const std::string& fmt);
}  // namespace xoos::error
