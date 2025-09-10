#pragma once

#include <filesystem>

#include <fmt/format.h>

namespace fs = std::filesystem;

template <>
struct fmt::formatter<fs::path> {
  static constexpr auto parse(format_parse_context& ctx) {  // NOLINT(readability-identifier-naming)
    return ctx.begin();
  }

  template <typename Context>
  constexpr auto format(const fs::path& path, Context& ctx) const {  // NOLINT(readability-identifier-naming)
    return fmt::format_to(ctx.out(), "{}", path.string());
  }
};
