#pragma once

#include <fmt/format.h>

#include <filesystem>

namespace fs = std::filesystem;

template <>
struct fmt::formatter<fs::path> {
  static constexpr auto parse(format_parse_context& ctx) {  // NOLINT(readability-identifier-naming)
    return ctx.begin();
  }

  template <typename Context>
  constexpr auto format(fs::path const& path, Context& ctx) const {  // NOLINT(readability-identifier-naming)
    return format_to(ctx.out(), "{}", path.string());
  }
};
