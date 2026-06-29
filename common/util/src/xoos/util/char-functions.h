#pragma once

#include <cstdint>

namespace xoos::character {

constexpr bool IsSpace(uint8_t c) noexcept {
  return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\v' || c == '\f';
}

constexpr bool IsNonTrivialSpace(uint8_t c) noexcept {
  return c == '\t' || c == '\n' || c == '\r' || c == '\v' || c == '\f';
}

}  // namespace xoos::character
