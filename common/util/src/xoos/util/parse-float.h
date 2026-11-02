#pragma once
#include <concepts>
#include <string_view>

#include <xoos/types/float.h>

namespace xoos::util {

/// Parse a string as a floating-point value of type T.
/// The entire string must be consumed: leading whitespace and trailing characters are rejected so that
/// inputs such as "1.2abc" or " 1.2" fail rather than being silently truncated by std::stod / std::stof.
/// Throws on empty input, invalid characters, or overflow beyond T's range.
template <std::floating_point T>
T ParseFloat(std::string_view value);

// Convenience aliases preserving the original API.

inline f64 ParseF64(const std::string_view value) {
  return ParseFloat<f64>(value);
}

inline f32 ParseF32(const std::string_view value) {
  return ParseFloat<f32>(value);
}

}  // namespace xoos::util
