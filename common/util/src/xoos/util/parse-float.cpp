#include "parse-float.h"

#include <concepts>
#include <cstddef>
#include <locale>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>

#include <xoos/error/error.h>
#include <xoos/types/float.h>

namespace xoos::util {

namespace {

/// Map floating-point types to their display names for error messages.
template <typename T>
constexpr std::string_view TypeName() {
  if constexpr (std::is_same_v<T, f32>) {
    return "f32";
  } else if constexpr (std::is_same_v<T, f64>) {
    return "f64";
  }
}

/// Convert a string to a floating-point value of type T using the matching std library call.
/// pos is set to the number of characters consumed so the caller can reject trailing characters.
template <std::floating_point T>
T StdParse(const std::string& value, std::size_t& pos) {
  if constexpr (std::is_same_v<T, f32>) {
    return std::stof(value, &pos);
  } else {
    return std::stod(value, &pos);
  }
}

}  // namespace

template <std::floating_point T>
T ParseFloat(const std::string_view value) {
  constexpr auto kName = TypeName<T>();
  if (value.empty()) {
    throw error::Error("Cannot parse empty string as {}", kName);
  }
  // std::stod / std::stof skip leading whitespace; reject it explicitly so behaviour matches ParseUnsigned /
  // ParseSigned.
  if (std::isspace(value.front(), std::locale::classic())) {
    throw error::Error("'{}' is not a valid {}", value, kName);
  }
  const std::string str(value);
  std::size_t pos = 0;
  T result{};
  try {
    result = StdParse<T>(str, pos);
  } catch (const std::out_of_range&) {
    throw error::Error("'{}' is out of range for {}", value, kName);
  } catch (const std::invalid_argument&) {
    throw error::Error("'{}' is not a valid {}", value, kName);
  }
  // Reject trailing characters (e.g. "1.2abc") so malformed fields fail fast rather than being silently truncated.
  if (pos != str.size()) {
    throw error::Error("'{}' is not a valid {}", value, kName);
  }
  return result;
}

// Explicit instantiations for all supported types.
template f32 ParseFloat<f32>(std::string_view);
template f64 ParseFloat<f64>(std::string_view);

}  // namespace xoos::util
