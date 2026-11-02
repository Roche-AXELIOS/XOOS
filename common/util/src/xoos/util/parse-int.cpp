#include "parse-int.h"

#include <algorithm>
#include <cctype>
#include <concepts>
#include <limits>
#include <locale>
#include <stdexcept>
#include <string>
#include <string_view>

#include <xoos/error/error.h>
#include <xoos/types/int.h>

namespace xoos::util {

namespace {

/// Map integer types to their display names for error messages.
template <typename T>
constexpr std::string_view TypeName() {
  if constexpr (std::is_same_v<T, u8>) {
    return "u8";
  } else if constexpr (std::is_same_v<T, u16>) {
    return "u16";
  } else if constexpr (std::is_same_v<T, u32>) {
    return "u32";
  } else if constexpr (std::is_same_v<T, u64>) {
    return "u64";
  } else if constexpr (std::is_same_v<T, s8>) {
    return "s8";
  } else if constexpr (std::is_same_v<T, s16>) {
    return "s16";
  } else if constexpr (std::is_same_v<T, s32>) {
    return "s32";
  } else if constexpr (std::is_same_v<T, s64>) {
    return "s64";
  }
}

void ValidateUnsignedDigits(const std::string_view value, const std::string_view type_name) {
  if (value.empty()) {
    throw error::Error("Cannot parse empty string as {}", type_name);
  }
  if (!std::ranges::all_of(value, [](const char c) { return std::isdigit(c, std::locale::classic()); })) {
    throw error::Error("'{}' is not a valid {}: contains non-digit characters", value, type_name);
  }
}

void ValidateSignedDigits(const std::string_view value, const std::string_view type_name) {
  if (value.empty()) {
    throw error::Error("Cannot parse empty string as {}", type_name);
  }
  const auto digits = value.front() == '-' ? value.substr(1) : value;
  if (digits.empty() ||
      !std::ranges::all_of(digits, [](const char c) { return std::isdigit(c, std::locale::classic()); })) {
    throw error::Error("'{}' is not a valid {}: contains non-digit characters", value, type_name);
  }
}

}  // namespace

template <std::unsigned_integral T>
T ParseUnsigned(const std::string_view value) {
  constexpr auto kName = TypeName<T>();
  ValidateUnsignedDigits(value, kName);
  try {
    const auto result = std::stoull(std::string(value));
    if constexpr (sizeof(T) < sizeof(u64)) {
      if (result > std::numeric_limits<T>::max()) {
        throw error::Error("'{}' is out of range for {}", value, kName);
      }
    }
    return static_cast<T>(result);
  } catch (const std::out_of_range&) {
    throw error::Error("'{}' is out of range for {}", value, kName);
  }
}

template <std::signed_integral T>
T ParseSigned(const std::string_view value) {
  constexpr auto kName = TypeName<T>();
  ValidateSignedDigits(value, kName);
  try {
    const auto result = std::stoll(std::string(value));
    if constexpr (sizeof(T) < sizeof(s64)) {
      if (result < std::numeric_limits<T>::min() || result > std::numeric_limits<T>::max()) {
        throw error::Error("'{}' is out of range for {}", value, kName);
      }
    }
    return static_cast<T>(result);
  } catch (const std::out_of_range&) {
    throw error::Error("'{}' is out of range for {}", value, kName);
  }
}

// Explicit instantiations for all supported types.
template u8 ParseUnsigned<u8>(std::string_view);
template u16 ParseUnsigned<u16>(std::string_view);
template u32 ParseUnsigned<u32>(std::string_view);
template u64 ParseUnsigned<u64>(std::string_view);
template s8 ParseSigned<s8>(std::string_view);
template s16 ParseSigned<s16>(std::string_view);
template s32 ParseSigned<s32>(std::string_view);
template s64 ParseSigned<s64>(std::string_view);

}  // namespace xoos::util
