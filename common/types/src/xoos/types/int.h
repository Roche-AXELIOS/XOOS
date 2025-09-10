#pragma once

#include <concepts>
#include <cstdint>
#include <type_traits>

namespace xoos {

using u8 = std::uint8_t;
using u16 = std::uint16_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using s8 = std::int8_t;
using s16 = std::int16_t;
using s32 = std::int32_t;
using s64 = std::int64_t;

// Using std::make_signed and std::make_unsigned to convert between signed and unsigned integers
// is more reliable and portable than using static_cast by itself
// https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2018/p0907r4.html
// this ensures that the conversion is done for the equivalent bit depth
template <std::unsigned_integral T>
constexpr auto ToSigned(const T value) {
  return static_cast<std::make_signed_t<T>>(value);
}

template <std::signed_integral T>
constexpr auto ToUnsigned(const T value) {
  return static_cast<std::make_unsigned_t<T>>(value);
}

}  // namespace xoos
