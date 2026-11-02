#pragma once

#include <xoos/types/int.h>

namespace xoos::demux {

/// @brief Bits used to encode a single DNA base (A=0, C=1, G=2, T=3).
constexpr u32 kBitsPerEncodedBase = 2;

/// @brief Number of valid encoded DNA bases (A, C, G, T).
constexpr u8 kValidEncodedBaseCount = 4;

/// @brief Maximum 2-bit encoded base value (T = 3).
constexpr u8 kMaxEncodedBaseValue = 3;

/// @brief Bitmask for a single 2-bit encoded base.
constexpr u8 kBaseMask = 0b11;

/// @brief Number of DNA bases packed into one byte (8 / kBitsPerEncodedBase).
constexpr u32 kBasesPerByte = 4;

/// @brief Number of distinct values a byte can take (2^8); used to size byte-indexed lookup tables.
constexpr u32 kByteValueCount = 256;

}  // namespace xoos::demux
