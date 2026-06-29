#pragma once

#include <xoos/error/error.h>
#include <xoos/types/int.h>

namespace xoos::demux {

/**
 * Bitap (shift-and) algorithm implementation for approximate string matching.
 * Optimized for DNA sequence matching with support for up to 6 edit distance.
 */
template <u32 MaxDist>
class Bitap {
 public:
  // Maximum length of the query string
  static constexpr s32 kQueryWindowSize = 64;
  // Region size for sequence representation (64 bases + padding)
  // Padding is require for when we perform simd conversion on the sequence
  static constexpr s32 kBufferSize = kQueryWindowSize + 16;

  explicit Bitap(std::string_view query, bool reverse_search = false);

  Bitap(const Bitap&) = default;
  Bitap(Bitap&&) = default;
  ~Bitap() = default;

  s32 Find(const std::array<u8, kBufferSize>& alphabet_ref, s32 begin, s32 end) const;
  s32 Find(std::string_view ref, s32 begin, s32 end) const;

  s32 GetQueryLength() const { return static_cast<s32>(_query.length()); }

  s32 ForwardScan(std::string_view ref, s32 min_pos = -1, s32 max_pos = -1) const;
  s32 ReverseScan(std::string_view ref, s32 min_pos = -1, s32 max_pos = -1) const;

 private:
  // reference to the query string
  const std::string_view _query;
  // Bit mask for each character in the alphabet (A, C, G, T)
  u64 _alphabet[4] = {0ULL, 0ULL, 0ULL, 0ULL};
  // if true, the query is reversed for reverse matching
  const bool _reverse_search;
  // maximum number index start positions you can search at a time
  const s32 _max_search_size;
  // Bit mask for the most significant bit (MSB)
  const u64 _msb;

  std::array<u8, kBufferSize> ToBitapAlphabetArray(std::string_view ref, s32 begin, s32 end) const;

  static u64 ComputeMSB(size_t query_length);
};

}  // namespace xoos::demux
