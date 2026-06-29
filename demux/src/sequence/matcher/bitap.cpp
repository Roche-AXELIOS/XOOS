#include "bitap.h"

#include <simd/simd-functions.h>
#include <xoos/util/sequence-functions.h>

#include <array>

/**
 * Find the marker in the read using the Bitap (shift-and) algorithm.
 *
 * This implementation is inspired by the example presented in https://microarch.org/micro53/papers/738300a951.pdf
 * (ETH Zurich), but includes several optimizations for our specific use case.
 */
namespace xoos::demux {
// Constant representing "not found" position
constexpr s32 kNotFound = std::numeric_limits<s32>::max();

/**
 * Helper function to convert a base character (A, C, G, T) to its corresponding index (0, 1, 2, 3).
 * @param base Base character (A, C, G, T)
 * @return Index (0 for A, 1 for C, 2 for G, 3 for T)
 */
inline u8 BaseToIndex(const char base) {
  const u8 bin_base = sequence::kBaseToBin[static_cast<u8>(base)];
  if (bin_base >= 4) {
    throw error::Error("Invalid character in sequence: " + std::string(1, base));
  }
  return bin_base;
}

template <u32 MaxDist>
u64 Bitap<MaxDist>::ComputeMSB(size_t query_length) {
  if (query_length > kQueryWindowSize) {
    throw error::Error("Query length must be between 0 and {}. Got {}", Bitap<0>::kQueryWindowSize, query_length);
  }
  return 1ULL << (query_length - 1);
}

/**
 * Construct Bitap matcher for approximate string matching.
 * @param query Pointer to query sequence
 * @param reverse_search Whether to perform reverse matching
 * @tparam MaxDist Maximum allowed edit distance
 */
template <u32 MaxDist>
Bitap<MaxDist>::Bitap(std::string_view query, const bool reverse_search)
    : _query(query),
      _reverse_search(reverse_search),
      _max_search_size(kQueryWindowSize - MaxDist - _query.length()),
      _msb(ComputeMSB(_query.length())) {
  if (_max_search_size < 1) {
    throw error::Error("Query length {} exceeds maximum supported length of {} accounting for edit distance {}",
                       _query.length(), kQueryWindowSize, MaxDist);
  }

  if (_reverse_search) {
    for (size_t i = query.length(); i-- > 0;) {
      _alphabet[BaseToIndex(query[i])] |= (1ULL << (query.length() - 1 - i));
    }
  } else {
    for (u32 i = 0; i < query.length(); ++i) {
      _alphabet[BaseToIndex(query[i])] |= (1ULL << i);
    }
  }

  // Invert bit masks for algorithm efficiency
  for (auto& mask : _alphabet) {
    mask = ~mask;
  }
}

/**
 * Find the query sequence in the reference string within specified range.
 * @param alphabet_ref Binary representation of reference sequence
 * @param begin Start position in reference inclusive
 * @param end End position in reference inclusive
 * @return Position of match, or -1 if not found
 */
template <u32 MaxDist>
s32 Bitap<MaxDist>::Find(const std::array<u8, kBufferSize>& alphabet_ref, const s32 begin, const s32 end) const {
  const s32 ref_length = 1 + end - begin;

  // maximum edit is determined by template parameter MaxDist
  constexpr s32 kNumLevels = MaxDist + 1;
  u64 current_state[kNumLevels];
  u64 previous_state[kNumLevels];
  s32 min_position[kNumLevels];
  for (s32 i = 0; i < kNumLevels; ++i) {
    current_state[i] = ~0ULL;
    min_position[i] = kNotFound;
  }

  // Process each character in the reference sequence
  for (s32 i = 0; i < ref_length; ++i) {
    const auto current_pattern_mask = _alphabet[alphabet_ref[i]];

    // Process exact match (distance 0)
    previous_state[0] = current_state[0];
    current_state[0] = current_state[0] << 1 | current_pattern_mask;

    if (!(current_state[0] & _msb)) {
      // Exact match found - return position immediately
      return _reverse_search ? (end - i) : begin + i;
    }

    // Using MaxDist as a template parameter ensures the loop bound is a compile-time constant.
    // This allows the compiler to optimize the loop, for example by unrolling it or eliminating bounds checks.
    // The main benefit is that the compiler can optimize fixed-size arrays and loop bounds for performance.
    for (u32 distance = 1; distance <= MaxDist; ++distance) {
      previous_state[distance] = current_state[distance];

      // Calculate operations: deletion, insertion, substitution, match
      const auto deletion = previous_state[distance - 1];
      const auto insertion = current_state[distance - 1] << 1;
      const auto substitution = previous_state[distance - 1] << 1;
      const auto match = (previous_state[distance] << 1) | current_pattern_mask;

      const auto new_state = deletion & insertion & substitution & match;

      // Record first occurrence at this edit distance
      if (min_position[distance] == kNotFound && !(new_state & _msb)) {
        min_position[distance] = i;
      }

      current_state[distance] = new_state;
    }
  }

  // Return best match within allowed edit distance
  for (u32 distance = 1; distance <= MaxDist; ++distance) {
    if (min_position[distance] != kNotFound) {
      return _reverse_search ? (end - min_position[distance]) : begin + min_position[distance];
    }
  }

  // No match found
  return -1;
}

/**
 * Transform a DNA sequence to internal representation (0,1,2,3 for A,T,G,C).
 * @param ref Pointer to reference sequence
 * @param begin Start position in reference
 * @param end End position in reference
 * @return representation of the sequence
 */
template <u32 MaxDist>
std::array<u8, Bitap<MaxDist>::kBufferSize> Bitap<MaxDist>::ToBitapAlphabetArray(std::string_view ref, const s32 begin,
                                                                                 const s32 end) const {
  // Validate range
  if (begin > end) {
    throw error::Error(
        "Bitap: begin position must be less than or equal to end position. Got begin = {} and end = {}. For matcher {}",
        begin, end, _query);
  }
  // The window size being searched
  const s32 ref_search_length = 1 + end - begin;
  if (ref_search_length > kQueryWindowSize) {
    throw error::Error("Bitap: reference search length must be between 0 and {} bases. Got length = {}. For matcher {}",
                       kQueryWindowSize, ref_search_length, _query);
  }
  // Transform sequence to edlib representation (max 64 bases + padding)
  alignas(64) std::array<u8, kBufferSize> alphabet_ref{};
  // Use SIMD function for efficient conversion of ATGC to 0,1,2,3 representation
  simd::TransformSequence(ref.data(), begin, end, alphabet_ref.data(), _reverse_search);
  return alphabet_ref;
}

/**
 * Find the query sequence in the reference string within specified range.
 * Uses a generic string representation of the reference sequence, but will A T G C to 0,1,2,3 internally.
 * @param ref Reference sequence
 * @param begin Start position in reference
 * @param end End position in reference
 * @return Position of match, or -1 if not found
 */
template <u32 MaxDist>
s32 Bitap<MaxDist>::Find(std::string_view ref, const s32 begin, const s32 end) const {
  return Find(ToBitapAlphabetArray(ref, begin, end), begin, end);
}

/**
 * @brief Finds position of a query sequence in the reference using Bitap algorithm with forward scanning.
 *
 * @param ref Reference sequence to search within
 * @param min_pos Minimum position to start searching from (defaults to 0 if < 0)
 * @param max_pos Maximum position to end searching at (defaults to end of sequence if < 0)
 * @return Position immediately after the matched sequence (match_pos + 1), or -1 if not found
 */
template <u32 MaxDist>
s32 Bitap<MaxDist>::ForwardScan(std::string_view ref, s32 min_pos, s32 max_pos) const {
  if (min_pos < 0) {
    min_pos = 0;
  }
  if (max_pos < 0) {
    max_pos = static_cast<s32>(ref.length() - 1);
  }
  // check if search positions are possible
  if (min_pos > max_pos) {
    throw error::Error(
        "Bitap::ForwardScan: min_pos must be less than or equal to max_pos. Got min_pos = {}, max_pos = {}.", min_pos,
        max_pos);
  }

  // increment amount after each search
  const s32 increment = _max_search_size;
  for (s32 begin = min_pos; begin <= max_pos; begin += increment) {
    const s32 end = begin + kQueryWindowSize - 1;
    if (const s32 pos = Find(ref, begin, std::min(end, max_pos)); pos != -1) {
      // return position after the found sequence so add 1
      return pos + 1;
    }
  }
  return -1;
}

/**
 * @brief Finds position of a query sequence in the reference using Bitap algorithm with reverse scanning.
 *
 * Performs a greedy search from 3' to 5' direction (right to left) by dividing the search range
 * into chunks of size `_max_search_size` and searching each chunk sequentially. Returns the position
 * of the first match found.
 *
 * @param ref Reference sequence to search within
 * @param min_pos Minimum position to start searching from (defaults to 0 if < 0)
 * @param max_pos Maximum position to end searching at (defaults to end of sequence if < 0)
 * @return Position of the matched sequence, or -1 if not found
 */
template <u32 MaxDist>
s32 Bitap<MaxDist>::ReverseScan(std::string_view ref, s32 min_pos, s32 max_pos) const {
  if (min_pos < 0) {
    min_pos = 0;
  }
  if (max_pos < 0) {
    max_pos = static_cast<s32>(ref.length() - 1);
  }
  // check if search positions are possible
  if (min_pos > max_pos) {
    throw error::Error(
        "Bitap::ReverseScan: min_pos must be less than or equal to max_pos. Got min_pos = {}, max_pos = {}.", min_pos,
        max_pos);
  }

  // increment amount after each search
  const s32 increment = _max_search_size;
  for (s32 end = max_pos; end >= min_pos; end -= increment) {
    const s32 begin = end - kQueryWindowSize + 1;
    if (const s32 pos = Find(ref, std::max(begin, min_pos), end); pos != -1) {
      return pos;
    }
  }
  return -1;
}

}  // namespace xoos::demux

// Explicit template instantiations
// This ensures that the compiler generates code for these specific template parameters
// In contrast to defining the function in the header file, this approach can reduce compilation times
// and binary size when the template is used with a limited set of parameters. Add more as needed.
template class xoos::demux::Bitap<2>;
template class xoos::demux::Bitap<3>;
template class xoos::demux::Bitap<4>;
