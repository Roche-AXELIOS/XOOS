#include "bitap.h"

#include <algorithm>

#include <xoos/util/sequence-functions.h>

#include "sequence/matcher/match-position-states.h"
#include "simd/simd-functions.h"

namespace xoos::demux {
// Constant representing "not found" position
constexpr s32 kNotFound = std::numeric_limits<s32>::max();
constexpr s32 kInclusiveRangeAdjustment = 1;
constexpr s32 kMinimumSearchWindowSize = 1;
constexpr size_t kBitapAlphabetAlignmentBytes = 64;

template <u32 MaxDist>
u64 Bitap<MaxDist>::ComputeMSB(size_t query_length) {
  if (query_length == 0 || query_length > kQueryWindowSize) {
    throw error::Error("Query length must be between 1 and {}. Got {}", Bitap<0>::kQueryWindowSize, query_length);
  }
  return 1ULL << (query_length - kInclusiveRangeAdjustment);
}

template <u32 MaxDist>
Bitap<MaxDist>::Bitap(const std::string_view query, const SearchDirection direction)
    : _query(query),
      _direction(direction),
      _max_search_size(kQueryWindowSize - static_cast<s32>(MaxDist + _query.length())),
      _msb(ComputeMSB(_query.length())) {
  if (_max_search_size < kMinimumSearchWindowSize) {
    throw error::Error("Query length {} exceeds maximum supported length of {} accounting for edit distance {}",
                       _query.length(),
                       kQueryWindowSize,
                       MaxDist);
  }

  if (_direction == SearchDirection::kBackward) {
    for (size_t i = 0; i < query.length(); ++i) {
      _alphabet[sequence::BaseToIndex(query[query.length() - 1 - i])] |= (1ULL << i);
    }
  } else {
    for (u32 i = 0; i < query.length(); ++i) {
      _alphabet[sequence::BaseToIndex(query[i])] |= (1ULL << i);
    }
  }

  // Invert bit masks for algorithm efficiency
  for (auto& mask : _alphabet) {
    mask = ~mask;
  }
}

template <u32 MaxDist>
s32 Bitap<MaxDist>::Find(const std::array<u8, kBufferSize>& alphabet_ref, const s32 begin, const s32 end) const {
  return SearchWindow<MatchPolicy::kFirst>(alphabet_ref, begin, end);
}

template <u32 MaxDist>
s32 Bitap<MaxDist>::FindWidest(std::string_view ref, const s32 begin, const s32 end) const {
  return FindWidestMatch(ref, begin, end).position;
}

template <u32 MaxDist>
typename Bitap<MaxDist>::Match Bitap<MaxDist>::FindWidestMatch(std::string_view ref,
                                                               const s32 begin,
                                                               const s32 end) const {
  return SearchWindowMatch<MatchPolicy::kWidest>(ToBitapAlphabetArray(ref, begin, end), begin, end);
}

template <u32 MaxDist>
template <typename Bitap<MaxDist>::MatchPolicy Policy>
void Bitap<MaxDist>::UpdateDistanceLevels(const u64 current_pattern_mask,
                                          const s32 position,
                                          std::array<u64, kNumLevels>& current_state,
                                          std::array<u64, kNumLevels>& previous_state,
                                          std::array<s32, kNumLevels>& best_position) const {
  // MaxDist as a template parameter ensures the loop bound is a compile-time constant,
  // allowing the compiler to unroll and eliminate bounds checks.
  for (u32 distance = 1; distance <= MaxDist; ++distance) {
    previous_state[distance] = current_state[distance];

    // Calculate operations: deletion, insertion, substitution, match
    const auto deletion = previous_state[distance - 1];
    const auto insertion = current_state[distance - 1] << 1;
    const auto substitution = previous_state[distance - 1] << 1;
    const auto match = (previous_state[distance] << 1) | current_pattern_mask;

    const auto new_state = deletion & insertion & substitution & match;
    const bool is_match = (new_state & _msb) == 0;

    if constexpr (Policy == MatchPolicy::kFirst) {
      // Latch on first occurrence at this edit distance.
      // NOSONAR(M23_002) -- flagged both always-true and always-false across template instantiations; the
      // condition is genuinely runtime-dependent for each concrete Bitap<MaxDist, Policy>.
      if ((best_position[distance] == kNotFound) && is_match) {  // NOSONAR(M23_002)
        best_position[distance] = position;
      }
    } else if (is_match) {
      // Overwrite on every hit — keeps the latest (widest) position.
      best_position[distance] = position;
    }

    current_state[distance] = new_state;
  }
}

template <u32 MaxDist>
template <typename Bitap<MaxDist>::MatchPolicy Policy>
s32 Bitap<MaxDist>::SearchWindow(const std::array<u8, kBufferSize>& alphabet_ref,
                                 const s32 begin,
                                 const s32 end) const {
  return SearchWindowMatch<Policy>(alphabet_ref, begin, end).position;
}

template <u32 MaxDist>
template <typename Bitap<MaxDist>::MatchPolicy Policy>
typename Bitap<MaxDist>::Match Bitap<MaxDist>::SearchWindowMatch(const std::array<u8, kBufferSize>& alphabet_ref,
                                                                 const s32 begin,
                                                                 const s32 end) const {
  const s32 ref_length = kInclusiveRangeAdjustment + end - begin;

  std::array<u64, kNumLevels> current_state{};
  std::array<u64, kNumLevels> previous_state{};
  std::array<s32, kNumLevels> best_position{};
  current_state.fill(~0ULL);
  best_position.fill(kNotFound);

  // Process each character in the reference sequence
  for (auto i = 0; i < ref_length; ++i) {
    const auto current_pattern_mask = _alphabet[alphabet_ref[static_cast<size_t>(i)]];

    // Process exact match (distance 0)
    previous_state[0] = current_state[0];
    current_state[0] = current_state[0] << 1 | current_pattern_mask;

    if ((current_state[0] & _msb) == 0) {
      if constexpr (Policy == MatchPolicy::kFirst) {
        // Early return on first exact match.
        return {.position = _direction == SearchDirection::kBackward ? (end - i) : begin + i, .edit_distance = 0};
      } else {
        best_position[0] = i;
      }
    }

    UpdateDistanceLevels<Policy>(current_pattern_mask, i, current_state, previous_state, best_position);
  }

  // Return the position at the lowest edit distance that had a match. For kFirst, distance zero
  // returned early above; starting from zero keeps this method's score contract uniform.
  for (u32 distance = 0; distance <= MaxDist; ++distance) {
    if (best_position[distance] != kNotFound) {
      return {.position = _direction == SearchDirection::kBackward ? (end - best_position[distance])
                                                                   : begin + best_position[distance],
              .edit_distance = distance};
    }
  }

  return {};
}

template <u32 MaxDist>
std::array<u8, Bitap<MaxDist>::kBufferSize> Bitap<MaxDist>::ToBitapAlphabetArray(const std::string_view ref,
                                                                                 const s32 begin,
                                                                                 const s32 end) const {
  // Validate range
  if (begin > end) {
    throw error::Error(
        "Bitap: begin position must be less than or equal to end position. Got begin = {} and end = {}. For matcher {}",
        begin,
        end,
        _query);
  }

  const s32 ref_search_length = kInclusiveRangeAdjustment + end - begin;
  if (ref_search_length > kQueryWindowSize) {
    throw error::Error("Bitap: reference search length must be between 0 and {} bases. Got length = {}. For matcher {}",
                       kQueryWindowSize,
                       ref_search_length,
                       _query);
  }

  alignas(kBitapAlphabetAlignmentBytes) std::array<u8, kBufferSize> alphabet_ref{};
  // Use SIMD function for efficient conversion of ATGC to 0,1,2,3 representation
  simd::TransformSequence(ref.data(), begin, end, alphabet_ref.data(), _direction == SearchDirection::kBackward);
  return alphabet_ref;
}

template <u32 MaxDist>
s32 Bitap<MaxDist>::Find(std::string_view ref, const s32 begin, const s32 end) const {
  return Find(ToBitapAlphabetArray(ref, begin, end), begin, end);
}

template <u32 MaxDist>
std::pair<s32, s32> Bitap<MaxDist>::FindStartEnd(std::string_view ref,
                                                 const s32 begin,
                                                 const s32 end,
                                                 const Bitap<MaxDist>& reverse_bitap) const {
  const auto query_len = static_cast<s32>(_query.length());
  const auto max_dist = static_cast<s32>(MaxDist);

  if (_direction == SearchDirection::kForward) {
    // Forward finds end, reverse recovers start
    const s32 end_pos = Find(ref, begin, end);
    if (end_pos == kNoMatchPosition) {
      return {kNoMatchPosition, kNoMatchPosition};
    }
    // Use FindWidest on the reverse pass to recover the leftmost (widest) start position.
    // Find (kFirst) would return the rightmost start, giving the narrowest span — but for
    // indel-containing matches the widest span better represents the true alignment extent.
    const s32 rev_begin = std::max(begin, end_pos - query_len - max_dist + 1);
    const s32 start_pos = reverse_bitap.FindWidest(ref, rev_begin, end_pos);
    return {start_pos, end_pos};
  }

  // Backward finds start, forward (reverse_bitap) recovers end
  const s32 start_pos = Find(ref, begin, end);
  if (start_pos == kNoMatchPosition) {
    return {kNoMatchPosition, kNoMatchPosition};
  }
  // Symmetrically, use FindWidest on the forward recovery pass for the rightmost (widest) end.
  const s32 fwd_end = std::min(end, start_pos + query_len + max_dist);
  const s32 end_pos = reverse_bitap.FindWidest(ref, start_pos, fwd_end);
  // if we find start_pos we always find end_pos
  return {start_pos, end_pos};
}

template <u32 MaxDist>
s32 Bitap<MaxDist>::ForwardScan(const std::string_view ref, s32 min_pos, s32 max_pos) const {
  min_pos = std::max(min_pos, 0);
  if (max_pos < 0) {
    max_pos = static_cast<s32>(ref.length() - kInclusiveRangeAdjustment);
  }
  // check if search positions are possible
  if (min_pos > max_pos) {
    throw error::Error(
        "Bitap::ForwardScan: min_pos must be less than or equal to max_pos. Got min_pos = {}, max_pos = {}.",
        min_pos,
        max_pos);
  }

  // increment amount after each search
  const s32 increment = _max_search_size;
  for (s32 begin = min_pos; begin <= max_pos; begin += increment) {
    const s32 end = begin + kQueryWindowSize - kInclusiveRangeAdjustment;
    if (const s32 pos = Find(ref, begin, std::min(end, max_pos)); pos != kNoMatchPosition) {
      // return position after the found sequence so add 1
      return pos + kPositionAfterMatchAdjustment;
    }
  }
  return kNoMatchPosition;
}

template <u32 MaxDist>
s32 Bitap<MaxDist>::ReverseScan(const std::string_view ref, s32 min_pos, s32 max_pos) const {
  min_pos = std::max(min_pos, 0);
  if (max_pos < 0) {
    max_pos = static_cast<s32>(ref.length() - kInclusiveRangeAdjustment);
  }
  // check if search positions are possible
  if (min_pos > max_pos) {
    throw error::Error(
        "Bitap::ReverseScan: min_pos must be less than or equal to max_pos. Got min_pos = {}, max_pos = {}.",
        min_pos,
        max_pos);
  }

  // increment amount after each search
  const s32 increment = _max_search_size;
  for (s32 end = max_pos; end >= min_pos; end -= increment) {
    const s32 begin = end - kQueryWindowSize + kInclusiveRangeAdjustment;
    if (const s32 pos = Find(ref, std::max(begin, min_pos), end); pos != kNoMatchPosition) {
      return pos;
    }
  }
  return kNoMatchPosition;
}

}  // namespace xoos::demux

constexpr auto kEditDistance2 = 2;
constexpr auto kEditDistance3 = 3;
constexpr auto kEditDistance4 = 4;

// Explicit template instantiations
// This ensures that the compiler generates code for these specific template parameters
// In contrast to defining the function in the header file, this approach can reduce compilation times
// and binary size when the template is used with a limited set of parameters. Add more as needed.
// NOSONAR(M23_051,M23_053) -- distinct explicit instantiations of a class template are not ODR/duplicate
// definitions; the analyzer misclassifies each Bitap<N> as a redefinition of the same entity.
template class xoos::demux::Bitap<kEditDistance2>;  // NOSONAR(M23_051,M23_053)
template class xoos::demux::Bitap<kEditDistance3>;  // NOSONAR(M23_051,M23_053)
template class xoos::demux::Bitap<kEditDistance4>;  // NOSONAR(M23_051,M23_053)
