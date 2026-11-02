#pragma once

#include <array>
#include <string_view>
#include <utility>

#include <xoos/error/error.h>
#include <xoos/types/int.h>

#include "sequence/matcher/match-position-states.h"
#include "sequence/matcher/search-direction.h"

namespace xoos::demux {

/**
 * @brief Bitap (shift-and) algorithm for approximate DNA string matching.
 *
 * Supports up to 6 edit distance (template parameter MaxDist). Inspired by the
 * implementation presented in https://microarch.org/micro53/papers/738300a951.pdf
 * (ETH Zurich), with several optimizations for our specific use case.
 *
 * SearchDirection controls which edge of the match position is returned:
 *   - kForward  → right edge (e.g. finding a 5' trim point after the adapter)
 *   - kBackward → left edge  (e.g. finding a 3' trim point before the adapter)
 *
 * ForwardScan/ReverseScan control the window iteration order across the reference
 * and can be used regardless of the initialized SearchDirection. This is useful for
 * conservative trimming — e.g. if you prefer to trim a 5' end closer to the insert
 * rather than the left end, you would use kForward with a ReverseScan.
 *
 * @tparam MaxDist Maximum allowed edit distance (substitution, insertion, deletion).
 */
template <u32 MaxDist>
class Bitap {
 public:
  /// @brief Maximum length of the query string (limited by u64 bitmask width).
  // NOSONAR(M23_051,M23_053) -- static constexpr member of a class template; each instantiation defines its
  // own member and this is not an ODR/duplicate-definition violation.
  static constexpr s32 kQueryWindowSize = 64;  // NOSONAR(M23_051,M23_053)

  /// @brief Region size for sequence representation (64 bases + SIMD padding).
  static constexpr s32 kBufferSize = kQueryWindowSize + 16;  // NOSONAR(M23_051,M23_053)

  /// @brief Number of symbols in the DNA alphabet (A, C, G, T).
  static constexpr size_t kAlphabetSize = 4;

  /**
   * @brief A match edge and its exact minimum edit distance within one Bitap window.
   *
   * The position follows the instance's SearchDirection convention: a forward matcher returns
   * the inclusive right edge; a backward matcher returns the inclusive left edge.
   */
  struct Match {
    s32 position{kNoMatchPosition};
    u32 edit_distance{MaxDist + 1};

    bool Found() const {
      return position != kNoMatchPosition;
    }
  };

  /**
   * @brief Construct a Bitap matcher for approximate string matching.
   *
   * @param query Query sequence to search for (length must be < kQueryWindowSize - MaxDist).
   * @param direction Search direction: kForward returns the right edge, kBackward returns
   *                  the left edge of the match.
   * @throws error::Error If the query length exceeds the supported window size.
   */
  explicit Bitap(std::string_view query, SearchDirection direction = SearchDirection::kForward);

  Bitap(const Bitap&) = default;
  Bitap(Bitap&&) noexcept = default;
  Bitap& operator=(Bitap&&) = delete;
  Bitap& operator=(const Bitap&) = delete;
  ~Bitap() = default;

  /**
   * @brief Find the query in a pre-converted reference alphabet array.
   *
   * Returns the first (earliest) match position at the best edit distance.
   *
   * @param alphabet_ref Reference sequence pre-converted to 0-3 representation via ToBitapAlphabetArray.
   * @param begin Start position in the original reference (inclusive).
   * @param end End position in the original reference (inclusive).
   * @return Position of match, or kNoMatchPosition if not found.
   */
  s32 Find(const std::array<u8, kBufferSize>& alphabet_ref, s32 begin, s32 end) const;

  /**
   * @brief Find the query in a raw DNA reference string within the specified range.
   *
   * Internally converts A/C/G/T to 0/1/2/3 and delegates to the array-based overload.
   *
   * @param ref Reference sequence (ASCII DNA).
   * @param begin Start position in reference (inclusive).
   * @param end End position in reference (inclusive).
   * @return Position of match, or kNoMatchPosition if not found.
   */
  s32 Find(std::string_view ref, s32 begin, s32 end) const;

  /**
   * @brief Find the widest (rightmost/latest) match at the best edit distance.
   *
   * Unlike Find(), which returns the first match and exits early, this variant scans the entire
   * [begin, end] window and returns the last position that matches at the lowest edit distance.
   * This yields the widest alignment span — useful for conservative 3' trimming where you want
   * the match that consumes the most reference bases.
   *
   * Because this always processes the full window (no early exit), callers should keep the
   * search range tight — e.g. query_length + MaxDist bases beyond the expected match region.
   *
   * @param ref Reference sequence (ASCII DNA).
   * @param begin Start position in reference (inclusive).
   * @param end End position in reference (inclusive).
   * @return Position of widest match, or kNoMatchPosition if not found.
   */
  s32 FindWidest(std::string_view ref, s32 begin, s32 end) const;

  /**
   * @brief Find the widest match and report its minimum edit distance.
   *
   * This is the score-aware companion to FindWidest. It is intended for callers that combine
   * several 64bp windows and must compare their candidates without weakening the edit threshold.
   */
  Match FindWidestMatch(std::string_view ref, s32 begin, s32 end) const;

  /**
   * @brief Find both start and end positions of the best match in the reference.
   *
   * Uses two passes: this instance locates one edge, then reverse_bitap (opposite direction,
   * same query) recovers the other edge within a narrowed window of query_len + MaxDist.
   *
   * - kForward instance: forward pass finds end, reverse_bitap recovers start.
   * - kBackward instance: backward pass finds start, reverse_bitap recovers end.
   *
   * @param ref Reference sequence (ASCII DNA).
   * @param begin Start position in reference (inclusive).
   * @param end End position in reference (inclusive).
   * @param reverse_bitap A Bitap constructed with the opposite direction on the same query.
   * @return {start, end} (both inclusive), or {kNoMatchPosition, kNoMatchPosition} if no match.
   */
  std::pair<s32, s32> FindStartEnd(std::string_view ref, s32 begin, s32 end, const Bitap& reverse_bitap) const;

  /// @brief Returns the query length as s32.
  s32 GetQueryLength() const {
    return static_cast<s32>(_query.length());
  }

  /// @brief Returns the immutable query sequence used by this matcher.
  std::string_view Query() const {
    return _query;
  }

  /**
   * @brief Scan the reference left-to-right in chunks, returning the first match found.
   *
   * @param ref Reference sequence to search within.
   * @param min_pos Minimum position to start searching from (clamped to 0 if negative).
   * @param max_pos Maximum position to end searching at (clamped to ref end if negative).
   * @return Position immediately after the matched sequence (match_pos + 1), or kNoMatchPosition.
   * @throws error::Error If min_pos > max_pos after clamping.
   */
  s32 ForwardScan(std::string_view ref, s32 min_pos, s32 max_pos) const;

  /**
   * @brief Scan the reference right-to-left in chunks, returning the first match found.
   *
   * Performs a greedy search from 3' to 5' direction by dividing the search range
   * into chunks of size _max_search_size and searching each chunk sequentially.
   *
   * @param ref Reference sequence to search within.
   * @param min_pos Minimum position to start searching from (clamped to 0 if negative).
   * @param max_pos Maximum position to end searching at (clamped to ref end if negative).
   * @return Position of the matched sequence, or kNoMatchPosition if not found.
   * @throws error::Error If min_pos > max_pos after clamping.
   */
  s32 ReverseScan(std::string_view ref, s32 min_pos, s32 max_pos) const;

 private:
  /// @brief The query sequence to search for.
  const std::string_view _query;

  /// @brief Bit mask for each character in the DNA alphabet (A=0, C=1, G=2, T=3).
  std::array<u64, kAlphabetSize> _alphabet = {0ULL, 0ULL, 0ULL, 0ULL};

  /// @brief Scanning direction: kBackward reverses the query for reverse matching.
  const SearchDirection _direction;

  /// @brief Maximum number of reference start positions searchable per window.
  const s32 _max_search_size;

  /// @brief Bit mask for the most significant bit used in match detection.
  const u64 _msb;

  /**
   * @brief Convert a DNA reference substring to internal 0-3 representation.
   *
   * @param ref Reference sequence (ASCII DNA).
   * @param begin Start position in reference (inclusive).
   * @param end End position in reference (inclusive).
   * @return Array of base indices suitable for the Find overload.
   * @throws error::Error If begin > end or the search length exceeds kQueryWindowSize.
   */
  std::array<u8, kBufferSize> ToBitapAlphabetArray(std::string_view ref, s32 begin, s32 end) const;

  /**
   * @brief Compute the MSB mask for a given query length.
   *
   * @param query_length Length of the query sequence.
   * @return Bitmask with the bit at position (query_length - 1) set.
   * @throws error::Error If query_length exceeds kQueryWindowSize.
   */
  static u64 ComputeMSB(size_t query_length);

  /// @brief Controls whether the internal search loop returns the first or last hit.
  enum class MatchPolicy : u8 {
    kFirst,
    kWidest
  };

  /// @brief Number of edit-distance levels tracked by the search (distances 0..MaxDist).
  /// NOSONAR(M23_051,M23_053) -- single definition (see below usages); the analyzer misreads each explicit
  /// Bitap<N> instantiation as redefining this static constexpr member, the same false positive already
  /// suppressed on the explicit-instantiation lines in bitap.cpp.
  static constexpr s32 kNumLevels = static_cast<s32>(MaxDist) + 1;  // NOSONAR(M23_051,M23_053)

  /**
   * @brief Advance the approximate-match states for edit distances 1..MaxDist at one reference position.
   *
   * Extracted from SearchWindow to keep that kernel's control flow flat. Updates current_state,
   * previous_state, and latches best_position per the match policy.
   *
   * @tparam Policy kFirst latches the earliest hit per level; kWidest keeps the latest.
   * @param current_pattern_mask Pattern mask for the current reference base.
   * @param position Index of the current reference base within the window.
   * @param current_state In/out approximate-match state per level.
   * @param previous_state In/out previous-column state per level.
   * @param best_position In/out best match position per level.
   */
  template <MatchPolicy Policy>
  void UpdateDistanceLevels(u64 current_pattern_mask,
                            s32 position,
                            std::array<u64, kNumLevels>& current_state,
                            std::array<u64, kNumLevels>& previous_state,
                            std::array<s32, kNumLevels>& best_position) const;

  /**
   * @brief Core bitap search over a single window, parameterized by match policy.
   *
   * kFirst returns the earliest match at the best edit distance (with early exit on exact match).
   * kWidest scans the full window and returns the latest (rightmost) match at the best edit distance.
   *
   * @tparam Policy kFirst or kWidest.
   * @param alphabet_ref Pre-converted reference array.
   * @param begin Start position in the original reference (inclusive).
   * @param end End position in the original reference (inclusive).
   * @return Position of match, or kNoMatchPosition if not found.
   */
  template <MatchPolicy Policy>
  s32 SearchWindow(const std::array<u8, kBufferSize>& alphabet_ref, s32 begin, s32 end) const;

  /** @brief Score-aware implementation of SearchWindow. */
  template <MatchPolicy Policy>
  Match SearchWindowMatch(const std::array<u8, kBufferSize>& alphabet_ref, s32 begin, s32 end) const;
};

}  // namespace xoos::demux
