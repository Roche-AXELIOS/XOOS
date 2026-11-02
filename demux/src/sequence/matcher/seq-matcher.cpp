#include "seq-matcher.h"

#include <algorithm>
#include <array>
#include <cassert>

#include <xoos/types/vec.h>
#include <xoos/util/container-functions.h>

#include "sequence/alignment/log-likelihood-scoring.h"
#include "sequence/sequence-enums.h"
#include "simd/simd-functions.h"
#include "utility/range-util.h"

namespace xoos::demux {

constexpr std::array kLengthMasks64 = {0ULL,
                                       0x3ULL,
                                       0xfULL,
                                       0x3fULL,
                                       0xffULL,
                                       0x3ffULL,
                                       0xfffULL,
                                       0x3fffULL,
                                       0xffffULL,
                                       0x3'ffffULL,
                                       0xf'ffffULL,
                                       0x3f'ffffULL,
                                       0xff'ffffULL,
                                       0x3ff'ffffULL,
                                       0xfff'ffffULL,
                                       0x3fff'ffffULL,
                                       0xffff'ffffULL,
                                       0x3'ffff'ffffULL,
                                       0xf'ffff'ffffULL,
                                       0x3f'ffff'ffffULL,
                                       0xff'ffff'ffffULL,
                                       0x3ff'ffff'ffffULL,
                                       0xfff'ffff'ffffULL,
                                       0x3fff'ffff'ffffULL,
                                       0xffff'ffff'ffffULL,
                                       0x3'ffff'ffff'ffffULL,
                                       0xf'ffff'ffff'ffffULL,
                                       0x3f'ffff'ffff'ffffULL,
                                       0xff'ffff'ffff'ffffULL,
                                       0x3ff'ffff'ffff'ffffULL,
                                       0xfff'ffff'ffff'ffffULL,
                                       0x3fff'ffff'ffff'ffffULL,
                                       0xffff'ffff'ffff'ffffULL};

constexpr size_t kLociListSizeForDeadPositionsErase = 140;
constexpr std::array<size_t, 10> kDeadLociPositions = {2, 4, 6, 10, 12, 15, 18, 25, 28, 34};

/**
 * Excision loci are the potential (start, end) positions for finding a barcode within a read.
 * It generates all (start, end) within [-wiggle_left, wiggle_right + barcode_length] such that
 * (end - start + 1) should not be greater than barcode_length.
 *
 * @param gt_seq_len - Barcode sequence length
 * @param max_edist - Maximum edit distance allowed while matching
 * @param max_wiggle_left - Max wiggle Left
 * @param max_wiggle_right - Max wiggle right
 *
 * @return vector of tuple of potential loci
 */
SeqMatcher::ExcisionLoci SeqMatcher::CreateRelativeExcisionLoci(const s32 gt_seq_len,
                                                                const s32 max_edist,
                                                                const s32 max_wiggle_left,
                                                                const s32 max_wiggle_right) {
  const auto start_positions = Range(-max_wiggle_left, gt_seq_len + max_wiggle_right + 1, 1);

  const auto min_len = std::max(1, gt_seq_len - max_edist);
  const auto max_len = gt_seq_len + max_edist;
  const auto lengths = Range(min_len, max_len + 1, 1);

  auto loci = std::vector<Loci>{};
  for (const auto spos : start_positions) {
    for (const auto len : lengths) {
      auto epos = spos + len;
      if (epos > start_positions.back()) {
        continue;
      }
      if (std::abs(gt_seq_len - (epos - spos)) > max_edist) {
        continue;
      }
      loci.emplace_back(kLengthMasks64[epos - spos], spos, epos, epos - spos);
    }
  }

  // Sort by closest to the relative-zero position and then by
  // distance to the ground truth (expected) sequence length
  auto loci_compare = [&gt_seq_len = std::as_const(gt_seq_len)](const Loci& a, const Loci& b) {
    const auto a_0 = std::abs(a.spos);
    const auto b_0 = std::abs(b.spos);
    if (a_0 == b_0) {
      const auto a_1 = std::abs(gt_seq_len - a.length);
      const auto b_1 = std::abs(gt_seq_len - b.length);
      return a_1 < b_1;
    }
    return a_0 < b_0;
  };

  std::ranges::stable_sort(loci, loci_compare);
  // Analysis of matching results revealed that some positions/edist combinations did not get any hit.
  // Removing them should yield identical results and slightly better performance.
  if (loci.size() == kLociListSizeForDeadPositionsErase) {
    util::container::VectorErase(loci, vec<size_t>(kDeadLociPositions.begin(), kDeadLociPositions.end()));
  }
  /*
      2            0    15
      4            0    16
      6            1    15
      10    1    16
      12    -1    15
      15    2    16
      18    2    15
      25    3    16
      28    3    15
      34    4    16
    */
  // If the prefilter does not have a hit, we can potentially skip a few positions - calculate the new index
  // position that should be used if that occurs.
  {
    auto i1 = 0UL;
    while (i1 < loci.size()) {
      const auto& loci1{loci[i1]};
      auto i2 = i1 + 1UL;
      // loop over i2 until you find different start position
      while (i2 < loci.size() && loci[i2].spos == loci1.spos) {
        ++i2;
      }
      while (i1 < i2) {
        loci[i1].skip = static_cast<s32>(i2 - i1);
        ++i1;
      }
    }
  }
  return loci;
}

SeqMatcher::SeqMatcher(
    const u32 seq_len, const s32 max_edist, const s32 max_wiggle_left, const s32 max_wiggle_right, SeqLutPtr lut)
    : _seq_len{seq_len},
      _max_edist{static_cast<u32>(max_edist)},
      _relative_excision_loci{
          CreateRelativeExcisionLoci(static_cast<s32>(seq_len), max_edist, max_wiggle_left, max_wiggle_right)},
      _lut{std::move(lut)},
      _nr_loci(_relative_excision_loci.size()) {
}

// Helper structure to hold immediate results
struct MatchValue {
  // 64-bit hash value to be used for unordered_map
  u64 hash64;
  // 29-bit hash value to be used for fast check
  u32 hash29;
  // last few bits from checksum - if non-zero, match is found
  u32 last_bits;
};

constexpr std::array<u32, 8> kSetBits = {1, 2, 4, 8, 16, 32, 64, 128};
// Mask to extract the last 3 bits
constexpr u32 kLast3BitsMask = 0x7;
// Number of bits used for indexing into the prefilter byte (3 bits for 8 combinations)
constexpr auto kNumBitsForBitPosInPrefilterByte = 3;

static inline bool Prefilter(MatchValue& match_value,
                             const u8* const two_bit,
                             const Loci& loci,
                             const u32 start,
                             const u8* const p_predata,
                             const s32 prefilter_mask) {
  const u32 start_pos{static_cast<u32>(loci.spos) + start};
  // start offset in 2-bit representation
  const u32 start2{start_pos >> 2};
  auto hash64 = simd::Load64(two_bit + start2);  // Load 64-bit worth of data.
  // Start = 0: lsb of hash at position 0, OK
  // Start = 1: lsb of hash at position 2, need to right shift in that data
  // etc.
  hash64 >>= (start_pos + start_pos) & kLast3BitsMask;
  // The lower 32 bits now contain the data we're after. Strip out the bits
  // we do not need.
  hash64 &= loci.mask64;

  // Switch to 32-bit integers, which is sufficient for barcodes.
  const auto hash32 = static_cast<u32>(hash64);
  // The last three bits get special treatment, we'll assign one bit for every possible combination of 3 bits
  const auto set_bits{kSetBits[hash32 & kLast3BitsMask]};
  // Switch to 29-bit representation, means that the largest LUT (16 bases) will have 0.5 GB size
  match_value.hash29 = hash32 >> kNumBitsForBitPosInPrefilterByte;
  // Use a subset of the 29 bits for a pre-filter with a CPU-cache friendly lookup to discard trivial cases
  match_value.last_bits = static_cast<u32>(p_predata[match_value.hash29 & static_cast<u32>(prefilter_mask)]) & set_bits;

  if (match_value.last_bits != 0U) {
    match_value.hash64 = hash64;
    return true;
  }
  return false;
}

MatchInfo SeqMatcher::FindBarcode(const ReadEnd read_end,
                                  u32 start_pos,
                                  const u8* const two_bit,
                                  const size_t seq_length) const {
  if (read_end == ReadEnd::k3p) {
    // If we are looking for a barcode from the 3' direction (right to left), then start_pos is actually
    // the end position, because of this we subtract _seq_len to determine the real start position
    start_pos = _seq_len >= start_pos ? 0 : start_pos - _seq_len;
  }

  // Because we do not know the exact start and end position of the barcode sequence, we will check
  // many different potential start and positions until we find an exact match, or we will aggregate
  // partial or ambiguous matches until we exhaust all attempts.
  // To eliminate overhead associated with allocation of vector to hold absolute positions, we're now
  // doing the legwork to determine the absolute positions from relative positions within this loop.

  const auto startpos{static_cast<s32>(start_pos)};
  const auto maxpos{static_cast<s32>(seq_length) - startpos};
  const auto minpos{-startpos};
  const auto prefilter_mask{_lut->PrefilterMask()};

  // result of match; is initialized to non-match and updated in the inner loop
  MatchInfo match_info;
  constexpr auto kMaxSize{200};
  assert(kMaxSize >= _nr_loci);

  // scratch information used by inner loop
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
  std::array<MatchValue, kMaxSize> match_values;
  const auto* const p_cache{_lut->PrefilterValues().data()};

  // index of the prefilter loop
  size_t prefilter_index{0};
  // how many candidates are under evaluation
  s32 nr_active{0};

  // loop over all combinations of positions/length.
  for (size_t index = 0; index < _nr_loci; ++index) {
    if (prefilter_index < _nr_loci) {
      // Apply the prefilter; do that until we found a hit or run out of positions
      do {
        // get start and end position for this candidate, test whether it falls inside the sequence
        if (const auto& loci = _relative_excision_loci[prefilter_index]; loci.spos >= minpos && loci.epos <= maxpos) {
          // Inside the sequence. Now prefilter; if the prefilter has no hit, it means that any length of
          // the fragment will not have a hit. We usually can safely skip a few samples now; information
          // about that was precomputed and stored in the 'skip' array.
          if (!Prefilter(
                  match_values[prefilter_index], two_bit, loci, static_cast<u32>(startpos), p_cache, prefilter_mask)) {
            // prefilter did not have any hit, so skip a few samples. We skip by setting last bits status
            // to zero. Used Duff's device to unroll the code that would otherwise be needed.
            switch (loci.skip) {  // Duffs device to unroll the loop where we mark the positions as no hit.
              case 5:             // NOLINT
                ++prefilter_index;
                match_values[prefilter_index].last_bits = 0;
              case 4:
                ++prefilter_index;
                match_values[prefilter_index].last_bits = 0;
              case 3:
                ++prefilter_index;
                match_values[prefilter_index].last_bits = 0;
              case 2:
                ++prefilter_index;
                match_values[prefilter_index].last_bits = 0;
              default:
                break;
            }
          } else {
            // We encountered a hit, we need to evaluate further. The required information was already
            // marked by the Prefilter() function and data required to make the final decision is "underway"
            // as we prefetched it.
            ++nr_active;
            // Yeah, I hate gotos too - in fact, never used them until now. But: if we did encountered a hit,
            // we should not waste time on wrapping up the prefilter loop, but immediately proceed to process it
            // ASAP.
            // required because we jumped out of the loop
            ++prefilter_index;
            // Yuck, but speeds things up a tad.
            // TODO: refactor the code to eliminate the need for this
            goto do_process;
          }
        } else {
          // Barcode would fall outside the range, so mark as irrelevant. This happens quite frequent, I discovered.
          match_values[prefilter_index].last_bits = 0;
          // Outside range, nothing found yet - we do not have to spend any time on further checking, so
          // skip that by setting the index variable to the current prefilter index - don't bother to flag
          // the bits values either as they will not be tested anymore.
          if (nr_active == 0) {
            index = prefilter_index;
          }
        }
        ++prefilter_index;
      } while (prefilter_index < _nr_loci && nr_active == 0);
      // If we exited the loop because we prefiltered all positions, we can terminate immediately if we still
      // did not find a hit.
      if (nr_active == 0) {
        break;
      }
    }

  do_process:
    // Second half of the loop checks whether prefilter found a candidate; if so, it checks the binary filter whether
    // it is a match and if so, uses the more expensive hash table to find the barcode information.
    const auto& mv{match_values[index]};
    if (mv.last_bits != 0U) {
      // This entry was found by the prefilter, so we need to process it.
      const auto& loci = _relative_excision_loci[index];
      auto [match, match_type] = _lut->SimpleFind(loci.length, mv.hash64);
      if (match_type != MatchType::kUnknown) {
        // Found a hash code, update our status
        match_info.Update(match, match_type, loci.spos + startpos, loci.epos + startpos);
        if (match_type == MatchType::kExact) {
          // stop if we found a complete match
          break;
        }
      }
      // we processed a candidate, decrease count.
      --nr_active;
    }
  }

  return match_info;
}

void SeqMatcher::FindNextBarcode(const ReadEnd read_end,
                                 const s32 pos,
                                 const u8* const two_bit,
                                 const size_t seq_length,
                                 std::vector<MatchInfo>& results) const {
  const s32 min_len = std::max(1, static_cast<s32>(_seq_len) - static_cast<s32>(_max_edist));
  const auto max_len = static_cast<s32>(_seq_len + _max_edist);
  const auto slen = static_cast<s32>(seq_length);
  const auto prefilter_mask = _lut->PrefilterMask();
  const auto* const p_cache = _lut->PrefilterValues().data();

  results.clear();

  // Lambda: at a given absolute start position, try all candidate lengths and report the single
  // best match. Preference: exact > approximate; among approximate, prefer the length closest
  // to nominal (_seq_len) which yields the highest LogLikelihoodScore.
  const auto nominal = static_cast<s32>(_seq_len);
  auto try_at_position = [this, &min_len, &max_len, &nominal, &two_bit, &p_cache, &prefilter_mask, &results](
                             const s32 start, const s32 max_end) {
    // Track the best candidate at this start
    bool found = false;
    BarcodeMatch best_match{};
    MatchType best_type = MatchType::kUnknown;
    s32 best_len = 0;
    s32 best_score = std::numeric_limits<s32>::min();

    for (s32 len = min_len; len <= max_len && start + len <= max_end; ++len) {
      // Compute 2-bit hash at this (start, length)
      const auto abs_start = static_cast<u32>(start);
      const u32 start2 = abs_start >> 2;
      auto hash64 = simd::Load64(two_bit + start2);
      hash64 >>= (abs_start * 2) & kLast3BitsMask;
      hash64 &= kLengthMasks64[static_cast<u32>(len)];

      // Prefilter check
      const auto hash32 = static_cast<u32>(hash64);
      const auto set_bits = kSetBits[hash32 & kLast3BitsMask];
      // TODO: switch to std::byte for 8 bit manipulation instead of char
      if (const auto hash29 = hash32 >> kNumBitsForBitPosInPrefilterByte;
          (p_cache[hash29 & static_cast<u32>(prefilter_mask)] & set_bits) == 0) {
        continue;
      }

      // Prefilter passed — full LUT lookup
      auto [match, match_type] = _lut->SimpleFind(static_cast<u32>(len), hash64);
      if (match_type == MatchType::kUnknown) {
        continue;
      }

      // Exact match is the best possible — emit immediately
      if (match_type == MatchType::kExact) {
        results.emplace_back(match, LociRange{static_cast<u32>(start), static_cast<u32>(start + len)}, match_type);
        return true;
      }

      // Among approximate matches, rank by LogLikelihoodScore (accounts for both edist and length)
      const s32 net_indels = len - nominal;
      const s32 abs_net = net_indels < 0 ? -net_indels : net_indels;
      const s32 num_insertions = net_indels > 0 ? net_indels : 0;
      const s32 num_deletions = net_indels < 0 ? -net_indels : 0;
      const s32 num_substitutions = static_cast<s32>(match.edist) - abs_net;
      const s32 num_matches = nominal - num_substitutions - num_deletions;
      const s32 score = (num_matches * scoring::kMatch) + (num_substitutions * scoring::kSubstitution) +
                        (num_insertions * scoring::kInsertion) + (num_deletions * scoring::kDeletion);

      if (!found || score > best_score) {
        best_match = match;
        best_type = match_type;
        best_len = len;
        best_score = score;
        found = true;
      }
    }

    if (found) {
      results.emplace_back(
          best_match, LociRange{static_cast<u32>(start), static_cast<u32>(start + best_len)}, best_type);
    }
    return found;
  };

  if (read_end == ReadEnd::k5p) {
    // Scan left-to-right from pos
    for (s32 start = pos; start + min_len <= slen; ++start) {
      // greedy: take first position with a match
      if (try_at_position(start, slen)) {
        break;
      }
    }
  } else {
    // Scan right-to-left: find the rightmost barcode ending at or before pos.
    // Iterate start positions from high to low so we find the closest match first.
    const s32 highest_start = std::min(pos - min_len, slen - min_len);
    // After finding a first match, continue scanning lower starts to find potentially
    // longer (better) matches ending at the same EPos. The max start difference for
    // same-EPos matches is (max_len - min_len).
    const s32 extra_scan = max_len - min_len;
    bool found_any = false;
    s32 first_match_start = 0;

    for (s32 start = highest_start; start >= 0; --start) {
      if (found_any && start < first_match_start - extra_scan) {
        // We've scanned far enough past the first match — any match here would have
        // a lower EPos, so stop.
        break;
      }
      if (try_at_position(start, pos) && !found_any) {
        found_any = true;
        first_match_start = start;
      }
    }
  }
}

SequenceTwoBit::SequenceTwoBit(const std::string_view& seq) : _length(seq.length()) {  // NOLINT
  // TODO: evaluate whether to include offest in the length memmory check
  if ((_length >> 2) > kMaxMemorySequence) {
    // TODO: throw exception here
  }
  simd::ConvertTo2Bit(reinterpret_cast<const u8*>(seq.data()), 0, static_cast<u32>(_length), _two_bit_data + kOffset);
}

const BarcodePool& SeqMatcher::Pool() const {
  return _lut->Pool();
}

}  // namespace xoos::demux
