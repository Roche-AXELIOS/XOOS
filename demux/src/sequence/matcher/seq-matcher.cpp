#include "seq-matcher.h"

#include <cassert>
#include <utility>

#include "simd/simd-functions.h"
#include "utility/range-util.h"

namespace xoos::demux {
constexpr uint64_t kLengthMasks64[] = {0ul,
                                       0x3ul,
                                       0xful,
                                       0x3ful,
                                       0xfful,
                                       0x3fful,
                                       0xffful,
                                       0x3ffful,
                                       0xfffful,
                                       0x3fffful,
                                       0xffffful,
                                       0x3ffffful,
                                       0xfffffful,
                                       0x3fffffful,
                                       0xffffffful,
                                       0x3ffffffful,
                                       0xfffffffful,
                                       0x3fffffffful,
                                       0xffffffffful,
                                       0x3ffffffffful,
                                       0xfffffffffful,
                                       0x3fffffffffful,
                                       0xffffffffffful,
                                       0x3ffffffffffful,
                                       0xfffffffffffful,
                                       0x3fffffffffffful,
                                       0xffffffffffffful,
                                       0x3ffffffffffffful,
                                       0xfffffffffffffful,
                                       0x3fffffffffffffful,
                                       0xffffffffffffffful,
                                       0x3ffffffffffffffful,
                                       0xfffffffffffffffful};

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
SeqMatcher::ExcisionLoci SeqMatcher::CreateRelativeExcisionLoci(int gt_seq_len, int max_edist, int max_wiggle_left,
                                                                int max_wiggle_right) {
  auto start_positions = Range(-max_wiggle_left, gt_seq_len + max_wiggle_right + 1, 1);

  auto min_len = std::max(int{1}, gt_seq_len - max_edist);
  auto max_len = gt_seq_len + max_edist;
  auto lengths = Range(min_len, max_len + 1, 1);

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
  auto loci_compare = [gt_seq_len](const Loci& a, const Loci& b) {
    auto a_0 = std::abs(a.spos);
    auto b_0 = std::abs(b.spos);
    if (a_0 == b_0) {
      auto a_1 = std::abs(gt_seq_len - a.length);
      auto b_1 = std::abs(gt_seq_len - b.length);
      return a_1 < b_1;
    } else {
      return a_0 < b_0;
    }
  };

  std::stable_sort(loci.begin(), loci.end(), loci_compare);
  // Analysis of matching results revealed that some positions/edist combinations did not get any hit.
  // Removing them should yield identical results and slightly better performance.
  if (loci.size() == 140) {
    loci.erase(loci.begin() + 34);
    loci.erase(loci.begin() + 28);
    loci.erase(loci.begin() + 25);
    loci.erase(loci.begin() + 18);
    loci.erase(loci.begin() + 15);
    loci.erase(loci.begin() + 12);
    loci.erase(loci.begin() + 10);
    loci.erase(loci.begin() + 6);
    loci.erase(loci.begin() + 4);
    loci.erase(loci.begin() + 2);
  }
  /*
      2	        0	15
      4	        0	16
      6	        1	15
      10	1	16
      12	-1	15
      15	2	16
      18	2	15
      25	3	16
      28	3	15
      34	4	16
    */
  // If the prefilter does not have a hit, we can potentially skip a few positions - calculate the new index
  // position that should be used if that occurs.
  for (auto i1 = 0ul; i1 < loci.size();) {
    const auto& loci1{loci[i1]};
    auto i2 = i1 + 1ul;
    // loop over i2 until you find different start position
    for (; i2 < loci.size() && loci[i2].spos == loci1.spos; ++i2) {
    }  // NOLINT
    for (; i1 < i2; ++i1) {
      loci[i1].skip = static_cast<int>(i2 - i1);
    }
  }
  return loci;
}

SeqMatcher::SeqMatcher(uint seq_len, int max_edist, int max_wiggle_left, int max_wiggle_right, SeqLutPtr lut)
    : _seq_len{seq_len},
      _relative_excision_loci{
          CreateRelativeExcisionLoci(static_cast<int>(seq_len), max_edist, max_wiggle_left, max_wiggle_right)},
      _lut{std::move(lut)},
      _nr_loci(_relative_excision_loci.size()) {}

// Helper structure to hold immediate results
struct MatchValue {
  uint64_t hash64;     // 64-bit hash value to be used for unordered_map
  uint32_t hash29;     // 29-bit hash value to be used for fast check
  uint32_t last_bits;  // last few bits from checksum - if non-zero, match is found
};

constexpr uint32_t kSetBits[] = {1, 2, 4, 8, 16, 32, 64, 128};

static inline bool Prefilter(MatchValue& match_value, const uint8_t* two_bit, const Loci& loci, uint start,
                             const uint8_t* const p_predata, int prefilter_mask) {
  uint startpos{loci.spos + start};
  uint start2{startpos >> 2};                    // start offset in 2-bit representation
  auto hash64 = simd::Load64(two_bit + start2);  // Load 64-bit worth of data.
  // Start = 0: lsb of hash at position 0, OK
  // Start = 1: lsb of hash at position 2, need to right shift in that data
  // etc.
  hash64 >>= (startpos + startpos) & 7;
  // The lower 32 bits now contain the data we're after. Strip out the bits
  // we do not need.
  hash64 &= loci.mask64;

  // Switch to 32-bit integers, which is sufficient for barcodes.
  auto hash32 = static_cast<uint32_t>(hash64);
  // The last three bits get special treatment, we'll assign one bit for every possible combination of 3 bits
  auto set_bits{kSetBits[hash32 & 7]};
  // Switch to 29-bit representation, means that the largest LUT (16 bases) will have 0.5 GB size
  match_value.hash29 = hash32 >> 3;
  // Use a subset of the 29 bits for a pre-filter with a CPU-cache friendly lookup to discard trivial cases
  match_value.last_bits = p_predata[match_value.hash29 & prefilter_mask] & set_bits;

  if (match_value.last_bits) {
    match_value.hash64 = hash64;
    return true;
  }
  return false;
}

MatchInfo SeqMatcher::FindBarcode(ReadEnd read_end, uint start_pos, const uint8_t* two_bit, size_t seq_length) const {
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

  const auto startpos{static_cast<int>(start_pos)};
  const auto maxpos{static_cast<int>(seq_length) - startpos};
  const auto minpos{-startpos};
  const auto prefilter_mask{_lut->PrefilterMask()};

  MatchInfo match_info;  // result of match; is initialized to non-match and updated in the inner loop
  constexpr auto kMaxSize{200};
  assert(kMaxSize >= _nr_loci);

  MatchValue match_values[kMaxSize];  // scratch information used by inner loop
  const auto* const p_cache{&_lut->PrefilterValues()[0]};

  size_t prefilter{0};  // index of the prefilter loop
  int nr_active{0};     // how many candidates are under evaluation

  // loop over all combinations of positions/length.
  for (size_t index = 0; index < _nr_loci; ++index) {
    if (prefilter < _nr_loci) {
      // Apply the prefilter; do that until we found a hit or run out of positions
      do {
        // get start and end position for this candidate, test whether it falls inside the sequence
        const auto& loci = _relative_excision_loci[prefilter];
        if (loci.spos >= minpos && loci.epos <= maxpos) {
          // Inside the sequence. Now prefilter; if the prefilter has no hit, it means that any length of
          // the fragment will not have a hit. We usually can safely skip a few samples now; information
          // about that was precomputed and stored in the 'skip' array.
          if (!Prefilter(match_values[prefilter], two_bit, loci, startpos, p_cache, prefilter_mask)) {
            // prefilter did not have any hit, so skip a few samples. We skip by setting last bits status
            // to zero. Used Duff's device to unroll the code that would otherwise be needed.
            switch (loci.skip) {  // Duffs device to unroll the loop where we mark the positions as no hit.
              case 5:             // NOLINT
                match_values[++prefilter].last_bits = 0;
              case 4:
                match_values[++prefilter].last_bits = 0;
              case 3:
                match_values[++prefilter].last_bits = 0;
              case 2:
                match_values[++prefilter].last_bits = 0;
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
            ++prefilter;      // required because we jumped out of the loop
            goto do_process;  // Yuck, but speeds things up a tad.
          }
        } else {
          // Barcode would fall outside the range, so mark as irrelevant. This happens quite frequent, I discovered.
          match_values[prefilter].last_bits = 0;
          // Outside range, nothing found yet - we do not have to spend any time on further checking, so
          // skip that by setting the index variable to the current prefilter index - don't bother to flag
          // the bits values either as they will not be tested anymore.
          if (!nr_active) {
            index = prefilter;
          }
        }
      } while (++prefilter < _nr_loci && !nr_active);
      // If we exited the loop because we prefiltered all positions, we can terminate immediately if we still
      // did not find a hit.
      if (!nr_active) {
        break;
      }
    }

  do_process:
    // Second half of the loop checks whether prefilter found a candidate; if so, it checks the binary filter whether
    // it is a match and if so, uses the more expensive hash table to find the barcode information.
    const auto& mv{match_values[index]};
    if (mv.last_bits) {
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
      --nr_active;  // we processed a candidate, decrease count.
    }
  }

  return match_info;
}

SequenceTwoBit::SequenceTwoBit(const std::string_view& seq) : _length(seq.length()) {  // NOLINT
  if ((_length >> 2) > kMaxMemorySequence) {
    // TODO: throw exception here
  }
  simd::ConvertTo2Bit(reinterpret_cast<const uint8_t*>(seq.data()), 0, _length, _two_bit_data + kOffset);
}

const BarcodePool& SeqMatcher::Pool() const { return _lut->Pool(); }

}  // namespace xoos::demux
