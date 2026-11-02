#include "trim-duplex-umi.h"

#include <array>
#include <cstdlib>
#include <numeric>
#include <utility>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>

#include "adapters/duplex/duplex-match.h"
#include "adapters/duplex/stem/trim-duplex-stem.h"
#include "adapters/duplex/trim-info-duplex.h"
#include "io/read-record.h"
#include "sequence/encoding-constants.h"
#include "sequence/matcher/match-info.h"
#include "sequence/matcher/match-position-states.h"
#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

// Should make dynamic based on the LUT used in the search which varies
// The k-mer size of the LUT, we are hardcoding to 10 because the UMI length is 12 with 2 indels allowed
constexpr s32 kUmiLutKmerLength = 10;
constexpr s64 kFivePrimeSearchOffsetStrideBases = 4;
constexpr s64 kFivePrimeExtraSearchOffset0 = 0;
constexpr s64 kFivePrimeExtraSearchOffset1 = kFivePrimeSearchOffsetStrideBases;
constexpr s64 kFivePrimeExtraSearchOffset2 = kFivePrimeSearchOffsetStrideBases + kFivePrimeSearchOffsetStrideBases;
constexpr s32 kThreePrimeTieBreakPenalty = 1;

// Positional constructor; adjacent FixedSequence params encode fixed read-structure positions
// documented in the class header and are constructed positionally by the assembler.
TrimDuplexUmi::TrimDuplexUmi(FixedSequence runway,  // NOSONAR(S5419)
                             FixedSequence start,
                             const FixedSequence& stem,
                             SeqMatcher umi_5p_matcher,
                             SeqMatcher umi_3p_matcher)
    : TrimDuplexStem(std::move(runway), std::move(start), stem),
      _umi_5p_matcher(std::move(umi_5p_matcher)),
      _umi_3p_matcher(std::move(umi_3p_matcher)),
      _umi_lut(CascadedLUTs(
          _umi_5p_matcher, _umi_3p_matcher, DuplexMatch::BarcodeType::kUMI5p, DuplexMatch::BarcodeType::kUMI3p)),
      _mask_umi5p(static_cast<s64>((1UL << (kBitsPerEncodedBase * _umi_lut.MaxLength(DuplexMatch::kUMI5p))) - 1UL)),
      _mask_umi3p(static_cast<s64>((1UL << (kBitsPerEncodedBase * _umi_lut.MaxLength(DuplexMatch::kUMI3p))) - 1UL)),
      // Array of masks for LUT search; the last four lanes are unused (mask 0).
      _lut_offset_umi_masks{_mask_umi5p, _mask_umi5p, _mask_umi3p, _mask_umi3p, 0, 0, 0, 0} {
  // Check if length of components are compatible with our length hack
  if (static_cast<s32>(Start().length() + RunwaySequence().length()) < _umi_info.umi_5p_offset_4bp) {
    throw error::Error("Start + Runway length is too short for UMI 5' 4bp offset.");
  }
  if (static_cast<s32>(Start().length() + RunwaySequence().length()) < _umi_info.umi_5p_offset_6bp) {
    throw error::Error("Start + Runway length is too short for UMI 5' 6bp offset.");
  }
}

/**
 * @brief Sorts the UMI results in the FixedReadRecord.
 * This function sorts the unfiltered matches by edit distance and then filters them into buckets based on their type.
 * Filters out matches that are out of bounds of the read length
 * @param record a FixedReadRecord containing the unfiltered matches to sort
 */
static void SortUMIs(FixedReadRecord& record, const s64 expect_pos_5p, const s64 expect_pos_3p) {
  // std::midpoint rounds toward its first argument; pass the smaller position first so the result floors.
  const auto average_search_pos =
      std::midpoint(std::min(expect_pos_5p, expect_pos_3p), std::max(expect_pos_5p, expect_pos_3p));

  auto& trim_info = record.trim_info_duplex;
  const auto nr_results = trim_info.unfiltered_matches_count;
  if (nr_results == 0) {
    return;
  }

  auto& results{trim_info.unfiltered_matches};
  // First step: sort the markers in increasing edit distance, this will put the most reliable markers first.
  std::sort(results.begin(), results.begin() + nr_results);
  const s64 max_valid_pos = record.consensus_seq_len - kUmiLutKmerLength;
  u32 count_5p = 0;
  u32 count_3p = 0;
  for (u32 i = 0; i < nr_results; ++i) {
    // Only include matches that are within the read boundaries
    if (std::cmp_greater(results[i].pos, max_valid_pos)) {
      continue;
    }
    const auto barcode_type = static_cast<DuplexMatch::BarcodeType>(results[i].match.barcode_type);
    if (barcode_type == DuplexMatch::kUMI5p) {
      // filter out matches that are beyond the average search position
      if (std::cmp_greater(results[i].pos, average_search_pos)) {
        continue;
      }
      trim_info.filtered_matches[DuplexMatch::kUMI5p][count_5p] = results[i];
      ++count_5p;
    } else if (barcode_type == DuplexMatch::kUMI3p) {
      // filter out matches that are before the average search position
      if (std::cmp_less(results[i].pos, average_search_pos)) {
        continue;
      }
      trim_info.filtered_matches[DuplexMatch::kUMI3p][count_3p] = results[i];
      ++count_3p;
    } else {
      throw error::Error("Unknown barcode type: {} while UMI trimming.",
                         static_cast<s32>(results[i].match.barcode_type));
    }
  }
  trim_info.filtered_matches_count[DuplexMatch::kUMI5p] = count_5p;
  trim_info.filtered_matches_count[DuplexMatch::kUMI3p] = count_3p;
}

/**
 * @brief Finds the 5' and 3' UMI in the read.
 * Uses Bitap to find the position of the 5' and 3' UMI, then uses a cascading LUT to find the actual UMI sequence.
 * @param record a FixedReadRecord containing the sequence to search
 * @return a pair of integers representing the positions of the 5' and 3' UMI to trim, and -1 if not found.
 */
namespace {
// Position adjustment (in bases) implied by the encoded indel length of a UMI match.
s32 IndelPositionDelta(const u16 encoded_length) {
  switch (encoded_length) {
    case DuplexMatch::kPlus2:
      return DuplexMatch::kPos2;
    case DuplexMatch::kPlus1:
      return DuplexMatch::kPos1;
    case DuplexMatch::kZero:
      return 0;
    case DuplexMatch::kMinus1:
      return DuplexMatch::kNeg1;
    case DuplexMatch::kMinus2:
      return DuplexMatch::kNeg2;
    default:
      throw error::Error("Unknown UMI indel size: {}", static_cast<s32>(encoded_length));
  }
}

// Among the filtered matches for one UMI end, pick the one whose implied trim position is closest to the
// bitap-anchored expectation, breaking ties in favour of a lower edit distance.
DuplexResult SelectBestUmiMatch(const std::array<DuplexResult, TrimInfoDuplex::kMaxFilteredMatches>& matches,
                                const u32 count,
                                const s32 expected_pos,
                                const BarcodePool& pool,
                                const bool use_sequence_offset,
                                const s32 tie_break_penalty) {
  auto best = matches[0];
  if (expected_pos == kNoMatchPosition) {
    return best;
  }
  const auto offset_for = [&pool, use_sequence_offset](const DuplexResult& r) {
    const auto& barcode = pool[r.match.barcode_id];
    return use_sequence_offset ? static_cast<s32>(barcode.sequence.length() - barcode.name.length() - 1)
                               : static_cast<s32>(barcode.name.length());
  };
  auto diff_best = std::abs(best.pos + offset_for(best) - expected_pos);
  for (u32 i = 1; i < count; ++i) {
    const auto& match = matches[i];
    const auto diff_current = std::abs(match.pos + offset_for(match) - expected_pos);
    // Prefer a closer match, but only when its net (distance + edit) cost does not increase.
    if (diff_current < diff_best &&
        diff_current + match.match.edist + tie_break_penalty < diff_best + best.match.edist) {
      best = match;
      diff_best = diff_current;
    }
  }
  return best;
}
}  // namespace

std::array<int64_t, TrimDuplexUmi::kSimdSearchLanes> TrimDuplexUmi::ComputeSearchOffsets(const FixedReadRecord& record,
                                                                                         const s32 pos_5p,
                                                                                         const s32 pos_3p) const {
  const s64 assumed_umi_len_3p_4bp = std::max(record.consensus_seq_len - static_cast<s32>(Stem().Size()) - kUmiLen1, 0);
  const s64 assumed_umi_len_3p_6bp = std::max(record.consensus_seq_len - static_cast<s32>(Stem().Size()) - kUmiLen2, 0);

  // we assign the search position to the found location, otherwise we try to search based on theoretical location
  const s64 search_pos_5p_4bp =
      pos_5p == kNoMatchPosition ? _assumed_umi_len_5p_4bp : std::max(0, pos_5p - _umi_info.umi_5p_offset_4bp);
  const s64 search_pos_5p_6bp =
      pos_5p == kNoMatchPosition ? _assumed_umi_len_5p_6bp : std::max(0, pos_5p - _umi_info.umi_5p_offset_6bp);
  const s64 search_pos_3p_4bp = pos_3p == kNoMatchPosition ? assumed_umi_len_3p_4bp : std::max(0, pos_3p - kUmiLen1);
  const s64 search_pos_3p_6bp = pos_3p == kNoMatchPosition ? assumed_umi_len_3p_6bp : std::max(0, pos_3p - kUmiLen2);

  // using SIMD search, it can search 8 positions at once, it searches every 8 bits (4 bases) the offset points to
  // We redundantly search the 5' 0, 4, 8, 12 positions more to capture cases where the 5' stem is not
  return {search_pos_5p_4bp,
          search_pos_5p_6bp,
          search_pos_3p_4bp,
          search_pos_3p_6bp,
          kFivePrimeExtraSearchOffset0,
          kFivePrimeExtraSearchOffset1,
          kFivePrimeExtraSearchOffset2,
          search_pos_3p_6bp};
}

std::pair<s32, s32> TrimDuplexUmi::FindUMI(FixedReadRecord& record) const {
  // Find the 5' UMI
  // Find adapter sequence to anchor search using bitap
  // Sometimes the position ends early due to mismatches so this can still let some of the UMI into the sequence
  // This is  fixed if we detect a valid 5' umi match bit it breaks if too far to search or if the umi is not found
  // NOTE: Ideally need to 1) get info if bitap match was missing trailing bases and 2) adjust the position
  const auto seq = std::string_view(record.consensus_buffer, ToUnsigned(record.consensus_seq_len));
  // Find adapter sequence to anchor search using bitap
  // search using the post-loop stem since it is more likely to exist
  s32 pos_3p = StemBitap().ReverseScan(seq, kNoMatchPosition, kNoMatchPosition);
  // if we found the post-loop stem, search for the pre-loop stem upstream
  // right to left to compensate for artefacts (repeated stems)
  s32 pos_5p = StartBitap().ReverseScan(seq, kNoMatchPosition, pos_3p);

  auto& trim_info = record.trim_info_duplex;

  // reset the UMI info in the record
  trim_info.umi_5p.reset();
  trim_info.umi_3p.reset();

  // reset the matches because we would have existing matches from the hairpin search
  trim_info.unfiltered_matches_count = 0;
  trim_info.filtered_matches_count[DuplexMatch::kUMI3p] = 0;
  trim_info.filtered_matches_count[DuplexMatch::kUMI5p] = 0;

  const std::array<int64_t, kSimdSearchLanes> offsets = ComputeSearchOffsets(record, pos_5p, pos_3p);

  // Convert the consensus sequence to 2-bit representation and store in the record.
  // NOSONAR(S5356) -- char and u8 are both byte types, so this pointer reinterpret is well-defined aliasing.
  record.InitTwoBit(reinterpret_cast<const u8*>(record.consensus_buffer),  // NOSONAR
                    static_cast<u32>(record.consensus_seq_len));

  FindMarker(_umi_lut, record, offsets.data(), _lut_offset_umi_masks.data(), _lut_offset_umi_types.data());
  // The 5' pair holds the 4bp/6bp search positions (ascending); the 3' pair holds them descending.
  static constexpr size_t kOffset5p4bp{0};
  static constexpr size_t kOffset5p6bp{1};
  static constexpr size_t kOffset3p4bp{2};
  static constexpr size_t kOffset3p6bp{3};
  // std::midpoint rounds toward its first argument, so pass the smaller offset first to floor the average.
  // The 3' pair is descending, so the 6bp offset is passed first.
  SortUMIs(record,
           std::midpoint(offsets[kOffset5p4bp], offsets[kOffset5p6bp]),
           std::midpoint(offsets[kOffset3p6bp], offsets[kOffset3p4bp]));
  if (trim_info.filtered_matches_count[DuplexMatch::kUMI5p] > 0) {
    // Among the filtered 5' matches, pick the one closest to the bitap-anchored position.
    const auto best_match_5p = SelectBestUmiMatch(trim_info.filtered_matches[DuplexMatch::kUMI5p],
                                                  trim_info.filtered_matches_count[DuplexMatch::kUMI5p],
                                                  pos_5p,
                                                  _umi_5p_matcher.Pool(),
                                                  /*use_sequence_offset=*/true,
                                                  /*tie_break_penalty=*/0);
    trim_info.umi_5p = best_match_5p.match.barcode_id;
    // adjust the position based on the length variable (actually refers to indel)
    pos_5p = best_match_5p.pos + _umi_info.length_5p - 1 + IndelPositionDelta(best_match_5p.match.length);
  }
  if (trim_info.filtered_matches_count[DuplexMatch::kUMI3p] > 0) {
    // Among the filtered 3' matches, pick the one closest to the bitap-anchored position.
    const auto best_match_3p = SelectBestUmiMatch(trim_info.filtered_matches[DuplexMatch::kUMI3p],
                                                  trim_info.filtered_matches_count[DuplexMatch::kUMI3p],
                                                  pos_3p,
                                                  _umi_3p_matcher.Pool(),
                                                  /*use_sequence_offset=*/false,
                                                  kThreePrimeTieBreakPenalty);
    trim_info.umi_3p = best_match_3p.match.barcode_id;
    // we shouldn't need to adjust based on the length variable (actually refers to indel)
    // NOTE: May be overtrimming here
    pos_3p = best_match_3p.pos;
  }

  // check if trimming positions are valid
  if (pos_5p != kNoMatchPosition && pos_3p != kNoMatchPosition && pos_3p <= pos_5p) {
    // trim positions of 3p UMI is before 5p UMI which is invalid
    if (trim_info.umi_5p.has_value() && trim_info.umi_3p.has_value()) {
      // choose the 3p umi in this case and discard the 5p umi position as it is less reliable
      trim_info.umi_5p.reset();
      return {kNoMatchPosition, pos_3p};
    }
    if (trim_info.umi_5p.has_value() && !trim_info.umi_3p.has_value()) {
      // chose the 5p umi in this case and discard the 3p umi position
      return {pos_5p, kNoMatchPosition};
    }
    // Either only the 3p UMI was found, or no UMI was found at all. In both cases discard the 5p umi position
    // and pick the 3p umi position as it is more reliable.
    return {kNoMatchPosition, pos_3p};
  }

  return {pos_5p, pos_3p};
}

}  // namespace xoos::demux
