#include "trim-duplex-umi.h"

#include <xoos/log/logging.h>
#include <xoos/types/int.h>

#include "adapters/duplex/duplex-match.h"
#include "adapters/duplex/stem/trim-duplex-stem.h"
#include "io/read-record.h"
#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

// Should make dynamic based on the LUT used in the search which varies
// The k-mer size of the LUT, we are hardcoding to 10 because the UMI length is 12 with 2 indels allowed
constexpr s32 kKmerLen = 10;

TrimDuplexUmi::TrimDuplexUmi(SeqMatcher sid_5p_matcher, SeqMatcher sid_3p_matcher, SeqMatcher start_matcher,
                             SeqMatcher end_matcher, const SeqMatcher& stem_5p_matcher,
                             const SeqMatcher& stem_3p_matcher, SeqMatcher umi_5p_matcher, SeqMatcher umi_3p_matcher,
                             std::string_view loop_sequence)
    : TrimDuplexStem(std::move(sid_5p_matcher), std::move(sid_3p_matcher), std::move(start_matcher),
                     std::move(end_matcher), stem_5p_matcher, stem_3p_matcher, loop_sequence),
      _umi_5p_matcher(std::move(umi_5p_matcher)),
      _umi_3p_matcher(std::move(umi_3p_matcher)),
      _umi_lut(CascadedLUTs(_umi_5p_matcher, _umi_3p_matcher, DuplexMatch::BarcodeType::kUMI5p,
                            DuplexMatch::BarcodeType::kUMI3p)),
      _mask_umi5p(static_cast<s64>((1ul << (2 * _umi_lut.MaxLength(DuplexMatch::kUMI5p))) - 1ul)),
      _mask_umi3p(static_cast<s64>((1ul << (2 * _umi_lut.MaxLength(DuplexMatch::kUMI3p))) - 1ul)),
      // Array of masks for LUT search
      _lut_offset_umi_masks(_mask_umi5p, _mask_umi5p, _mask_umi3p, _mask_umi3p, 0, 0, 0, 0),
      // Array of barcode types for LUT search
      _lut_offset_umi_types(DuplexMatch::BarcodeType::kUMI5p, DuplexMatch::BarcodeType::kUMI5p,
                            DuplexMatch::BarcodeType::kUMI3p, DuplexMatch::BarcodeType::kUMI3p,
                            DuplexMatch::BarcodeType::kUnknown, DuplexMatch::BarcodeType::kUnknown,
                            DuplexMatch::BarcodeType::kUnknown, DuplexMatch::BarcodeType::kUnknown) {}

/**
 * @brief Sorts the UMI results in the FixedReadRecord.
 * This function sorts the unfiltered matches by edit distance and then filters them into buckets based on their type.
 * Filters out matches that are out of bounds of the read length
 * @param record a FixedReadRecord containing the unfiltered matches to sort
 */
static void SortUMIs(FixedReadRecord& record, const s64 expect_pos_5p, const s64 expect_pos_3p) {
  const auto average_search_pos = (expect_pos_5p + expect_pos_3p) / 2;

  auto& trim_info = record.trim_info_duplex;
  auto nr_results = trim_info.unfiltered_matches_count;
  if (nr_results == 0) {
    return;
  }

  auto& results{trim_info.unfiltered_matches};
  // First step: sort the markers in increasing edit distance, this will put the most reliable markers first.
  std::sort(results.begin(), results.begin() + nr_results);
  s32 count_5p = 0;
  s32 count_3p = 0;
  for (u32 i = 0; i < nr_results; ++i) {
    // Only include matches that are within the read boundaries
    if (results[i].pos > record.consensus_seq_len - kKmerLen) {
      continue;
    }
    if (results[i].match.barcode_type == DuplexMatch::kUMI5p) {
      // filter out matches that are beyond the average search position
      if (results[i].pos > average_search_pos) {
        continue;
      }
      trim_info.filtered_matches[DuplexMatch::kUMI5p][count_5p++] = results[i];
    } else if (results[i].match.barcode_type == DuplexMatch::kUMI3p) {
      // filter out matches that are before the average search position
      if (results[i].pos < average_search_pos) {
        continue;
      }
      trim_info.filtered_matches[DuplexMatch::kUMI3p][count_3p++] = results[i];
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
std::pair<s32, s32> TrimDuplexUmi::FindUMI(FixedReadRecord& record) const {
  // Find the 5' UMI
  // Find adapter sequence to anchor search using bitap
  // Sometimes the position ends early due to mismatches so this can still let some of the UMI into the sequence
  // This is  fixed if we detect a valid 5' umi match bit it breaks if too far to search or if the umi is not found
  // TODO: Ideally need to 1) get info if bitap match was missing trailing bases and 2) adjust the position
  const auto seq = std::string_view(record.consensus_buffer, record.consensus_seq_len);
  // Find adapter sequence to anchor search using bitap
  // search using the 3' stem since it is more likely to exist
  s32 pos_3p = stem_3p_bitap.ReverseScan(seq);
  // if we found the 3' stem, search for the 5' stem upstream
  // right to left to compensate for artefacts (repeated stems)
  s32 pos_5p = stem_5p_bitap.ReverseScan(seq, -1, pos_3p);

  auto& trim_info = record.trim_info_duplex;

  // reset the UMI info in the record
  trim_info.umi_5p.reset();
  trim_info.umi_3p.reset();

  // reset the matches because we would have existing matches from the hairpin search
  trim_info.unfiltered_matches_count = 0;
  trim_info.filtered_matches_count[DuplexMatch::kUMI3p] = 0;
  trim_info.filtered_matches_count[DuplexMatch::kUMI5p] = 0;

  // TODO: the 6 is from the runway length which could in theory be retrieved
  const s64 assumed_umi_len_5p_4bp =
      static_cast<s32>(stem_5p_matcher.Pool().front().sequence.length() + 6 - _umi_info.umi_5p_offset_4bp);
  const s64 assumed_umi_len_5p_6bp =
      static_cast<s32>(stem_5p_matcher.Pool().front().sequence.length() + 6 - _umi_info.umi_5p_offset_6bp);
  const s64 assumed_umi_len_3p_4bp =
      static_cast<s32>(record.consensus_seq_len - stem_3p_matcher.Pool().front().sequence.length() - kUmiLen1);
  const s64 assumed_umi_len_3p_6bp =
      static_cast<s32>(record.consensus_seq_len - stem_3p_matcher.Pool().front().sequence.length() - kUmiLen2);

  // we assign the search position to the found location, otherwise we try to search based on theoretical location
  const s64 search_pos_5p_4bp =
      pos_5p == -1 ? assumed_umi_len_5p_4bp : std::max(0, pos_5p - _umi_info.umi_5p_offset_4bp);
  const s64 search_pos_5p_6bp =
      pos_5p == -1 ? assumed_umi_len_5p_6bp : std::max(0, pos_5p - _umi_info.umi_5p_offset_6bp);
  const s64 search_pos_3p_4bp = pos_3p == -1 ? assumed_umi_len_3p_4bp : std::max(0, pos_3p - kUmiLen1);
  const s64 search_pos_3p_6bp = pos_3p == -1 ? assumed_umi_len_3p_6bp : std::max(0, pos_3p - kUmiLen2);

  // using SIMD search, it can search 8 positions at once, it searches every 8 bits (4 bases) the offset points to
  // We redundantly search the 5' 0, 4, 8, 12 positions more to capture cases where the 5' stem is not
  int64_t offsets[8] = {search_pos_5p_4bp, search_pos_5p_6bp, search_pos_3p_4bp, search_pos_3p_6bp, 0, 4, 8,
                        search_pos_3p_6bp};

  // Convert the consensus sequence to 2-bit representation and store in the record
  record.InitTwoBit(reinterpret_cast<const u8*>(record.consensus_buffer), static_cast<u32>(record.consensus_seq_len));

  FindMarker(_umi_lut, record, offsets, _lut_offset_umi_masks, _lut_offset_umi_types);
  SortUMIs(record, (search_pos_5p_4bp + search_pos_5p_6bp) / 2, (search_pos_3p_4bp + search_pos_3p_6bp) / 2);
  if (trim_info.filtered_matches_count[DuplexMatch::kUMI5p] > 0) {
    // The best match is the first one after sorting by edit distance
    auto best_match_5p = trim_info.filtered_matches[DuplexMatch::kUMI5p][0];
    // if we have an expected position, try to find the match that is closest to it
    if (pos_5p != -1) {
      // iterate over the matches and find the one that is closest to the expected position
      for (u32 i = 1; i < trim_info.filtered_matches_count[DuplexMatch::kUMI5p]; ++i) {
        const auto match = trim_info.filtered_matches[DuplexMatch::kUMI5p][i];
        // the difference between the expected position and the found position adjusted for the UMI length
        auto best = _umi_5p_matcher.Pool()[best_match_5p.match.barcode_id];
        auto diff_best =
            std::abs(best_match_5p.pos + static_cast<s32>(best.sequence.length() - best.name.length() - 1) - pos_5p);
        auto current = _umi_5p_matcher.Pool()[match.match.barcode_id];
        auto diff_current =
            std::abs(match.pos + static_cast<s32>(current.sequence.length() - current.name.length() - 1) - pos_5p);
        if (diff_current < diff_best) {
          // if difference is smaller, take this match but only if the net error if the difference exists is not larger
          if (diff_current + match.match.edist < diff_best + best_match_5p.match.edist) {
            best_match_5p = match;
          }
        }
      }
    }

    auto umi_5p = best_match_5p.match.barcode_id;
    trim_info.umi_5p = umi_5p;

    // adjust the position based on the length variable (actually refers to indel)
    auto mod5p = best_match_5p.pos + _umi_info.length_5p - 1;
    switch (best_match_5p.match.length) {
      case DuplexMatch::kPlus2:
        mod5p += 2;
        break;
      case DuplexMatch::kPlus1:
        mod5p += 1;
        break;
      case DuplexMatch::kZero:
        break;
      case DuplexMatch::kMinus1:
        mod5p -= 1;
        break;
      case DuplexMatch::kMinus2:
        mod5p -= 2;
        break;
      default:
        throw error::Error("Unknown UMI indel size: {}", static_cast<s32>(best_match_5p.match.length));
    }
    pos_5p = mod5p;
  }
  if (trim_info.filtered_matches_count[DuplexMatch::kUMI3p] > 0) {
    // The best match is the first one after sorting by edit distance
    auto best_match_3p = trim_info.filtered_matches[DuplexMatch::kUMI3p][0];
    // if we have an expected position, try to find the match that is closest to it
    if (pos_3p != -1) {
      // iterate over the matches and find the one that is closest to the expected position
      for (u32 i = 1; i < trim_info.filtered_matches_count[DuplexMatch::kUMI3p]; ++i) {
        const auto match = trim_info.filtered_matches[DuplexMatch::kUMI3p][i];
        // the difference between the expected position and the found position adjusted for the UMI length
        auto best = _umi_3p_matcher.Pool()[best_match_3p.match.barcode_id];
        auto diff_best = std::abs(best_match_3p.pos + static_cast<s32>(best.name.length()) - pos_3p);
        auto current = _umi_3p_matcher.Pool()[match.match.barcode_id];
        auto diff_current = std::abs(match.pos + static_cast<s32>(current.name.length()) - pos_3p);
        if (diff_current < diff_best) {
          // if difference is smaller, take this match but only if the net error if the difference exists is not larger
          if (diff_current + match.match.edist + 1 < diff_best + best_match_3p.match.edist) {
            best_match_3p = match;
          }
        }
      }
    }

    auto umi_3p = best_match_3p.match.barcode_id;
    trim_info.umi_3p = umi_3p;
    // we shouldn't need to adjust based on the length variable (actually refers to indel)
    // TODO: May be overtrimming here
    pos_3p = best_match_3p.pos - 1;
  }
  // check if trimming positions are valid
  if (pos_5p != -1 && pos_3p != -1 && pos_3p <= pos_5p) {
    // trim positions of 3p UMI is before 5p UMI which is invalid
    if (trim_info.umi_5p.has_value() && trim_info.umi_3p.has_value()) {
      // choose the 3p umi in this case and discard the 5p umi position as it is less reliable
      trim_info.umi_5p.reset();
      return {-1, pos_3p};
    }
    if (trim_info.umi_5p.has_value() && !trim_info.umi_3p.has_value()) {
      // chose the 5p umi in this case and discard the 3p umi position
      return {pos_5p, -1};
    }
    if (!trim_info.umi_5p.has_value() && trim_info.umi_3p.has_value()) {
      // we only found the 3p umi so discard the 5p umi position
      return {-1, pos_3p};
    }
    // if no UMI was found at all still return the found positions but pick the 3p umi position as it is more reliable
    return {-1, pos_3p};
  }

  return {pos_5p, pos_3p};
}

}  // namespace xoos::demux
