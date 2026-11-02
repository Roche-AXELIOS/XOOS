
#include <hwy/highway.h>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>

#include "adapters/duplex/hairpin-finder.h"
#include "metrics/duplex-metrics.h"
#include "sequence/matcher/bitap.h"
#include "sequence/matcher/match-position-states.h"

// Duplex HD sequences consist of two identical sequences, one in forward and one in reverse/complement orientation.
// The sequences are separated by a hairpin loop. The hairpin loop is a short palindromic sequence.
// Adjacent to the loop, we find the SID associated with the sequence; for Duplex, we have an SID at the
// 5p (just before the loop) and at the 3p (just after the loop). Together, SID5p + loop + SID3p make up the hairpin.
//
// Fortunately, there are plenty of sequences we can use for finding symmetry in the read.
// Notably, with a full read the 5p half is identical to the 3p half (aside from the fact that the 3p half is
// reverse-complemented). This means that we can look for almost any sequence in the 5p half and expect to find it in
// the 3p half; this allows us to pick a longer sequence (eg 16 bases).
//
// However, we need to be careful with the position of the sequence. Reads can be truncated, which will make it
// impossible to find a match at the 3p end if it is not present. So, we can solve this by picking a sequence at the
// 3p end and looking for it at the 5p end. This way, it would maximize the likelihood of finding the sequence.
//
// For this approach to succeed, we're counting on the search to be done extremely fast - this is where AVX512 comes
// into play. We can perform 8 searches for the price of one, and we can do it that in a couple of clock cycles because
// I can implement it such that all the data lives in registers. By repeating the search several times, we are able to
// evaluate 256 positions at the 5p end of the read, which usually should be enough to find a match.
//
// If not successful, we can do additional searches using a different approach - with a CPU, every read can be
// treated individually.

namespace xoos::demux {

namespace {

constexpr std::array<u8, 3> kEditDistancePenaltyByEdist = {0, 2, 4};
constexpr u8 kMaxRecoverable5pError = 2;
constexpr u8 kMaxAcceptedSidPairError = 7;
constexpr u8 kRejectedReadError = kMaxAcceptedSidPairError + 1;
constexpr u8 kSimdSearchLaneCount = 8;
constexpr u8 kCandidateAgreementThreshold = 10;
constexpr u8 kMinLocalSymmetryScore = 16;
/// @brief Maximum edit distance supported by the SID LUT-based search.
constexpr u8 kSidMaxEditDistance = kEditDistancePenaltyByEdist.size() - 1;
/// @brief The loop search window is sized as this multiple of the loop length.
constexpr s32 kLoopSearchWindowMultiplier = 3;

}  // namespace

// Take the candidate SIDs found by the LUTs and sort them into two buckets, one with 3p and one with
// 5p SIDs. Also sort them based on edit distance. Returns true if we found SIDs.
static bool SortSIDs(FixedReadRecord& record) {
  auto& trim_info{record.trim_info_duplex};
  const auto nr_results{trim_info.unfiltered_matches_count};
  if (nr_results == 0) {
    return false;
  }

  auto& results{trim_info.unfiltered_matches};
  // First step: sort the markers in increasing edit distance, this will put the most reliable markers first.
  std::sort(results.begin(), results.begin() + nr_results);
  // Second step: copy the sorted markers into buckets, one for each type.
  for (u32 i = 0; i < DuplexMatch::kUnknown; ++i) {
    u32 count = 0;
    for (u32 j = 0; j < nr_results; ++j) {
      if (results[j].match.barcode_type == i) {
        trim_info.filtered_matches[i][count] = results[j];
        ++count;
      }
    }
    trim_info.filtered_matches_count[i] = count;
  }
  return true;
}

// Write out results, notably when the midadapter was found.
static void WriteResults(FixedReadRecord& record,
                         const CascadedLUTs& cascaded_luts,
                         const s32 min_error,
                         const s32 best_match_index_5,
                         const s32 best_match_index_3) {
  auto& trim_info{record.trim_info_duplex};

  record.error_metric = min_error;
  trim_info.matches[DuplexMatch::kSID5p] = trim_info.filtered_matches[DuplexMatch::kSID5p][best_match_index_5];
  record.trim_info_duplex.matches[DuplexMatch::kSID3p] =
      trim_info.filtered_matches[DuplexMatch::kSID3p][best_match_index_3];
  const auto len{cascaded_luts.Length(DuplexMatch::kSID3p, trim_info.matches[DuplexMatch::kSID3p].match)};
  // Note: increased the mid adapter size by 1 at 5p and 3p to account for overhang.
  trim_info.midadapter_range = LociRange(trim_info.matches[DuplexMatch::kSID5p].pos - 1,
                                         static_cast<u32>(trim_info.matches[DuplexMatch::kSID3p].pos + len));
}

void HairpinFinder::FindSidsAroundLoop(FixedReadRecord& record, const s32 loop_start, const s32 loop_end) const {
  // We search two pairs of SID anchors in a single SIMD pass (4 of 8 available lanes).
  //
  // Lanes 0–1 ("indel-aware"): both anchored from the actual matched loop_start.
  // 5p is shifted left by kSidMaxEditDistance to absorb leftward SID drift from insertions.
  // 3p uses loop_start + nominal loop_length to project where the 3p SID should be, mirroring
  // how the nominal pair projects from loop_end. This keeps both anchors in the indel-aware
  // pair derived from the same reference point (loop_start).
  const s32 sid_5p_indel = loop_start - _sid_5p_length - kSidMaxEditDistance;
  const s32 sid_3p_indel = loop_start + _loop_length;

  // Lanes 2–3 ("nominal"): anchored from loop_end using the nominal loop length, matching the
  // window placement that covers the majority of reads where the loop has no indels. Because the
  // SIMD gather quantizes offsets to 4-aligned bases, even a 1-base shift between these two
  // strategies can move the quantized window start by 4, recovering SIDs at the right edge.
  const s32 sid_5p_nominal = loop_end - _loop_length - _sid_5p_length;
  const s32 sid_3p_nominal = loop_end;

  if (sid_5p_indel < 0 && sid_5p_nominal < 0) {
    return;
  }
  const std::array<s64, kSimdSearchLaneCount> sid_search_starts = {
      std::max<s64>(sid_5p_indel, 0), sid_3p_indel, std::max<s64>(sid_5p_nominal, 0), sid_3p_nominal, 0, 0, 0, 0};
  FindMarker(_cascaded_luts, record, sid_search_starts.data(), _mask_plan_a_full.data(), _types_plan_a_full.data());
  FilterResults(record, _cascaded_luts);
}

// This function calculates the likely 3p and 5p adapters by doing pairwise comparison and calculating
// an error metric that depends on their edit distance and the position of the loop.
s32 HairpinFinder::CalculateMinimumErrorHairpin(FixedReadRecord& record,
                                                const CascadedLUTs& cascaded_luts,
                                                s32& best_match_index_5,
                                                s32& best_match_index_3) const {
  const auto& trim_info{record.trim_info_duplex};

  // Phase 1: Score all 5p/3p SID pairs.
  const auto score = ScoreSidPairs(trim_info, cascaded_luts, record.loop_end_pos, _loop_length);
  best_match_index_5 = score.best_5p;
  best_match_index_3 = score.best_3p;

  // Phase 2: Reject reads with absent R1/5' overhang.
  if (best_match_index_5 >= 0 && trim_info.filtered_matches[DuplexMatch::kSID5p][best_match_index_5].pos == 0) {
    best_match_index_5 = kNoMatchPosition;
    return kRejectedReadError;
  }

  // Phase 3: Bitap fallback for missing 3p.
  return TryBitapFallback3p(record, score, best_match_index_3).value_or(score.min_error);
}

HairpinFinder::SidPairScore HairpinFinder::ScoreSidPairs(const TrimInfoDuplex& trim_info,
                                                         const CascadedLUTs& cascaded_luts,
                                                         const s32 loop_end_pos,
                                                         const s32 loop_length) {
  s32 min_error = std::numeric_limits<s32>::max();
  s32 error5p = min_error;
  u32 min_5p_index = 0;
  s32 best_5 = kNoMatchPosition;
  s32 best_3 = kNoMatchPosition;

  // Expected distances from loop_end_pos to SID boundaries:
  // - SID 5p end to loop_end = loop_length - 1 (SID 5p ends at loop_start, loop_end is inclusive)
  // - loop_end to SID 3p start = 1 (SID 3p starts at loop_end + 1, which is 1 past the inclusive end)
  const s32 loop_end_to_sid5p_gap = loop_length - 1;
  constexpr s32 kLoopEndToSid3pGap = 1;

  for (u32 i = 0; i < trim_info.filtered_matches_count[DuplexMatch::kSID5p]; ++i) {
    const auto& sid5p = trim_info.filtered_matches[DuplexMatch::kSID5p][i];
    const auto len5p{cascaded_luts.Length(DuplexMatch::kSID5p, sid5p.match)};
    const auto distance_5p{std::abs(loop_end_pos - static_cast<s32>(len5p + sid5p.pos) - loop_end_to_sid5p_gap)};

    error5p = kEditDistancePenaltyByEdist[sid5p.match.edist];
    if (error5p > 0) {
      error5p += distance_5p;
    }
    if (error5p > min_error) {
      break;
    }
    if (error5p < min_error) {
      min_5p_index = i;
    }
    for (u32 j = 0; j < trim_info.filtered_matches_count[DuplexMatch::kSID3p]; ++j) {
      const auto& sid3p = trim_info.filtered_matches[DuplexMatch::kSID3p][j];
      if (sid5p.match.barcode_id == sid3p.match.barcode_id) {
        const s32 raw_3p = sid3p.pos - loop_end_pos - kLoopEndToSid3pGap;
        // Double the penalty; negative offsets are extra penalised.
        const s32 distance_3p = 2 * ((raw_3p < 0) ? kLoopEndToSid3pGap - raw_3p : raw_3p);
        const auto total_error = error5p + distance_3p + kEditDistancePenaltyByEdist[sid3p.match.edist];
        if (total_error < min_error) {
          min_error = total_error;
          best_5 = static_cast<s32>(i);
          best_3 = static_cast<s32>(j);
        }
      }
    }
  }
  return {min_error, error5p, min_5p_index, best_5, best_3};
}

std::optional<s32> HairpinFinder::TryBitapFallback3p(FixedReadRecord& record,
                                                     const SidPairScore& score,
                                                     s32& best_match_index_3) const {
  auto& trim_info{record.trim_info_duplex};

  if (trim_info.filtered_matches_count[DuplexMatch::kSID3p] != 0 ||
      !std::cmp_less_equal(score.error5p, kMaxRecoverable5pError)) {
    return std::nullopt;
  }

  const auto& best5p{trim_info.filtered_matches[DuplexMatch::kSID5p][score.min_5p_index]};
  const auto& bitap_fw{_sid_3p_bitap[best5p.match.barcode_id]};
  const auto& bitap_bw{_sid_3p_bw[best5p.match.barcode_id]};

  // We need enough remaining bases to find the 3p SID: at least the SID length + a small margin.
  const s32 min_remaining_bases = _sid_3p_length + 2;
  const s32 fallback_search_span = 2 * _sid_3p_length;

  if (record.loop_end_pos > (static_cast<s32>(record.SeqLen()) - min_remaining_bases)) {
    return std::nullopt;
  }

  const auto [pos3p_start, pos3p_end] =
      bitap_fw.FindStartEnd(record.Seq(),
                            record.loop_end_pos,
                            std::min(record.loop_end_pos + fallback_search_span, static_cast<s32>(record.SeqLen() - 1)),
                            bitap_bw);

  best_match_index_3 = 0;
  if (pos3p_start == kNoMatchPosition) {
    return kRejectedReadError;
  }

  auto& best3p{trim_info.filtered_matches[DuplexMatch::kSID3p][best_match_index_3]};
  best3p.match.barcode_id = best5p.match.barcode_id;
  best3p.match.barcode_type = DuplexMatch::kSID3p;
  best3p.match.edist = 0;
  const s32 delta = (pos3p_end - pos3p_start + 1) - _sid_3p_length;
  best3p.match.length = DuplexMatch::LengthFromDelta(delta);
  best3p.pos = static_cast<u16>(pos3p_start);
  return score.error5p;
}

// This is the filter function to be used for plan A and plan B.
void HairpinFinder::FilterResults(FixedReadRecord& record, const CascadedLUTs& cascaded_luts) const {
  if (!SortSIDs(record)) {
    return;
  }
  auto best_match_5 = kNoMatchPosition;
  auto best_match_3 = kNoMatchPosition;
  const auto min_error = CalculateMinimumErrorHairpin(record, cascaded_luts, best_match_5, best_match_3);
  // If both 3p and 5p have an edit distance of 2, the read should not be trusted. This implies that
  // a max error of 7 should be the max value; experiments with the synthetic dataset proved that using
  // a max error of 7 will not generate any false SIDs.
  if (std::cmp_less_equal(min_error, kMaxAcceptedSidPairError)) {
    WriteResults(record, cascaded_luts, min_error, best_match_5, best_match_3);
  } else {  // reset the search
    record.trim_info_duplex.unfiltered_matches_count = 0;
    record.loop_end_pos = kNoMatchPosition;
    record.loop_start_pos = kNoMatchPosition;
  }
}

void HairpinFinder::FindHairpinByStringSearch(FixedReadRecord& record) const {
  // Easiest method: look for the loop sequence in the read and look for SIDs next to it if found.
  // The loop sequence is not particularly unique, so prepare to keep searching if needed.

  size_t loop_pos = 0;

  const std::string_view seq(record.Seq(), record.Seq() + record.SeqLen());
  while (loop_pos != std::string_view::npos && !record.trim_info_duplex.midadapter_range.has_value()) {
    loop_pos = seq.find(_loop_sequence, loop_pos);
    if (loop_pos != std::string_view::npos) {
      record.loop_start_pos = static_cast<s32>(loop_pos);
      // end position of kmer
      loop_pos += _loop_sequence.length() - 1;
      const auto min_loop_offset = static_cast<size_t>(_loop_length) + static_cast<size_t>(_sid_5p_length);
      if (loop_pos >= min_loop_offset) {
        record.loop_end_pos = static_cast<s32>(loop_pos);
        FindSidsAroundLoop(record, record.loop_start_pos, record.loop_end_pos);
      }
    }
  }
}

void HairpinFinder::FindHairpinByGlobalSymmetry(FixedReadRecord& record) const {
  // This plan looks for "global" symmetry in the data, i.e. it uses SIMD processing to find symmetry in the
  // data. If we don't find symmetry (e.g. because of a truncated read that does not have enough
  // data), no problem, we have a plan C for that. But finding symmetry will put us on a fast path that usually
  // works. The calculated symmetry position is stored in sym; if sym is negative, we did not find symmetry.
  const auto sym = FindSymmetryPosition(record);

  // The minimum symmetry position must be large enough that `start = sym - search_start_offset` stays non-negative
  // with margin. search_start_offset = sid_5p_length + loop_length; add a small safety margin of 3.
  const s32 search_start_offset = _sid_5p_length + _loop_length;
  const s32 min_sym_position = search_start_offset + 3;

  if (sym > min_sym_position) {
    // Before looking for the SIDs, let's try to find the loop located in between the SIDs.
    const s32 start = sym - search_start_offset;
    // Clamp `end` to the read bounds so the bitap search never reads outside the sequence buffer.
    // `FindSymmetryPosition` only guarantees `sym < SeqLen`, so `end` can spill past `SeqLen - 1`
    // when symmetry lands near the 3p tail of the read.
    const s32 end = std::min(start + _loop_search_window - 1, static_cast<s32>(record.SeqLen()) - 1);
    const auto [loop_start, loop_end] = _loop_fw.FindStartEnd(record.Seq(), start, end, _loop_bw);
    if (loop_end == kNoMatchPosition) {
      return;
    }
    record.loop_start_pos = loop_start;
    record.loop_end_pos = loop_end;

    // When looking for the markers, I'd like to have the highest possible likelihood of finding the match.
    // This partly depends on the 4-alignment of the data. For 3p: if pos_3p is 0-aligned, we'd like to start
    // searching at position pos_3p - 4.
    //         if pos_3p is 1-aligned, we'd like to start searching at position pos_3p - 3.
    //         if pos_3p is 2-aligned, we'd like to start searching at position pos_3p - 2.
    //         if pos_3p is 3-aligned, we'd like to start searching at position pos_3p - 5.
    // Similar things hold for the 5p adapter.
    FindSidsAroundLoop(record, loop_start, loop_end);
  }
}

void HairpinFinder::FindHairpinByLocalSymmetry(FixedReadRecord& record, const TrimResults& results) const {
  // Plan C uses a more fine-grained search for symmetry (it does an exhaustive search for symmetry).
  // If symmetry is found, it performs a search for the loop k-mer; if it finds that, we'll
  // continue with a data path similar to plan A/B.
  // Because of the exhaustive search for symmetry, the success rate for 1+ rates is improved considerably.

  // Backtrack from the candidate position to give some room before the expected loop start.
  const s32 local_search_backtrack = _loop_length + 1;
  // The local symmetry candidate is more precise, so the primary window can be narrower.
  const s32 local_search_window_length =
      std::min(kLoopSearchWindowMultiplier * _loop_length, Bitap<2>::kQueryWindowSize);

  const u32 begin_pos = 0;
  u32 offset = results.length;
  do {
    offset = static_cast<u32>(FindHairpinSliding(record, offset));
  } while (offset != 0);

  // Find the best matches for the point of symmetry. Calculate two positions; one with the maximum
  // score of shd + edit distance, one with maximum shd score
  std::array<u32, 2> pos_3p = {0, 0};
  u32 max_score_total = 0;
  u32 max_score_shd = 0;
  for (u32 i = begin_pos; i < results.length; ++i) {
    const u32 score = record.match_values[i].shd_score + record.match_values[i].raw_score;
    if (score >= max_score_total) {
      max_score_total = score;
      pos_3p[0] = i;
    }
    if (record.match_values[i].shd_score > max_score_shd) {
      max_score_shd = record.match_values[i].shd_score;
      pos_3p[1] = i;
    }
  }

  const u32 candidate_distance = pos_3p[0] > pos_3p[1] ? pos_3p[0] - pos_3p[1] : pos_3p[1] - pos_3p[0];
  const u32 nr_candidates = candidate_distance < kCandidateAgreementThreshold ? 1 : static_cast<u32>(pos_3p.size());
  for (u32 i = 0; i < nr_candidates; ++i) {
    if (record.match_values[pos_3p[i]].shd_score >= kMinLocalSymmetryScore) {
      // likely pos 3 position
      record.loop_end_pos = static_cast<s32>(pos_3p[i]);

      // We have a likely loop-end candidate. Now we need to find the 3p and 5p SIDs to confirm the hairpin. We might
      // have done an earlier search on a full read, so we should change the end marker positions if that's the case.
      auto start = static_cast<s32>(pos_3p[i]) - local_search_backtrack;
      auto end = start + local_search_window_length - 1;
      if (end >= static_cast<s32>(results.length)) {  // do not search beyond the end of the string
        end = static_cast<s32>(results.length) - 1;
        start = end - (_loop_search_window - 1);
      }
      const auto [loop_start, loop_end] = _loop_fw.FindStartEnd(record.Seq(), start, end, _loop_bw);
      if (loop_end != kNoMatchPosition) {
        record.loop_start_pos = loop_start;
        record.loop_end_pos = loop_end;
        FindSidsAroundLoop(record, loop_start, loop_end);
      }
    }
    // return upon success
    if (record.trim_info_duplex.midadapter_range.has_value()) {
      return;
    }
  }

  // Special case: if still not successful, it might be because the scoring failed to work correctly because
  // it is at the end of a read (using a sliding window of 32, so if the k-mer is less than 32 away from the
  // end, we'll have issues).
  if (!record.trim_info_duplex.midadapter_range.has_value()) {
    const auto end = static_cast<s32>(results.length) - 1;
    const auto start = std::max(0, end - (_loop_search_window - 1));
    const auto [loop_start, loop_end] = _loop_fw.FindStartEnd(record.Seq(), start, end, _loop_bw);
    if (loop_end != kNoMatchPosition) {
      record.loop_start_pos = loop_start;
      record.loop_end_pos = loop_end;
      FindSidsAroundLoop(record, loop_start, loop_end);
    }
  }
}

void HairpinFinder::FindHairpin(FixedReadRecord& record, DuplexMetrics& metrics) const {  // NOLINT
  TrimResults result;
  result.length = record.SeqLen();

  auto& trim_info{record.trim_info_duplex};
  enum class HairpinFindMethod {
    kMiss = 0,
    kString,
    kGlobal,
    kLocal
  };
  using enum HairpinFindMethod;

  auto hairpin_found = kMiss;
  // Look for kmer using std::string::find (edit distance = 0 - fastest approach)
  FindHairpinByStringSearch(record);

  if (!trim_info.midadapter_range.has_value()) {
    // Use global symmetry to find likely kmer position
    FindHairpinByGlobalSymmetry(record);

    if (!trim_info.midadapter_range.has_value()) {
      // Do a brute-force search for the mid adapter using local symmetry. Because earlier plans usually work,
      // reads, we can afford to be a bit more wasteful in the search.
      FindHairpinByLocalSymmetry(record, result);
      if (trim_info.midadapter_range.has_value()) {
        hairpin_found = kLocal;
      }
    } else {
      // found hairpin by global search
      hairpin_found = kGlobal;
    }
  } else {
    // found hairpin by string search
    hairpin_found = kString;
  }

  // Finally, only write out the read if we found the midadapter.
  if (trim_info.midadapter_range.has_value()) {
    const auto& match5p{trim_info.matches[DuplexMatch::kSID5p]};
    const auto sid{match5p.match.barcode_id};
    trim_info.duplex_status = TrimInfoDuplex::DuplexStatus::kMidAdapterFound;
    switch (hairpin_found) {
      case kString:
        metrics.midadapter_counts.found_by_string_compare[sid] += 1;
        break;
      case kGlobal:
        metrics.midadapter_counts.found_by_global_symmetry[sid] += 1;
        break;
      case kLocal:
        metrics.midadapter_counts.found_by_local_symmetry[sid] += 1;
        break;
      default:
        break;
    }
  } else {
    trim_info.duplex_status = TrimInfoDuplex::DuplexStatus::kZeroPlus;
    record.loop_end_pos = kNoMatchPosition;
    record.loop_start_pos = kNoMatchPosition;
  }

  if (record.loop_end_pos == kNoMatchPosition) {
    metrics.unassigned_length_distr.AddCountToHistogram(result.length, 1);
    metrics.no_hairpin_length_distr.AddCountToHistogram(result.length, 1);
    metrics.unassigned_counts.no_hairpin_found += 1;
    metrics.unassigned_counts.raw_bases += record.SeqLen();
  }
}

}  // namespace xoos::demux
