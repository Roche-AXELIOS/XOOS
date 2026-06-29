#include <hwy/highway.h>
#include <xoos/error/error.h>
#include <xoos/log/logging.h>

#include <array>

#include "adapters/duplex/demux-and-trim-duplex.h"
#include "metrics/duplex-metrics.h"
#include "sequence/matcher/bitap.h"

// Duplex HD sequences consists of two identical sequences, one in forward and one in reverse/complement orientation.
// The sequences are separated by a hairpin loop. The hairpin loop is a short sequence that is palindromic and only
// consists of a few bases (7).
// Adjacent to the hairpin, we find the SID associated with the sequence; for Duplex, the SID is 12 bases long and
// we have an SID at the 5p (just before the hairpin) and at the 3p (just after the hairpin).
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

// Take the candidate SIDs found by the LUTs and sort them into two buckets, one with 3p and one with
// 5p SIDs. Also sort them based on edit distance. Returns true if we found SIDs.
static bool SortSIDs(FixedReadRecord& record) {
  auto& trim_info{record.trim_info_duplex};
  auto nr_results{trim_info.unfiltered_matches_count};
  if (nr_results == 0) {
    return false;
  }

  auto& results{trim_info.unfiltered_matches};
  // First step: sort the markers in increasing edit distance, this will put the most reliable markers first.
  std::sort(results.begin(), results.begin() + nr_results);
  // Second step: copy the sorted markers into buckets, one for each type.
  for (int i = 0; i < DuplexMatch::kUnknown; ++i) {
    int count = 0;
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
static void WriteResults(FixedReadRecord& record, const CascadedLUTs& cascaded_luts, int min_error,
                         int best_match_index_5, int best_match_index_3) {
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

// This function calculates the likely 3p and 5p adapters by doing pairwise comparison and calculating
// an error metric that depends on their edit distance and the position of the hairpin.
int TrimDuplex::CalculateMinimumErrorHairpin(FixedReadRecord& record, const CascadedLUTs& cascaded_luts,
                                             int& best_match_index_5, int& best_match_index_3) const {
  auto& trim_info{record.trim_info_duplex};

  // We have a list of sorted matches. The most likely case: both SID3p and SID5p have been found and the
  // markers with the smallest edit distance have the same SIDs. Now loop over the SID5p markers, attempt to find
  // a corresponding SID3p marker, calculate an error metric and keep track of the markers that have the minimum
  // error. The combination with the minimal error wins.
  constexpr std::array kEditDistanceError = {0, 2, 4};
  int min_error = 10000;
  int error5p = min_error;
  u32 min_5p_index = 0;
  // Main loop: assume that 5p is present, handle absence of 3p within this loop
  for (u32 i = 0; i < trim_info.filtered_matches_count[DuplexMatch::kSID5p]; ++i) {
    // Loop over all 5p SIDs.
    auto& sid5p = trim_info.filtered_matches[DuplexMatch::kSID5p][i];
    // Calculate the length of the 5p SID.
    const auto len5p{cascaded_luts.Length(DuplexMatch::kSID5p, sid5p.match)};
    // We can now calculate an error term using the distance between the end of the hairpin marker and the end of the
    // SID5p. Using the length, that distance should be nominally 6.
    const auto distance_5p{std::abs(record.hairpin_pos - static_cast<int>(len5p) - static_cast<int>(sid5p.pos) - 6)};

    // now take into account the edit distance, for now let it count double.
    error5p = kEditDistanceError[sid5p.match.edist];
    if (error5p > 0) {
      error5p += distance_5p;  // only penalize for hairpin distance with edit distance > 0
    }
    if (error5p > min_error) {
      break;  // poor match, skip
    }
    if (error5p < min_error) {
      min_5p_index = i;  // usually not used
    }
    for (u32 j = 0; j < trim_info.filtered_matches_count[DuplexMatch::kSID3p]; ++j) {
      auto& sid3p = trim_info.filtered_matches[DuplexMatch::kSID3p][j];
      if (sid5p.match.barcode_id == sid3p.match.barcode_id) {  // same SID, so they'd be a pair
        // Calculate the error metric for this pair
        int distance_3p = sid3p.pos - record.hairpin_pos - 1;  // nominally, should be 1
        // Really penalize negative errors
        distance_3p = (distance_3p < 0) ? 1 - distance_3p : distance_3p;
        distance_3p += distance_3p;  // penalize any distance at the 3p adapter!
        auto total_error = error5p + distance_3p + kEditDistanceError[sid3p.match.edist];
        if (total_error < min_error) {
          min_error = total_error;
          best_match_index_5 = static_cast<int>(i);
          best_match_index_3 = static_cast<int>(j);
        }
      }
    }
  }

  // this is the case of R1 and 5' overhang being absent
  if (best_match_index_5 >= 0 && trim_info.filtered_matches[DuplexMatch::kSID5p][best_match_index_5].pos == 0) {
    best_match_index_5 = -1;
    // Set the min error to 20, that will cause the read to be skipped.
    return 20;
  }

  // Handle the case that we did not find any 3p adapter - disabled for now
  if (trim_info.filtered_matches_count[DuplexMatch::kSID3p] == 0 &&
      error5p <= 2) {  // we have a good 5p adapter at the correct position with respect to the hairpin. If we can't
    // find the 3p adapter, it is likely because the edit distance of the 3p adapter is too high.
    // So look for the 3p adapter using the bitap object of the 5p SID; this should find the 3p adapter.
    auto& best5p{trim_info.filtered_matches[DuplexMatch::kSID5p][min_5p_index]};
    auto& bitap{_sid_3p_bitap[best5p.match.barcode_id]};
    if (record.hairpin_pos <= (static_cast<int>(record.SeqLen()) - 14)) {
      // More than 14 bases left after hairpin, worth finding
      auto pos3p = bitap.Find(record.Seq(), record.hairpin_pos,
                              std::min(record.hairpin_pos + 24, static_cast<int>(record.SeqLen() - 1)));
      best_match_index_3 = 0;
      if (pos3p != -1) {  // we found a SID3p adapter with the same SID, we're done
        auto& best3p{trim_info.filtered_matches[DuplexMatch::kSID3p][best_match_index_3]};
        best3p.match.barcode_id = best5p.match.barcode_id;
        best3p.match.length = DuplexMatch::kZero;  // placeholder
        best3p.match.edist = 0;                    // placeholder, we're really lying here
        best3p.match.barcode_type = DuplexMatch::kSID3p;

        // What's the nominal length?
        const auto len{cascaded_luts.Length(DuplexMatch::kSID3p, best3p.match)};
        best3p.pos = 1 + pos3p - len;
        return error5p;  // this used to be return error5p; but I disabled this feature for now
      } else {           // can't find 3p adapter with acceptable edit distance - give up. Set the min error to 20,
        // that will cause the read to be skipped.
        return 20;
      }
    }
  }

  return min_error;
}

// This is the filter function to be used for plan A and plan B.
void TrimDuplex::FilterResults(FixedReadRecord& record, const CascadedLUTs& cascaded_luts) const {
  if (!SortSIDs(record)) {
    return;
  }
  int best_match_5 = -1, best_match_3 = -1;
  int min_error = CalculateMinimumErrorHairpin(record, cascaded_luts, best_match_5, best_match_3);
  // If both 3p and 5p have an edit distance of 2, the read should not be trusted. This implies that
  // a max error of 7 should be the max value; experiments with the synthetic dataset proved that using
  // a max error of 7 will not generate any false SIDs.
  if (min_error <= 7) {
    WriteResults(record, cascaded_luts, min_error, best_match_5, best_match_3);
  } else {  // reset the search
    record.trim_info_duplex.unfiltered_matches_count = 0;
    record.hairpin_pos = -1;
  }
}

void TrimDuplex::FindHairpinByStringSearch(FixedReadRecord& record) const {
  // Easiest method: look for the kmer "ATGTATG" in the sequence and look for SIDs next to it if found.
  // That sequence is not particularly unique, so prepare to keep searching if needed.

  size_t hairpin_pos = 0;

  std::string_view seq(record.Seq(), record.Seq() + record.SeqLen());
  while (hairpin_pos != std::string_view::npos && !record.trim_info_duplex.midadapter_range.has_value()) {
    hairpin_pos = seq.find(_loop_sequence, hairpin_pos);
    if (hairpin_pos != std::string_view::npos) {
      hairpin_pos += _loop_sequence.length() - 1;  // end position of kmer
      const size_t offset_hairpin = _loop_sequence.length() + _sid_5p_matcher.Pool().front().sequence.length();
      if (hairpin_pos >= offset_hairpin) {
        record.hairpin_pos = static_cast<int>(hairpin_pos);
        std::array<int64_t, 8> offsets = {
            record.hairpin_pos - static_cast<int>(offset_hairpin), record.hairpin_pos, 0, 0, 0, 0, 0, 0};

        // TODO: we can refactor to use std::array iterators instead of raw pointers if it doesn't affect performance.
        FindMarker(_cascaded_luts, record, offsets.data(), _mask_plan_a_full.data(), _types_plan_a_full.data());
        FilterResults(record, _cascaded_luts);
      }
    }
  }
}

void TrimDuplex::FindHairpinByGlobalSymmetry(FixedReadRecord& record) const {
  // This plan looks for "global" symmetry in the data, i.e. it uses SIMD processing to find symmetry in the
  // data. If we don't find symmetry (e.g. because of a truncated read that does not have enough
  // data), no problem, we have a plan C for that. But finding symmetry will put us on a fast path that usually
  // works. The calculated symmetry position is stored in sym; if sym is negative, we did not find symmetry.
  int sym = FindSymmetryPosition(record);

  // We are looking for the marker using symmetry; if the symmetry position is smaller than 22, the hardcoded
  // offset calculation will not work. Hence the '22', which is the "safe" value (19+3).
  if (sym > 22) {
    // Before looking for the SIDs, let's try to find the hairpin, which are 7 bases located in between the SIDs.
    // Let's make the search length 32 bases, that is pretty extensive.
    int start = sym - 19;
    int end = start + 31;
    int hairpin_pos = _kHairpin.Find(record.Seq(), start, end);
    record.hairpin_pos = hairpin_pos;

    const auto offset_hairpin =
        static_cast<int>(_loop_sequence.length() + _sid_5p_matcher.Pool().front().sequence.length());
    if (hairpin_pos >= offset_hairpin) {
      std::array<int64_t, 8> offsets = {record.hairpin_pos - offset_hairpin, record.hairpin_pos, 0, 0, 0, 0, 0, 0};

      // When looking for the markers, I'd like to have the highest possible likelihood of finding the match.
      // This partly depends on the 4-alignment of the data. For 3p: if pos_3p is 0-aligned, we'd like to start
      // searching at position pos_3p - 4.
      //         if pos_3p is 1-aligned, we'd like to start searching at position pos_3p - 3.
      //         if pos_3p is 2-aligned, we'd like to start searching at position pos_3p - 2.
      //         if pos_3p is 3-aligned, we'd like to start searching at position pos_3p - 5.
      // Similar things hold for the 5p adapter.

      FindMarker(_cascaded_luts, record, offsets.data(), _mask_plan_a_full.data(), _types_plan_a_full.data());
      FilterResults(record, _cascaded_luts);
    }
  }
}

void TrimDuplex::FindHairpinByLocalSymmetry(FixedReadRecord& record, TrimResults& results) const {
  // Plan C uses a more fine-grained search for symmetry (it does an exhaustive search for symmetry).
  // If symmetry is found, it performs a search for the 7-base middle k-mer; if it finds that, we'll
  // continue with a data path similar to plan A/B.
  // Because of the exhaustive search for symmetry, the success rate for 1+ rates is improved considerably.

  u32 begin_pos = 0;
  u32 offset = results.length;
  do {
    offset = static_cast<u32>(FindHairpinSliding(record, offset));
  } while (offset != 0);

  // Find the best matches for the point of symmetry. Calculate two positions; one with the maximum
  // score of shd + edit distance, one with maximum shd score
  std::array pos_3p = {0, 0};
  u32 max_score_total = 0, max_score_shd = 0;
  for (u32 i = begin_pos; i < results.length; ++i) {
    u32 score = record.match_values[i].shd_score + record.match_values[i].raw_score;
    if (score >= max_score_total) {
      max_score_total = score;
      pos_3p[0] = static_cast<int>(i);
    }
    if (record.match_values[i].shd_score > max_score_shd) {
      max_score_shd = record.match_values[i].shd_score;
      pos_3p[1] = static_cast<int>(i);
    }
  }

  constexpr int kMinScore = 16;
  int nr_candidates = std::abs(pos_3p[0] - pos_3p[1]) < 10 ? 1 : 2;
  for (int i = 0; i < nr_candidates; ++i) {
    if (record.match_values[pos_3p[i]].shd_score >= kMinScore) {
      record.hairpin_pos = static_cast<int>(pos_3p[i]);  // likely pos 3 position

      // We found a good match for the hairpin. Now we need to find the 3p and 5p adapter. Note that we might have
      // done an earlier search on a full read, so we should change the end marker positions if that's the case.
      int start = static_cast<int>(pos_3p[i]) - 8;
      int end = start + 20;
      if (end >= static_cast<int>(results.length)) {  // do not search beyond the end of the string
        end = static_cast<int>(results.length) - 1;
        start = end - 31;
      }
      int hairpin_pos = _kHairpin.Find(record.Seq(), start, end);
      if (hairpin_pos != -1) {
        record.hairpin_pos = hairpin_pos;  // found actual start of 3p
      }

      const int offset_hairpin =
          static_cast<int>(_loop_sequence.length() + _sid_5p_matcher.Pool().front().sequence.length());
      if (hairpin_pos >= offset_hairpin) {
        std::array<int64_t, 8> offsets = {record.hairpin_pos - offset_hairpin, record.hairpin_pos, 0, 0, 0, 0, 0, 0};

        FindMarker(_cascaded_luts, record, offsets.data(), _mask_plan_a_full.data(), _types_plan_a_full.data());
        FilterResults(record, _cascaded_luts);
      }
    }
    // return upon success
    if (record.trim_info_duplex.midadapter_range.has_value()) {
      return;
    }
  }

  // special case: if still not successful, it might be because the scoring failed to work correctly because
  // it is at the end of a read (I am using a sliding window of 32, so if the kmer is less than 32 away from the end,
  // we'll have issues.
  if (!record.trim_info_duplex.midadapter_range.has_value()) {
    int end = static_cast<int>(results.length) - 1;
    int start = std::max(0, end - 31);
    int hairpin_pos = _kHairpin.Find(record.Seq(), start, end);
    if (hairpin_pos != -1) {
      record.hairpin_pos = hairpin_pos;
    }
    const auto offset_hairpin =
        static_cast<int>(_loop_sequence.length() + _sid_5p_matcher.Pool().front().sequence.length());
    if (hairpin_pos >= offset_hairpin) {
      std::array<int64_t, 8> offsets = {record.hairpin_pos - offset_hairpin, record.hairpin_pos, 0, 0, 0, 0, 0, 0};

      FindMarker(_cascaded_luts, record, offsets.data(), _mask_plan_a_full.data(), _types_plan_a_full.data());
      FilterResults(record, _cascaded_luts);
    }
  }
}

void TrimDuplex::FindHairpin(FixedReadRecord& record, DuplexMetrics& metrics) const {  // NOLINT
  TrimResults result;
  result.length = record.SeqLen();

  auto& trim_info{record.trim_info_duplex};
#if 0
  // Please keep, allows me to debug a specific read.
  std::string_view seq{record.Name(), record.NameLen()};
  bool interesting{seq.find("fm-38:") != std::string::npos};
  if (interesting) {
    Logging:Debug(seq);
  }
#endif
  enum { kMiss = 0, kString, kGlobal, kLocal };

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
      hairpin_found = kGlobal;  // found hairpin by global search
    }
  } else {
    hairpin_found = kString;  // found hairpin by string search
  }

  // Finally, only write out the read if we found the midadapter.
  if (trim_info.midadapter_range.has_value()) {
    auto& match5p{trim_info.matches[DuplexMatch::kSID5p]};
    auto sid{match5p.match.barcode_id};
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
    record.hairpin_pos = -1;
  }

  if (record.hairpin_pos == -1) {
    metrics.unassigned_length_distr.AddCountToHistogram(result.length, 1);
    metrics.no_hairpin_length_distr.AddCountToHistogram(result.length, 1);
    metrics.unassigned_counts.no_hairpin_found += 1;
  }
}

}  // namespace xoos::demux
