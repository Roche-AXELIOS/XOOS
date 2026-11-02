#include "trim3p-ys.h"

#include <algorithm>
#include <cassert>
#include <string_view>

#include <xoos/enum/enum-util.h>
#include <xoos/error/error.h>

#include "sequence/matcher/match-position-states.h"
#include "sequence/sequence-enums.h"

namespace xoos::demux {
Trim3pInfoYs::Trim3pInfoYs(const u32 insert_end) : insert_end(insert_end) {
}

// Preferred seed position within the 3' stem.  Among the nine possible 5bp
// windows in the 13bp stem, the window at offset 3 (GTCGT) has the lowest
// hit rate in real data (~60% of reads vs 70% for TTTGT), keeping the
// number of localized Bitap calls low.
static constexpr s32 kPreferredSeedOffset = 3;
static constexpr s32 kSeedLength = 5;

Trim3pYs::Trim3pYs(SeqMatcher sid_3p_matcher, const std::string& stem_3p)
    : _sid_3p_matcher(std::move(sid_3p_matcher)),
      _flank_3p_bitap(stem_3p, SearchDirection::kBackward),
      _seed_offset(std::min(kPreferredSeedOffset, std::max(0, static_cast<s32>(stem_3p.size()) - kSeedLength))),
      _stem_seed(stem_3p.substr(_seed_offset, kSeedLength)) {  // NOSONAR: initializer list required
  if (_stem_seed.empty()) {
    throw error::Error("3' stem '{}' is too short to extract a seed", stem_3p);
  }
}

bool Trim3pYs::TryStemAtPosition(const s32 match_pos,
                                 const u8* const seq2,
                                 const u32 length,
                                 const u32 insert_start,
                                 Trim3pInfoYs& trim_info) const {
  assert(match_pos >= 0 && "TryStemAtPosition called with negative match_pos");
  const auto current_pos = static_cast<u32>(match_pos);

  auto sid_match0 = _sid_3p_matcher.FindBarcode(ReadEnd::k3p, current_pos, seq2, length);
  if (AbortTrim(sid_match0, insert_start)) {
    return true;
  }
  if (!sid_match0.IsUnknown()) {
    trim_info.sid_match = sid_match0;
    trim_info.insert_end = sid_match0.SPos();

    return true;
  }
  return false;
}

Trim3pInfoYs Trim3pYs::Trim(const u8* const seq2,
                            const u32 length,
                            const char* const seq,
                            const u32 insert_start) const {
  Trim3pInfoYs trim_info(length);

  if (length == 0) {
    return trim_info;
  }

  const auto seq_view = std::string_view(seq, length);
  const auto min_pos = static_cast<s32>(insert_start);
  const auto query_len = _flank_3p_bitap.GetQueryLength();
  constexpr s32 kMaxDist = 4;

  // Seed-guided left-to-right search.
  //
  // The original algorithm used a single backward Bitap::Find over the last
  // 64bp of the read.  Two problems arise on SBX data:
  //
  //  1. Post-adapter noise can push the 3' stem beyond the 64bp window,
  //     causing under-trimming.
  //  2. Duplicated 3' adapters place two stem copies in the read; the R→L
  //     Bitap finds the rightmost (wrong) copy, so the SID matcher locks
  //     onto the inter-stem SID and leaves the first adapter in the output.
  //
  // Scanning L→R for a short exact seed, then running a targeted Bitap in a
  // small window around each hit, solves both: the leftmost stem is found
  // first, and detection extends beyond the 64bp window.
  const auto seed_len = static_cast<s32>(_stem_seed.size());
  const s32 search_end = static_cast<s32>(length) - seed_len;

  for (s32 i = min_pos; i <= search_end; ++i) {
    if (seq_view.substr(static_cast<size_t>(i), static_cast<size_t>(seed_len)) != _stem_seed) {
      continue;
    }

    const s32 estimated_stem_start = i - _seed_offset;
    const s32 window_begin = std::max(min_pos, estimated_stem_start - kMaxDist);
    const s32 window_end = std::min(static_cast<s32>(length) - 1, estimated_stem_start + query_len + kMaxDist - 1);

    if (window_begin > window_end) {
      continue;
    }

    const s32 match_pos = _flank_3p_bitap.Find(seq_view, window_begin, window_end);
    if (match_pos == kNoMatchPosition) {
      continue;
    }

    if (TryStemAtPosition(match_pos, seq2, length, insert_start, trim_info)) {
      return trim_info;
    }
  }

  // Fallback: original Bitap over the last 64bp.
  //
  // ~7% of real stems have the seed disrupted (ed >= 1 in the seed region).
  // The backward Bitap catches these.  This only runs when the seed scan
  // found nothing, so duplicated-adapter reads (always caught above) are
  // unaffected.  The SID matcher always runs — either from the Bitap-found
  // stem or from the end of the read — mirroring the original algorithm.
  {
    auto fallback_pos = length;
    const s32 search_start = std::max(0, static_cast<s32>(length) - Bitap<4>::kQueryWindowSize);
    const auto end = static_cast<s32>(length - 1);
    const s32 match_pos = _flank_3p_bitap.Find(seq_view, search_start, end);
    if (match_pos != kNoMatchPosition) {
      fallback_pos = static_cast<u32>(match_pos);
    }

    auto sid_match0 = _sid_3p_matcher.FindBarcode(ReadEnd::k3p, fallback_pos, seq2, length);
    if (AbortTrim(sid_match0, insert_start)) {
      return trim_info;
    }
    if (!sid_match0.IsUnknown()) {
      trim_info.sid_match = sid_match0;
      trim_info.insert_end = sid_match0.SPos();
    }
  }

  return trim_info;
}

bool Trim3pYs::AbortTrim(const MatchInfo& match_info, const u32 insert_start) {
  // Abort on an ambiguous match, or a found 3' match overlapping the 5' insert start.
  // An unknown (not-found) match is permitted; completeness is decided downstream.
  return (!match_info.IsUnknown() && match_info.SPos() <= insert_start) || match_info.IsAmbiguous();
}

}  // namespace xoos::demux
