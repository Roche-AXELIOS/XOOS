#include "trim3p-ys.h"

#include <xoos/enum/enum-util.h>

#include <algorithm>
#include <string_view>
#include <utility>

#include "utility/string-util.h"

namespace xoos::demux {
Trim3pInfoYs::Trim3pInfoYs(const u32 insert_end) : insert_end(insert_end), sid_match{} {}

/**
 * Construct a 3' trimmer for YS adapters.
 *
 * @param enable_partial           Whether trimming is allowed when the SID match is unknown or partial.
 * @param sid_3p_matcher           Matcher used to locate the 3' end of the sample identifier (SID).
 * @param sid_spacer_3p_matcher    Matcher used to locate the 3' boundary between the SID and spacer sequence.
 */
Trim3pYs::Trim3pYs(const bool enable_partial, SeqMatcher sid_3p_matcher, SeqMatcher sid_spacer_3p_matcher)
    : _enable_partial(enable_partial),
      _sid_3p_matcher(std::move(sid_3p_matcher)),
      _sid_spacer_3p_matcher(std::move(sid_spacer_3p_matcher)),
      // TODO: get from adapter bundle instead of hardcoding
      _flank_3p_bitap("TTTGTCGTGTAGG", true) {}

Trim3pInfoYs Trim3pYs::Trim(const u8* const seq2, const u32 length, const char* const seq,
                            const u32 insert_start) const {
  Trim3pInfoYs trim_info(length);
  auto current_pos = length;

  if (length > 0) {
    const s32 search_start = std::max(0, static_cast<s32>(length) - Bitap<4>::kQueryWindowSize);
    // search is inclusive thus -1
    const s32 end = static_cast<s32>(length - 1);
    const s32 match_pos = _flank_3p_bitap.Find(seq, search_start, end);
    if (match_pos != -1) {
      // Reverse Bitap returns match start position for 3' end; use directly without offset.
      current_pos = static_cast<u32>(match_pos);
    }
  }
  {
    auto sid_match0 = _sid_3p_matcher.FindBarcode(ReadEnd::k3p, current_pos, seq2, length);
    if (AbortTrim(sid_match0, insert_start)) {
      return trim_info;
    }
    if (!sid_match0.IsUnknown()) {
      current_pos = sid_match0.SPos();
      trim_info.sid_match = sid_match0;
      trim_info.insert_end = current_pos;
    }
  }
  {
    // TODO: Try to get this from adapter design bundle
    const std::string overhang_seq_3p{"CGC"};
    constexpr u32 kOverhangSize = 3;
    auto overhang_match = FindExactMatch3p(overhang_seq_3p.c_str(), kOverhangSize, 2, current_pos, seq, length);
    if (overhang_match) {
      current_pos = overhang_match->spos;
      trim_info.insert_end = current_pos;
    } else {
      current_pos -= std::min(static_cast<u32>(overhang_seq_3p.length()), current_pos);
      trim_info.insert_end = current_pos;
    }
  }
  return trim_info;
}

bool Trim3pYs::AbortTrim(const MatchInfo& match_info, const u32 insert_start) const {
  return (!match_info.IsUnknown() && match_info.SPos() <= insert_start) ||
         (match_info.IsUnknown() && !_enable_partial) || match_info.IsAmbiguous();
}

}  // namespace xoos::demux
