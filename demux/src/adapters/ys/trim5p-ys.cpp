#include "trim5p-ys.h"

#include <xoos/enum/enum-util.h>

#include <algorithm>
#include <string_view>
#include <utility>

#include "sequence/matcher/bitap.h"
#include "utility/string-util.h"

namespace xoos::demux {

Trim5pYs::Trim5pYs(bool enable_partial, SeqMatcher runway_5p_matcher, SeqMatcher sid_5p_matcher,
                   SeqMatcher sid_spacer_5p_matcher)
    : _enable_partial{enable_partial},
      _runway_5p_matcher{std::move(runway_5p_matcher)},
      _sid_5p_matcher(std::move(sid_5p_matcher)),
      _sid_spacer_5p_matcher{std::move(sid_spacer_5p_matcher)},
      // TODO: get from adapter bundle instead of hardcoding
      _flank_5p_bitap("CAACAAGAGTCTTTT", false) {}

/**
 * This is the main function that should be used to trim 5' adapters. It uses a 2-bit encoding (seq2) to allow
 * for fast LUT operations; it is using the original 8-bit encoding too for an exact string match.
 *
 * TODO: Why are we allowing for so much hardcoded strings that aren't using the bundle
 *       Alternatively, we could hardcode (if the optimizer can use this faster) but this should be via code generation
 *
 * @param seq2 2-bit encoded sequence for fast LUT operations
 * @param length length of the sequence
 * @param seq 8-bit encoded original sequence for exact string matching
 * @return Trim5pInfoYs structure containing trim information including SID match and insert start position
 */
Trim5pInfoYs Trim5pYs::Trim(const u8* seq2, size_t length, const char* seq) const {
  uint current_pos{0};

  Trim5pInfoYs trim_info;

  if (length > 0) {
    const s32 end = std::min(Bitap<4>::kQueryWindowSize - 1, static_cast<s32>(length - 1));
    const s32 match_pos = _flank_5p_bitap.Find(seq, 0, end);
    if (match_pos != -1) {
      // Forward Bitap returns match end position; +1 gives position after flank.
      current_pos = static_cast<u32>(match_pos + 1);
    }
  }

  // Find the sid
  std::optional<MatchInfo> sid_match;
  {
    auto sid_match0 = _sid_5p_matcher.FindBarcode(ReadEnd::k5p, current_pos, seq2, length);
    if (AbortTrim(sid_match0)) {
      return trim_info;
    }

    if (!sid_match0.IsUnknown()) {
      current_pos = sid_match0.EPos();
      sid_match = sid_match0;
      trim_info.sid_match = sid_match;
      trim_info.insert_start = current_pos;
    }
  }

  {
    // handle overhang
    const std::string overhang_seq_5p{"GCGT"};
    const int overhang_size = 4;
    auto overhang_match = FindExactMatch5p(overhang_seq_5p.c_str(), overhang_size, 2, current_pos, seq, length);
    if (overhang_match) {
      current_pos = overhang_match->epos;
      trim_info.insert_start = current_pos;
    } else {
      current_pos += std::min(static_cast<u32>(overhang_seq_5p.length()), static_cast<u32>(length) - current_pos);
      trim_info.insert_start = current_pos;
    }
  }

  return trim_info;
}

bool Trim5pYs::AbortTrim(const MatchInfo& match_info) const {
  return (match_info.IsUnknown() && !_enable_partial) || match_info.IsAmbiguous();
}
}  // namespace xoos::demux
