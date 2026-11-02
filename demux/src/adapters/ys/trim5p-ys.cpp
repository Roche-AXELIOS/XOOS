#include "trim5p-ys.h"

#include <algorithm>
#include <string>
#include <string_view>
#include <utility>

#include <xoos/enum/enum-util.h>
#include <xoos/types/int.h>

#include "sequence/matcher/bitap.h"
#include "sequence/sequence-enums.h"

namespace xoos::demux {

Trim5pYs::Trim5pYs(const FixedSequence& runway_5p, SeqMatcher sid_5p_matcher)
    : _sid_5p_matcher(std::move(sid_5p_matcher)),
      _runway_5p_sequence(runway_5p.Seq()),
      _runway_5p_bitap(_runway_5p_sequence, SearchDirection::kForward) {
}

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
  u32 current_pos{0};

  Trim5pInfoYs trim_info;

  if (length > 0) {
    const s32 end = std::min(Bitap<4>::kQueryWindowSize - 1, static_cast<s32>(length - 1));
    const s32 match_pos = _runway_5p_bitap.Find(seq, 0, end);
    if (match_pos != -1) {
      // Forward Bitap returns the match end position; +1 gives the position just after the
      // runway, which is where the SID begins.
      current_pos = static_cast<u32>(match_pos + 1);
    }
  }

  // Find the sid
  auto sid_match0 = _sid_5p_matcher.FindBarcode(ReadEnd::k5p, current_pos, seq2, length);
  if (AbortTrim(sid_match0)) {
    return trim_info;
  }

  if (!sid_match0.IsUnknown()) {
    current_pos = sid_match0.EPos();
    const std::optional sid_match = sid_match0;
    trim_info.sid_match = sid_match;
    trim_info.insert_start = current_pos;
  }
  return trim_info;
}

bool Trim5pYs::AbortTrim(const MatchInfo& match_info) {
  // Abort only on an ambiguous match. An unknown (not-found) match is permitted: the read is
  // trimmed to whatever landmarks were found, and completeness/assignment are decided downstream.
  return match_info.IsAmbiguous();
}
}  // namespace xoos::demux
