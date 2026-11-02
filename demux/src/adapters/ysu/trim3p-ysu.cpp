#include "trim3p-ysu.h"

#include <algorithm>
#include <string_view>
#include <utility>

#include <xoos/enum/enum-util.h>
#include <xoos/types/int.h>

#include "sequence/sequence-enums.h"

namespace xoos::demux {

Trim3pInfoYsu::Trim3pInfoYsu(const u32 insert_end) : insert_end(insert_end) {
}

Trim3pYsu::Trim3pYsu(SeqMatcher sid_3p_matcher, SeqMatcher umi_3p_matcher, std::optional<std::string> bait_3p)
    : _sid_3p_matcher(std::move(sid_3p_matcher)),
      _umi_3p_matcher(std::move(umi_3p_matcher)),
      _bait_3p{std::move(bait_3p)} {
}

Trim3pInfoYsu Trim3pYsu::Trim(const u8* const seq2,
                              const u32 seq_length,
                              const char* const seq,
                              const u32 insert_start) const {
  Trim3pInfoYsu trim_info(seq_length);
  u32 current_pos = seq_length;

  // Search for the bait sequence adjacent to the end of 3'sid and update
  // current position. This would better handle situations where there are more
  // extra bases after the SID than expected.
  if (_bait_3p.has_value()) {
    const u32 search_length = std::min(64U, seq_length);
    const std::string_view tail(seq + seq_length - search_length, search_length);
    const auto pos = tail.rfind(_bait_3p.value());
    if (pos != std::string_view::npos) {
      current_pos = seq_length - search_length + static_cast<u32>(pos);
    }
  }

  {
    auto sid_match0 = _sid_3p_matcher.FindBarcode(ReadEnd::k3p, current_pos, seq2, seq_length);
    if (AbortTrim(sid_match0, insert_start)) {
      return trim_info;
    }

    // Report the SID independently of the UMI. When only the SID is found, the read is still
    // trimmed to the SID edge ("no-star" trimming) — the UMI width is not imputed. Assignment
    // is SID-only; a missing SID+UMI combination is handled by the post-assignment UMI gate.
    if (!sid_match0.IsUnknown()) {
      current_pos = sid_match0.SPos();
      trim_info.sid_match = sid_match0;
      trim_info.insert_end = current_pos;
    }
  }

  {
    auto umi_match = _umi_3p_matcher.FindBarcode(ReadEnd::k3p, current_pos, seq2, seq_length);
    if (AbortTrim(umi_match, insert_start)) {
      return trim_info;
    }

    // Report the UMI independently. When found, extend the trim past the UMI edge.
    if (!umi_match.IsUnknown()) {
      current_pos = umi_match.SPos();
      trim_info.umi_match = umi_match;
      trim_info.insert_end = current_pos;
    }
  }

  return trim_info;
}

bool Trim3pYsu::AbortTrim(const MatchInfo& match_info, const u32 insert_start) {
  // Abort when a found 3' match would overlap the 5' insert start, or on an ambiguous match.
  // An unknown (not-found) match is permitted: the read is trimmed to whatever landmarks were
  // found, and completeness/assignment are decided downstream.
  return (!match_info.IsUnknown() && match_info.SPos() <= insert_start) || match_info.IsAmbiguous();
}

}  // namespace xoos::demux
