#include "trim5p-ysu.h"

#include <algorithm>
#include <string_view>
#include <utility>

#include <xoos/enum/enum-util.h>
#include <xoos/types/int.h>

#include "sequence/sequence-enums.h"

namespace xoos::demux {

Trim5pYsu::Trim5pYsu(SeqMatcher sid_5p_matcher, SeqMatcher umi_5p_matcher, std::optional<std::string> bait_5p)
    : _sid_5p_matcher(std::move(sid_5p_matcher)),
      _umi_5p_matcher(std::move(umi_5p_matcher)),
      _bait_5p{std::move(bait_5p)} {
}

Trim5pInfoYsu Trim5pYsu::Trim(const u8* const seq2, const u32 length, const char* const seq) const {
  u32 current_pos{0};

  Trim5pInfoYsu trim_info;

  // Search for the bait sequence adjacent to the start of 5'sid and update
  // current position. This would better handle situations where there are more
  // extra bases before the SID than expected.
  if (_bait_5p.has_value()) {
    const u32 search_length = std::min(64U, length);
    const std::string_view head(seq, search_length);
    const auto pos = head.find(_bait_5p.value());
    if (pos != std::string_view::npos) {
      current_pos = static_cast<u32>(pos) + static_cast<u32>(_bait_5p->size());
    }
  }

  {
    auto sid_match0 = _sid_5p_matcher.FindBarcode(ReadEnd::k5p, current_pos, seq2, length);
    if (AbortTrim(sid_match0)) {
      return trim_info;
    }

    // Report the SID independently of the UMI. When only the SID is found, the read is still
    // trimmed to the SID edge ("no-star" trimming) — the UMI width is not imputed. Assignment
    // is SID-only; a missing SID+UMI combination is handled by the post-assignment UMI gate.
    if (!sid_match0.IsUnknown()) {
      current_pos = sid_match0.EPos();
      trim_info.sid_match = sid_match0;
      trim_info.insert_start = current_pos;
    }
  }

  {
    auto umi_match = _umi_5p_matcher.FindBarcode(ReadEnd::k5p, current_pos, seq2, length);
    if (AbortTrim(umi_match)) {
      return trim_info;
    }

    // Report the UMI independently. When found, extend the trim past the UMI edge.
    if (!umi_match.IsUnknown()) {
      current_pos = umi_match.EPos();
      trim_info.umi_match = umi_match;
      trim_info.insert_start = current_pos;
    }
  }

  return trim_info;
}

bool Trim5pYsu::AbortTrim(const MatchInfo& match_info) {
  // Abort only on an ambiguous match. An unknown (not-found) match is permitted: the read is
  // trimmed to whatever landmarks were found, and completeness/assignment are decided downstream.
  return match_info.IsAmbiguous();
}
}  // namespace xoos::demux
