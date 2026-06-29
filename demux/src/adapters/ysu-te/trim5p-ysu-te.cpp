#include "trim5p-ysu-te.h"

#include <xoos/enum/enum-util.h>

#include <algorithm>
#include <string_view>
#include <utility>

#include "utility/string-util.h"

namespace xoos::demux {

Trim5pYsuTe::Trim5pYsuTe(bool enable_partial, SeqMatcher runway_5p_matcher, SeqMatcher sid_5p_matcher,
                         SeqMatcher sid_umi_spacer_5p_matcher, SeqMatcher umi_5p_matcher)
    : _enable_partial{enable_partial},
      _runway_5p_matcher{std::move(runway_5p_matcher)},
      _sid_5p_matcher(std::move(sid_5p_matcher)),
      _sid_umi_spacer_5p_matcher{std::move(sid_umi_spacer_5p_matcher)},
      _umi_5p_matcher(std::move(umi_5p_matcher)) {}

Trim5pInfoYsuTe Trim5pYsuTe::Trim(const uint8_t* seq2, size_t length, const char* seq) const {
  uint current_pos{0};

  Trim5pInfoYsuTe trim_info;

  // search for the bait sequence adjacent to the start of 5'sid and update
  // current position. This would better handle situations where there are more
  // extra bases before the SID than expected.
  const std::string bait_seq_5p{"GAGTCTTTT"};
  const auto search_length = std::min(64ul, length);
  std::string_view head(seq, search_length);
  auto pos = head.find(bait_seq_5p);
  if (pos != std::string_view::npos) {
    current_pos = pos + bait_seq_5p.size();
  }

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
    auto umi_match = _umi_5p_matcher.FindBarcode(ReadEnd::k5p, current_pos, seq2, length);
    if (AbortTrim(umi_match)) {
      return trim_info;
    }

    if (!umi_match.IsUnknown()) {
      current_pos = umi_match.EPos();
      trim_info.sid_match = sid_match;
      trim_info.umi_match = umi_match;
      trim_info.insert_start = current_pos;
      const std::string overhang_seq_5p{"GTT"};
      const int overhang_size = 3;
      auto overhang_match = FindExactMatch5p(overhang_seq_5p.c_str(), overhang_size, 2, current_pos, seq, length);
      if (overhang_match) {
        current_pos = overhang_match->epos;
        trim_info.insert_start = current_pos;
      } else {
        current_pos += std::min(static_cast<uint>(overhang_seq_5p.length()), static_cast<uint>(length) - current_pos);
        trim_info.insert_start = current_pos;
      }
    }
  }

  return trim_info;
}

bool Trim5pYsuTe::AbortTrim(const MatchInfo& match_info) const {
  return (match_info.IsUnknown() && !_enable_partial) || match_info.IsAmbiguous();
}
}  // namespace xoos::demux
