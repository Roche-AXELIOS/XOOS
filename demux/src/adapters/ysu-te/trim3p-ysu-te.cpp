#include "trim3p-ysu-te.h"

#include <xoos/enum/enum-util.h>

#include <algorithm>
#include <string_view>
#include <utility>

#include "utility/string-util.h"

namespace xoos::demux {

Trim3pInfoYsuTe::Trim3pInfoYsuTe(uint insert_end) : insert_end(insert_end), sid_match{}, umi_match{} {}

Trim3pYsuTe::Trim3pYsuTe(bool enable_partial, SeqMatcher sid_3p_matcher, SeqMatcher umi_sid_spacer_3p_matcher,
                         SeqMatcher umi_3p_matcher)
    : _enable_partial(enable_partial),
      _sid_3p_matcher(std::move(sid_3p_matcher)),
      _umi_sid_spacer_3p_matcher(std::move(umi_sid_spacer_3p_matcher)),
      _umi_3p_matcher(std::move(umi_3p_matcher)) {}

Trim3pInfoYsuTe Trim3pYsuTe::Trim(const uint8_t* seq2, size_t length, const char* seq, uint insert_start) const {
  Trim3pInfoYsuTe trim_info(static_cast<uint>(length));
  std::optional<MatchInfo> sid_match;
  uint current_pos = length;

  // search for the bait sequence adjacent to the end of 3'sid and update
  // current position. This would better handle situations where there are more
  // extra bases after the SID than expected.
  const std::string bait_seq_3p{"TTTGTCGT"};
  const auto search_length = std::min(64ul, length);
  std::string_view tail(seq + length - search_length, search_length);
  auto pos = tail.rfind(bait_seq_3p);
  if (pos != std::string::npos) {
    current_pos = length - 1 - search_length + pos;
  }

  {
    auto sid_match0 = _sid_3p_matcher.FindBarcode(ReadEnd::k3p, current_pos, seq2, length);
    if (AbortTrim(sid_match0, insert_start)) {
      return trim_info;
    }

    if (!sid_match0.IsUnknown()) {
      current_pos = sid_match0.SPos();
      sid_match = sid_match0;
      trim_info.sid_match = sid_match;
      trim_info.insert_end = current_pos;
    }
  }

  {
    auto umi_match = _umi_3p_matcher.FindBarcode(ReadEnd::k3p, current_pos, seq2, length);
    if (AbortTrim(umi_match, insert_start)) {
      return trim_info;
    }

    if (!umi_match.IsUnknown()) {
      current_pos = umi_match.SPos();
      trim_info.sid_match = sid_match;
      trim_info.umi_match = umi_match;
      trim_info.insert_end = current_pos;
      const std::string overhang_seq_3p{"AAC"};
      const int overhang_size = 3;
      auto overhang_match = FindExactMatch3p(overhang_seq_3p.c_str(), overhang_size, 2, current_pos, seq, length);
      if (overhang_match) {
        current_pos = overhang_match->spos;
        trim_info.insert_end = current_pos;
      } else {
        current_pos -= std::min(static_cast<uint>(overhang_seq_3p.length()), current_pos);
        trim_info.insert_end = current_pos;
      }
    }
  }

  return trim_info;
}

bool Trim3pYsuTe::AbortTrim(const MatchInfo& match_info, uint insert_start) const {
  return (!match_info.IsUnknown() && match_info.SPos() <= insert_start) ||
         (match_info.IsUnknown() && !_enable_partial) || match_info.IsAmbiguous();
}

}  // namespace xoos::demux
