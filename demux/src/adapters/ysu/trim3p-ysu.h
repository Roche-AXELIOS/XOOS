#pragma once

#include <optional>
#include <string>

#include <xoos/types/int.h>

#include "sequence/matcher/match-info.h"
#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

/**
 * Information about 3' adapter trimming
 */
struct Trim3pInfoYsu {
  explicit Trim3pInfoYsu(u32 insert_end);

  u32 insert_end;
  std::optional<MatchInfo> sid_match;
  std::optional<MatchInfo> umi_match;
};

/**
 * Responsible for trimming 3' adapters from reads with YSU adapters.
 */
class Trim3pYsu {
 public:
  Trim3pYsu(SeqMatcher sid_3p_matcher, SeqMatcher umi_3p_matcher, std::optional<std::string> bait_3p);

  // This is the main function that should be used to trim 3' adapters. It uses a 2-bit encoding (seq2) to allow
  // for fast LUT operations; it is using the original 8-bit encoding too for an exact string match.
  Trim3pInfoYsu Trim(const u8* seq2, u32 seq_length, const char* seq, u32 insert_start) const;

 private:
  /**
   * Determine if the 3' trimming process should be aborted based on the currently
   * matched sequences in the 3' adapter.
   */
  static bool AbortTrim(const MatchInfo& match_info, u32 insert_start);

  SeqMatcher _sid_3p_matcher;
  SeqMatcher _umi_3p_matcher;
  std::optional<std::string> _bait_3p;
};

}  // namespace xoos::demux
