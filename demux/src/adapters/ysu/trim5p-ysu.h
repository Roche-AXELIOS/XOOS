#pragma once

#include <optional>
#include <string>

#include <xoos/types/int.h>

#include "sequence/matcher/match-info.h"
#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {
/**
 * Information about 5' adapter trimming
 */
struct Trim5pInfoYsu {
  u32 insert_start{0};
  std::optional<MatchInfo> sid_match;
  std::optional<MatchInfo> umi_match;
};

/**
 * Responsible for trimming 5' adapters from reads with YSU adapters.
 */
class Trim5pYsu {
 public:
  Trim5pYsu(SeqMatcher sid_5p_matcher, SeqMatcher umi_5p_matcher, std::optional<std::string> bait_5p);

  // This is the main function that should be used to trim 5' adapters. It uses a 2-bit encoding (seq2) to allow
  // for fast LUT operations; it is using the original 8-bit encoding too for an exact string match.
  Trim5pInfoYsu Trim(const u8* seq2, u32 length, const char* seq) const;

 private:
  static bool AbortTrim(const MatchInfo& match_info);

  SeqMatcher _sid_5p_matcher;
  SeqMatcher _umi_5p_matcher;
  std::optional<std::string> _bait_5p;
};

}  // namespace xoos::demux
