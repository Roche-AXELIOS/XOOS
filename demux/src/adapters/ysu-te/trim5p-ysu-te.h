#pragma once

#include "sequence/matcher/match-info.h"
#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {
/**
 * Information about 5' adapter trimming
 */
struct Trim5pInfoYsuTe {
  uint insert_start;
  std::optional<MatchInfo> sid_match;
  std::optional<MatchInfo> umi_match;

  Trim5pInfoYsuTe() : insert_start{0}, sid_match{}, umi_match{} {}
};

/**
 * Responsible for trimming 5' adapters from reads with YSU-TE adapters.
 */
class Trim5pYsuTe {
 public:
  Trim5pYsuTe(bool enable_partial, SeqMatcher runway_5p_matcher, SeqMatcher sid_5p_matcher,
              SeqMatcher sid_umi_spacer_5p_matcher, SeqMatcher umi_5p_matcher);

  // This is the main function that should be used to trim 5' adapters. It uses a 2-bit encoding (seq2) to allow
  // for fast LUT operations; it is using the original 8-bit encoding too for an exact string match.
  Trim5pInfoYsuTe Trim(const uint8_t* seq2, size_t length, const char* seq) const;

 private:
  bool AbortTrim(const MatchInfo& match_info) const;

 private:
  bool _enable_partial;
  SeqMatcher _runway_5p_matcher;
  SeqMatcher _sid_5p_matcher;
  SeqMatcher _sid_umi_spacer_5p_matcher;
  SeqMatcher _umi_5p_matcher;
};

}  // namespace xoos::demux
