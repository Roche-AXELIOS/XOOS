#pragma once

#include "sequence/matcher/bitap.h"
#include "sequence/matcher/match-info.h"
#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {
/**
 * Information about 5' adapter trimming
 */
struct Trim5pInfoYs {
  u32 insert_start;
  std::optional<MatchInfo> sid_match;

  Trim5pInfoYs() : insert_start{0}, sid_match{} {}
};

/**
 * Responsible for trimming 5' adapters from reads with YS adapters.
 */
class Trim5pYs {
 public:
  Trim5pYs(bool enable_partial, SeqMatcher runway_5p_matcher, SeqMatcher sid_5p_matcher,
           SeqMatcher sid_spacer_5p_matcher);

  Trim5pInfoYs Trim(const u8* seq2, size_t length, const char* seq) const;

 private:
  bool AbortTrim(const MatchInfo& match_info) const;
  bool _enable_partial;
  SeqMatcher _runway_5p_matcher;
  SeqMatcher _sid_5p_matcher;
  SeqMatcher _sid_spacer_5p_matcher;
  // Bitap matcher for 5' flank sequence with up to 4 edit distance
  const Bitap<4> _flank_5p_bitap;
};

}  // namespace xoos::demux
