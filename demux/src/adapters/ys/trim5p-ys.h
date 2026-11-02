#pragma once

#include "sequence/fixed-sequence.h"
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

  Trim5pInfoYs() : insert_start{0} {
  }
};

/**
 * Responsible for trimming 5' adapters from reads with YS adapters.
 */
class Trim5pYs {
 public:
  Trim5pYs(const FixedSequence& runway_5p, SeqMatcher sid_5p_matcher);

  Trim5pInfoYs Trim(const u8* seq2, size_t length, const char* seq) const;

 private:
  static bool AbortTrim(const MatchInfo& match_info);
  SeqMatcher _sid_5p_matcher;
  // 5' runway sequence: the fixed component immediately preceding the SID.
  std::string _runway_5p_sequence;
  // Bitap matcher for the 5' runway with up to 4 edit distance. A fuzzy match anchors the
  // position just after the runway (i.e. the SID start) before the SID matcher runs.
  const Bitap<4> _runway_5p_bitap;
};

}  // namespace xoos::demux
