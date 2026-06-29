#pragma once

#include "sequence/matcher/bitap.h"
#include "sequence/matcher/match-info.h"
#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

/**
 * Information about 3' adapter trimming
 */
struct Trim3pInfoYs {
  explicit Trim3pInfoYs(u32 insert_end);

  u32 insert_end;
  std::optional<MatchInfo> sid_match;
};

/**
 * Responsible for trimming 3' adapters from reads with YS adapters.
 */
class Trim3pYs {
 public:
  Trim3pYs(bool enable_partial, SeqMatcher sid_3p_matcher, SeqMatcher sid_spacer_3p_matcher);

  // This is the main function that should be used to trim 3' adapters. It uses a 2-bit encoding (seq2) to allow
  // for fast LUT operations; it is using the original 8-bit encoding too for an exact string match.
  Trim3pInfoYs Trim(const u8* seq2, u32 length, const char* seq, u32 insert_start) const;

 private:
  /**
   * Determine if the 3' trimming process should be aborted based on the currently
   * matched sequences in the 3' adapter.
   */
  bool AbortTrim(const MatchInfo& match_info, u32 insert_start) const;

  bool _enable_partial;
  SeqMatcher _sid_3p_matcher;
  SeqMatcher _sid_spacer_3p_matcher;
  // Bitap matcher for 3' flank sequence with up to 4 edit distance
  const Bitap<4> _flank_3p_bitap;
};

}  // namespace xoos::demux
