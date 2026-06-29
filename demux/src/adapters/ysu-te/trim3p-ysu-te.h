#pragma once

#include "sequence/matcher/match-info.h"
#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

/**
 * Information about 3' adapter trimming
 */
struct Trim3pInfoYsuTe {
  explicit Trim3pInfoYsuTe(uint insert_end);

  uint insert_end;
  std::optional<MatchInfo> sid_match;
  std::optional<MatchInfo> umi_match;
};

/**
 * Responsible for trimming 3' adapters from reads with YSU-TE adapters.
 */
class Trim3pYsuTe {
 public:
  Trim3pYsuTe(bool enable_partial, SeqMatcher sid_3p_matcher, SeqMatcher umi_sid_spacer_3p_matcher,
              SeqMatcher umi_3p_matcher);

  // This is the main function that should be used to trim 3' adapters. It uses a 2-bit encoding (seq2) to allow
  // for fast LUT operations; it is using the original 8-bit encoding too for an exact string match.
  Trim3pInfoYsuTe Trim(const uint8_t* seq2, size_t seq_length, const char* seq, uint insert_start) const;

 private:
  /**
   * Determine if the 3' trimming process should be aborted based on the currently
   * matched sequences in the 3' adapter.
   */
  bool AbortTrim(const MatchInfo& match_info, uint insert_start) const;

  bool _enable_partial;
  SeqMatcher _sid_3p_matcher;
  SeqMatcher _umi_sid_spacer_3p_matcher;
  SeqMatcher _umi_3p_matcher;
};

}  // namespace xoos::demux
