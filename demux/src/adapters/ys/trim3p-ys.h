#pragma once

#include <string>

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
  /**
   * @param sid_3p_matcher  Barcode matcher for the 3' SID sequence.
   * @param stem_3p         3' stem sequence from the adapter bundle, used to build the Bitap fuzzy matcher.
   */
  Trim3pYs(SeqMatcher sid_3p_matcher, const std::string& stem_3p);

  // This is the main function that should be used to trim 3' adapters. It uses a 2-bit encoding (seq2) to allow
  // for fast LUT operations; it is using the original 8-bit encoding too for an exact string match.
  Trim3pInfoYs Trim(const u8* seq2, u32 length, const char* seq, u32 insert_start) const;

 private:
  /**
   * Try to match the stem+SID at a Bitap-confirmed position.
   * Returns true if a definitive result was reached (match found or trim aborted).
   */
  bool TryStemAtPosition(s32 match_pos, const u8* seq2, u32 length, u32 insert_start, Trim3pInfoYs& trim_info) const;

  /**
   * Determine if the 3' trimming process should be aborted based on the currently
   * matched sequences in the 3' adapter.
   */
  static bool AbortTrim(const MatchInfo& match_info, u32 insert_start);

  SeqMatcher _sid_3p_matcher;
  // Bitap matcher for 3' flank sequence with up to 4 edit distance
  const Bitap<4> _flank_3p_bitap;
  // Offset of _stem_seed within the stem (clamped if the stem is shorter than the preferred offset)
  const s32 _seed_offset;
  // Exact seed extracted from the stem for the L→R scan
  const std::string _stem_seed;
};

}  // namespace xoos::demux
