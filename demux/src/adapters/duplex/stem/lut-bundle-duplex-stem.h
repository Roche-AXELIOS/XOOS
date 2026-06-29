#pragma once

#include "adapters/duplex/lut-bundle-duplex.h"
#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

class LutBundleDuplexStem : public LutBundleDuplex {
 public:
  LutBundleDuplexStem(SeqMatcher sid_5p, SeqMatcher start_sequence, SeqMatcher sid_3p, SeqMatcher stop_sequence,
                      SeqMatcher stem_5p, SeqMatcher stem_3p, SeqMatcher loop_sequence_matcher);

  // stem accessors
  const SeqMatcher& Stem5pMatcher() const;
  const SeqMatcher& Stem3pMatcher() const;

 private:
  SeqMatcher _stem_5p_matcher;
  SeqMatcher _stem_3p_matcher;
};

}  // namespace xoos::demux
