#pragma once

#include "adapters/duplex/stem/lut-bundle-duplex-stem.h"
#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

class LutBundleDuplexUmi : public LutBundleDuplexStem {
 public:
  LutBundleDuplexUmi(SeqMatcher sid_5p, SeqMatcher start_sequence, SeqMatcher sid_3p, SeqMatcher stop_sequence,
                     SeqMatcher stem_5p, SeqMatcher stem_3p, SeqMatcher umi_5p, SeqMatcher umi_3p,
                     SeqMatcher loop_sequence_matcher);

  // UMI accessors
  const SeqMatcher& Umi5pMatcher() const;
  const BarcodePool& Umi5pPool() const;

  const SeqMatcher& Umi3pMatcher() const;
  const BarcodePool& Umi3pPool() const;

 private:
  SeqMatcher _umi_5p_matcher;
  SeqMatcher _umi_3p_matcher;
};

}  // namespace xoos::demux
