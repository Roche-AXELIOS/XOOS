#pragma once

#include "adapters/duplex/stem/lut-bundle-duplex-stem.h"
#include "sequence/fixed-sequence.h"
#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

/**
 * LUT bundle for a stem-duplex adapter carrying UMIs. Components are addressed by their read-structure
 * position:
 *   [runway_1][sid_2][start_3][stem_4][loop_5][umi_6][sid_7][anchor_8][stem_9][umi_10]
 * The stem-duplex components are inherited; this adds the two UMI matchers, which differ only in
 * orientation and are distinguished positionally.
 */
class LutBundleDuplexUmi : public LutBundleDuplexStem {
 public:
  LutBundleDuplexUmi(FixedSequence runway_1,
                     SeqMatcher sid_2,
                     FixedSequence start_3,
                     FixedSequence stem_4,
                     FixedSequence loop_5,
                     SeqMatcher umi_6,
                     SeqMatcher sid_7,
                     const FixedSequence& stem_9,
                     SeqMatcher umi_10);

  // UMI accessors
  const SeqMatcher& Umi6Matcher() const;

  const SeqMatcher& Umi10Matcher() const;

 private:
  SeqMatcher _umi_6_matcher;
  SeqMatcher _umi_10_matcher;
};

template <>
LutBundleDuplexUmi CreateLutBundle(const AdapterDesignBundle& designs);  // NOSONAR(M23_184)

}  // namespace xoos::demux
