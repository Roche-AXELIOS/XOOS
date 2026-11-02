#pragma once

#include "adapters/duplex/lut-bundle-duplex.h"
#include "sequence/fixed-sequence.h"
#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

/**
 * LUT bundle for a stem-duplex adapter. Components are addressed by their read-structure position:
 *   [runway_1][sid_2][start_3][stem_4][loop_5][sid_6][anchor_7][stem_8]
 * The base duplex components are inherited; this retains the pre-loop stem used by the trimmer. The
 * post-loop stem remains a positional constructor argument for derived UMI bundles but is not exposed
 * because this adapter never consumes it.
 */
class LutBundleDuplexStem : public LutBundleDuplex {
 public:
  LutBundleDuplexStem(FixedSequence runway_1,
                      SeqMatcher sid_2,
                      FixedSequence start_3,
                      FixedSequence stem_4,
                      FixedSequence loop_5,
                      SeqMatcher sid_6,
                      [[maybe_unused]] const FixedSequence& stem_8);

  const FixedSequence& Stem4() const;

 private:
  FixedSequence _stem_4;
};

template <>
LutBundleDuplexStem CreateLutBundle(const AdapterDesignBundle& designs);  // NOSONAR(M23_184)

}  // namespace xoos::demux
