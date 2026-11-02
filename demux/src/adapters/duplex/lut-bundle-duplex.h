#pragma once

#include <string_view>

#include "lut-bundle/lut-bundle.h"
#include "sequence/fixed-sequence.h"
#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

/**
 * LUT bundle for a duplex adapter. Components are addressed by their position in the adapter read
 * structure, not by a 5'/3' side tag:
 *   [runway_1][sid_2][start_3][loop_4][anchor_5][sid_6]
 * The two SID components differ only in orientation (baked in by the assembler) and are distinguished
 * positionally. The fixed components are held as FixedSequence.
 */
class LutBundleDuplex {
 public:
  LutBundleDuplex(
      FixedSequence runway_1, SeqMatcher sid_2, FixedSequence start_3, FixedSequence loop_4, SeqMatcher sid_6);

  const FixedSequence& Runway1() const;

  const SeqMatcher& Sid2Matcher() const;
  const BarcodePool& Sid2Pool() const;  // NOSONAR(S5536) -- used by FlowManager

  const FixedSequence& StartSequence() const;

  std::string_view LoopSequence() const;

  const SeqMatcher& Sid6Matcher() const;

 private:
  FixedSequence _runway_1;
  SeqMatcher _sid_2_matcher;
  FixedSequence _start_3;
  FixedSequence _loop_4;
  SeqMatcher _sid_6_matcher;
};

template <>
LutBundleDuplex CreateLutBundle(const AdapterDesignBundle& designs);  // NOSONAR(M23_184)

}  // namespace xoos::demux
