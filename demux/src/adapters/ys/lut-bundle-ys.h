#pragma once

#include <string>

#include "lut-bundle/lut-bundle.h"
#include "sequence/fixed-sequence.h"
#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

/**
 * A bundle of LUTs to be used in trimming YS data. These LUTs are meant to be loaded from disk, but
 * they might be created in memory for test cases.
 *
 * Components are addressed by their position in the adapter read structure, not by a 5'/3' side tag:
 *   [runway_1][sid_2][spacer][insert][spacer'][sid_4][bait_5]
 * The fixed components (runway, bait) are held as FixedSequence -- no LUT is built for them since
 * consumers only read the raw sequence. The variable SID components remain SeqMatcher.
 */
class LutBundleYs {
 public:
  LutBundleYs(FixedSequence runway_1, SeqMatcher sid_2, SeqMatcher sid_4, std::string bait_5);

  const FixedSequence& Runway1() const;

  const SeqMatcher& Sid2Matcher() const;
  const BarcodePool& Sid2Pool() const;  // NOSONAR(S5536) -- used by FlowManager

  const SeqMatcher& Sid4Matcher() const;

  const std::string& Bait5() const;

 private:
  FixedSequence _runway_1;
  SeqMatcher _sid_2_matcher;
  SeqMatcher _sid_4_matcher;
  std::string _bait_5;
};

template <>
LutBundleYs CreateLutBundle(const AdapterDesignBundle& designs);  // NOSONAR(M23_184)

}  // namespace xoos::demux
