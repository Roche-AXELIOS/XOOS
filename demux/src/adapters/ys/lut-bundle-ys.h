#pragma once

#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

namespace fs = std::filesystem;

/**
 * A bundle of LUT to be used in trimming YS data. These LUTs are meant
 * to be loaded from disk, but they might be created in memory for tests cases.
 */
class LutBundleYs {
 public:
  LutBundleYs(SeqMatcher runway_5p, SeqMatcher sid_5p, SeqMatcher sid_spacer_5p, SeqMatcher sid_3p,
              SeqMatcher sid_spacer_3p);

  const SeqMatcher& Runway5pMatcher() const;

  const SeqMatcher& Sid5pMatcher() const;
  const BarcodePool& Sid5pPool() const;

  const SeqMatcher& SidSpacer5pMatcher() const;
  const BarcodePool& SidSpacer5pPool() const;

  const SeqMatcher& Sid3pMatcher() const;
  const BarcodePool& Sid3pPool() const;

  const SeqMatcher& SidSpacer3pMatcher() const;
  const BarcodePool& SidSpacer3pPool() const;

 private:
  SeqMatcher _runway_5p_matcher;
  SeqMatcher _sid_5p_matcher;
  SeqMatcher _sid_spacer_5p_matcher;

  SeqMatcher _sid_3p_matcher;
  SeqMatcher _sid_spacer_3p_matcher;
};
}  // namespace xoos::demux
