#pragma once

#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

namespace fs = std::filesystem;

/**
 * A bundle of LUT to be used in trimming YSU-TE data. These LUTs are meant
 * to be loaded from disk, but they might be created in memory for tests cases.
 */
class LutBundleYsuTe {
 public:
  LutBundleYsuTe(SeqMatcher runway_5p, SeqMatcher sid_5p, SeqMatcher sid_umi_spacer_5p, SeqMatcher umi_5p,
                 SeqMatcher sid_3p, SeqMatcher umi_sid_spacer_3p, SeqMatcher umi_3p);

  const SeqMatcher& Runway5pMatcher() const;

  const SeqMatcher& Sid5pMatcher() const;
  const BarcodePool& Sid5pPool() const;

  const SeqMatcher& SidUmiSpacer5pMatcher() const;
  const BarcodePool& SidUmiSpacer5pPool() const;

  const SeqMatcher& Umi5pMatcher() const;
  const BarcodePool& Umi5pPool() const;

  const SeqMatcher& Sid3pMatcher() const;
  const BarcodePool& Sid3pPool() const;

  const SeqMatcher& SidUmiSpacer3pMatcher() const;
  const BarcodePool& SidUmiSpacer3pPool() const;

  const SeqMatcher& Umi3pMatcher() const;
  const BarcodePool& Umi3pPool() const;

 private:
  SeqMatcher _runway_5p_matcher;
  SeqMatcher _sid_5p_matcher;
  SeqMatcher _sid_umi_spacer_5p_matcher;
  SeqMatcher _umi_5p_matcher;

  SeqMatcher _sid_3p_matcher;
  SeqMatcher _umi_sid_spacer_3p_matcher;
  SeqMatcher _umi_3p_matcher;
};
}  // namespace xoos::demux
