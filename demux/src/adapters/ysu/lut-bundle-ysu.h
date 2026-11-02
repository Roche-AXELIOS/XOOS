#pragma once

#include <optional>
#include <string>

#include "lut-bundle/lut-bundle.h"
#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

/**
 * A bundle of LUTs to be used in trimming YSU data. These LUTs are meant to be loaded from disk, but
 * they might be created in memory for test cases.
 *
 * Components are addressed by their position in the adapter read structure, not by a 5'/3' side tag:
 *   [sid_1][umi_2][insert][umi_3][sid_4]
 * Each SID may carry an Infer'd bait seed (the adjacent fixed trimmed to the configured bait length),
 * held alongside its matcher.
 */
class LutBundleYsu {
 public:
  LutBundleYsu(SeqMatcher sid_1,
               SeqMatcher umi_2,
               SeqMatcher umi_3,
               SeqMatcher sid_4,
               std::optional<std::string> bait_1,
               std::optional<std::string> bait_4);

  const SeqMatcher& Sid1Matcher() const;
  const BarcodePool& Sid1Pool() const;  // NOSONAR(S5536) -- used by FlowManager

  const SeqMatcher& Umi2Matcher() const;

  const SeqMatcher& Umi3Matcher() const;

  const SeqMatcher& Sid4Matcher() const;

  const std::optional<std::string>& Bait1() const;
  const std::optional<std::string>& Bait4() const;

 private:
  SeqMatcher _sid_1_matcher;
  SeqMatcher _umi_2_matcher;
  SeqMatcher _umi_3_matcher;
  SeqMatcher _sid_4_matcher;

  std::optional<std::string> _bait_1;
  std::optional<std::string> _bait_4;
};

template <>
LutBundleYsu CreateLutBundle(const AdapterDesignBundle& designs);  // NOSONAR(M23_184)

}  // namespace xoos::demux
