#include "lut-bundle-ys.h"

#include <xoos/error/error.h>

#include "adapter-design/adapter-component-indices.h"
#include "adapter-design/adapter-design-bundle.h"
// Clangd's include-cleaner doesn't detect usage from template specializations.
// This include is required for explicit specialization below.
#include "lut-bundle/lut-bundle.h"  // IWYU pragma: keep

namespace xoos::demux {

LutBundleYs::LutBundleYs(FixedSequence runway_1, SeqMatcher sid_2, SeqMatcher sid_4, std::string bait_5)
    : _runway_1{std::move(runway_1)},
      _sid_2_matcher{std::move(sid_2)},
      _sid_4_matcher{std::move(sid_4)},
      _bait_5{std::move(bait_5)} {
}

const FixedSequence& LutBundleYs::Runway1() const {
  return _runway_1;
}

const SeqMatcher& LutBundleYs::Sid2Matcher() const {
  return _sid_2_matcher;
}

const BarcodePool& LutBundleYs::Sid2Pool() const {
  return _sid_2_matcher.Pool();
}

const SeqMatcher& LutBundleYs::Sid4Matcher() const {
  return _sid_4_matcher;
}

const std::string& LutBundleYs::Bait5() const {
  return _bait_5;
}

template <>
LutBundleYs CreateLutBundle(const AdapterDesignBundle& designs) {
  // Read-order component layout: [runway, sid, spacer, sid, bait]. Components are addressed by
  // position; the trailing bait is just the last fixed literal.
  using Index = adapter_component_index::Ys;
  const std::string bait{designs.FixedAt(Index::kBait).Seq()};
  if (bait.empty()) {
    throw error::Error("Adapter design '{}': bait must not be empty", designs.design.name);
  }
  return LutBundleYs{
      designs.FixedAt(Index::kRunway), *designs.MatcherAt(Index::kSidFwd), *designs.MatcherAt(Index::kSidRev), bait};
}

}  // namespace xoos::demux
