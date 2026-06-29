#include "lut-bundle-ys.h"

#include <utility>

#include "adapter-design/adapter-design-bundle.h"
// Clangd's include-cleaner doesn't detect usage from template specializations.
// This include is required for explicit specialization below.
#include "lut-bundle/lut-bundle.h"  // IWYU pragma: keep

namespace xoos::demux {

LutBundleYs::LutBundleYs(SeqMatcher runway_5p, SeqMatcher sid_5p, SeqMatcher sid_spacer_5p, SeqMatcher sid_3p,
                         SeqMatcher sid_spacer_3p)
    : _runway_5p_matcher{std::move(runway_5p)},
      _sid_5p_matcher{std::move(sid_5p)},
      _sid_spacer_5p_matcher{std::move(sid_spacer_5p)},
      _sid_3p_matcher{std::move(sid_3p)},
      _sid_spacer_3p_matcher{std::move(sid_spacer_3p)} {}

const SeqMatcher& LutBundleYs::Runway5pMatcher() const { return _runway_5p_matcher; }

const SeqMatcher& LutBundleYs::Sid5pMatcher() const { return _sid_5p_matcher; }

const BarcodePool& LutBundleYs::Sid5pPool() const { return _sid_5p_matcher.Pool(); }

const SeqMatcher& LutBundleYs::SidSpacer5pMatcher() const { return _sid_spacer_5p_matcher; }

const BarcodePool& LutBundleYs::SidSpacer5pPool() const { return _sid_spacer_5p_matcher.Pool(); }

const SeqMatcher& LutBundleYs::Sid3pMatcher() const { return _sid_3p_matcher; }

const BarcodePool& LutBundleYs::Sid3pPool() const { return _sid_3p_matcher.Pool(); }

const SeqMatcher& LutBundleYs::SidSpacer3pMatcher() const { return _sid_spacer_3p_matcher; }

const BarcodePool& LutBundleYs::SidSpacer3pPool() const { return _sid_spacer_3p_matcher.Pool(); }

template <>
LutBundleYs CreateLutBundle(const AdapterDesignBundle& designs) {
  using enum BarcodeType;
  return LutBundleYs{*designs.GetMatcher5p(kRunway), *designs.GetMatcher5p(kSid), *designs.GetMatcher5p(kStem),
                     *designs.GetMatcher3p(kSid), *designs.GetMatcher3p(kStem)};
}

}  // namespace xoos::demux
