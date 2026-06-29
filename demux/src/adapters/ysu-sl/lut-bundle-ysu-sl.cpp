#include "lut-bundle-ysu-sl.h"

#include <utility>

#include "adapter-design/adapter-design-bundle.h"
// Clangd's include-cleaner doesn't detect usage from template specializations.
// This include is required for explicit specialization below.
#include "lut-bundle/lut-bundle.h"  // IWYU pragma: keep

namespace xoos::demux {

LutBundleYsuSl::LutBundleYsuSl(SeqMatcher runway_5p, SeqMatcher sid_5p, SeqMatcher sid_umi_spacer_5p, SeqMatcher umi_5p,
                               SeqMatcher sid_3p, SeqMatcher umi_sid_spacer_3p, SeqMatcher umi_3p)
    : _runway_5p_matcher{std::move(runway_5p)},
      _sid_5p_matcher{std::move(sid_5p)},
      _sid_umi_spacer_5p_matcher{std::move(sid_umi_spacer_5p)},
      _umi_5p_matcher{std::move(umi_5p)},
      _sid_3p_matcher{std::move(sid_3p)},
      _umi_sid_spacer_3p_matcher{std::move(umi_sid_spacer_3p)},
      _umi_3p_matcher{std::move(umi_3p)} {}

const SeqMatcher& LutBundleYsuSl::Runway5pMatcher() const { return _runway_5p_matcher; }

const SeqMatcher& LutBundleYsuSl::Sid5pMatcher() const { return _sid_5p_matcher; }

const BarcodePool& LutBundleYsuSl::Sid5pPool() const { return _sid_5p_matcher.Pool(); }

const SeqMatcher& LutBundleYsuSl::SidUmiSpacer5pMatcher() const { return _sid_umi_spacer_5p_matcher; }

const BarcodePool& LutBundleYsuSl::SidUmiSpacer5pPool() const { return _sid_umi_spacer_5p_matcher.Pool(); }

const SeqMatcher& LutBundleYsuSl::Umi5pMatcher() const { return _umi_5p_matcher; }

const BarcodePool& LutBundleYsuSl::Umi5pPool() const { return _umi_5p_matcher.Pool(); }

const SeqMatcher& LutBundleYsuSl::Sid3pMatcher() const { return _sid_3p_matcher; }

const BarcodePool& LutBundleYsuSl::Sid3pPool() const { return _sid_3p_matcher.Pool(); }

const SeqMatcher& LutBundleYsuSl::SidUmiSpacer3pMatcher() const { return _umi_sid_spacer_3p_matcher; }

const BarcodePool& LutBundleYsuSl::SidUmiSpacer3pPool() const { return _umi_sid_spacer_3p_matcher.Pool(); }

const SeqMatcher& LutBundleYsuSl::Umi3pMatcher() const { return _umi_3p_matcher; }

const BarcodePool& LutBundleYsuSl::Umi3pPool() const { return _umi_3p_matcher.Pool(); }

template <>
LutBundleYsuSl CreateLutBundle(const AdapterDesignBundle& designs) {
  return LutBundleYsuSl{
      *designs.GetMatcher5p(BarcodeType::kRunway), *designs.GetMatcher5p(BarcodeType::kSid),
      *designs.GetMatcher5p(BarcodeType::kStem),   *designs.GetMatcher5p(BarcodeType::kUmi),
      *designs.GetMatcher3p(BarcodeType::kSid),    *designs.GetMatcher3p(BarcodeType::kStem),
      *designs.GetMatcher3p(BarcodeType::kUmi),
  };
}

}  // namespace xoos::demux
