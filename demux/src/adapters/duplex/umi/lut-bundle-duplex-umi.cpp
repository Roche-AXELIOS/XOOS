#include "lut-bundle-duplex-umi.h"

#include <utility>

#include "adapter-design/adapter-design-bundle.h"
#include "adapters/duplex/stem/lut-bundle-duplex-stem.h"
#include "lut-bundle/lut-bundle.h"

namespace xoos::demux {

// TODO: rewrite to use a struct
LutBundleDuplexUmi::LutBundleDuplexUmi(SeqMatcher sid_5p, SeqMatcher start_sequence, SeqMatcher sid_3p,
                                       SeqMatcher stop_sequence, SeqMatcher stem_5p, SeqMatcher stem_3p,
                                       SeqMatcher umi_5p, SeqMatcher umi_3p, SeqMatcher loop_sequence_matcher)
    : LutBundleDuplexStem(std::move(sid_5p), std::move(start_sequence), std::move(sid_3p), std::move(stop_sequence),
                          std::move(stem_5p), std::move(stem_3p), std::move(loop_sequence_matcher)),
      _umi_5p_matcher{std::move(umi_5p)},
      _umi_3p_matcher{std::move(umi_3p)} {}

const SeqMatcher& LutBundleDuplexUmi::Umi5pMatcher() const { return _umi_5p_matcher; }

const BarcodePool& LutBundleDuplexUmi::Umi5pPool() const { return _umi_5p_matcher.Pool(); }

const SeqMatcher& LutBundleDuplexUmi::Umi3pMatcher() const { return _umi_3p_matcher; }

const BarcodePool& LutBundleDuplexUmi::Umi3pPool() const { return _umi_3p_matcher.Pool(); }

template <>
LutBundleDuplexUmi CreateLutBundle<LutBundleDuplexUmi>(const AdapterDesignBundle& designs) {
  return LutBundleDuplexUmi{*designs.GetMatcher5p(BarcodeType::kSid),  *designs.GetMatcher5p(BarcodeType::kAnchor),
                            *designs.GetMatcher3p(BarcodeType::kSid),  *designs.GetMatcher3p(BarcodeType::kAnchor),
                            *designs.GetMatcher5p(BarcodeType::kStem), *designs.GetMatcher3p(BarcodeType::kStem),
                            *designs.GetMatcher5p(BarcodeType::kUmi),  *designs.GetMatcher3p(BarcodeType::kUmi),
                            *designs.GetMatcher5p(BarcodeType::kLoop)};
}

}  // namespace xoos::demux
