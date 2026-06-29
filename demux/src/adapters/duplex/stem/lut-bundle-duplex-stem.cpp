#include "lut-bundle-duplex-stem.h"

#include <utility>

#include "adapter-design/adapter-design-bundle.h"
#include "adapters/duplex/lut-bundle-duplex.h"
#include "lut-bundle/lut-bundle.h"

namespace xoos::demux {

// SeqMatcher sid_5p, SeqMatcher start_sequence, SeqMatcher sid_3p,
//                                  SeqMatcher stop_sequence, SeqMatcher loop_sequence
LutBundleDuplexStem::LutBundleDuplexStem(SeqMatcher sid_5p, SeqMatcher start_sequence, SeqMatcher sid_3p,
                                         SeqMatcher stop_sequence, SeqMatcher stem_5p, SeqMatcher stem_3p,
                                         SeqMatcher loop_sequence_matcher)
    : LutBundleDuplex(std::move(sid_5p), std::move(start_sequence), std::move(sid_3p), std::move(stop_sequence),
                      std::move(loop_sequence_matcher)),
      _stem_5p_matcher{std::move(stem_5p)},
      _stem_3p_matcher{std::move(stem_3p)} {}

const SeqMatcher& LutBundleDuplexStem::Stem5pMatcher() const { return _stem_5p_matcher; }

const SeqMatcher& LutBundleDuplexStem::Stem3pMatcher() const { return _stem_3p_matcher; }

template <>
LutBundleDuplexStem CreateLutBundle<LutBundleDuplexStem>(const AdapterDesignBundle& designs) {
  using enum BarcodeType;
  return LutBundleDuplexStem{*designs.GetMatcher5p(kSid),  *designs.GetMatcher5p(kAnchor),
                             *designs.GetMatcher3p(kSid),  *designs.GetMatcher3p(kAnchor),
                             *designs.GetMatcher5p(kStem), *designs.GetMatcher3p(kStem),
                             *designs.GetMatcher5p(kLoop)};
}

}  // namespace xoos::demux
