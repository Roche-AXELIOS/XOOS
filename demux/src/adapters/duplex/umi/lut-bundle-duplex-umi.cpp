#include "lut-bundle-duplex-umi.h"

#include "adapter-design/adapter-component-indices.h"
#include "adapter-design/adapter-design-bundle.h"
#include "adapters/duplex/stem/lut-bundle-duplex-stem.h"
#include "lut-bundle/lut-bundle.h"

namespace xoos::demux {

// Positional constructor; adjacent FixedSequence params encode fixed read-structure positions
// documented in the class header and are constructed positionally by the assembler.
LutBundleDuplexUmi::LutBundleDuplexUmi(FixedSequence runway_1,  // NOSONAR(S5419) NOSONAR(S107)
                                       SeqMatcher sid_2,
                                       FixedSequence start_3,
                                       FixedSequence stem_4,
                                       FixedSequence loop_5,
                                       SeqMatcher umi_6,
                                       SeqMatcher sid_7,
                                       const FixedSequence& stem_9,
                                       SeqMatcher umi_10)
    : LutBundleDuplexStem(std::move(runway_1),
                          std::move(sid_2),
                          std::move(start_3),
                          std::move(stem_4),
                          std::move(loop_5),
                          std::move(sid_7),
                          stem_9),
      _umi_6_matcher{std::move(umi_6)},
      _umi_10_matcher{std::move(umi_10)} {
}

const SeqMatcher& LutBundleDuplexUmi::Umi6Matcher() const {
  return _umi_6_matcher;
}

const SeqMatcher& LutBundleDuplexUmi::Umi10Matcher() const {
  return _umi_10_matcher;
}

template <>
LutBundleDuplexUmi CreateLutBundle<LutBundleDuplexUmi>(const AdapterDesignBundle& designs) {
  // Read-order component layout: [runway, sid, start, stem, loop, umi, sid, anchor, stem, umi].
  // Components are addressed by position; the two sid and the two umi matchers each differ only in
  // orientation, and the dead stop anchor (index 7) is not read.
  using Index = adapter_component_index::DuplexUmi;
  return LutBundleDuplexUmi{designs.FixedAt(Index::kRunway),
                            *designs.MatcherAt(Index::kSidFwd),
                            designs.FixedAt(Index::kStart),
                            designs.FixedAt(Index::kStemPre),
                            designs.FixedAt(Index::kLoop),
                            *designs.MatcherAt(Index::kUmiFwd),
                            *designs.MatcherAt(Index::kSidRev),
                            designs.FixedAt(Index::kStemPost),
                            *designs.MatcherAt(Index::kUmiRev)};
}

}  // namespace xoos::demux
