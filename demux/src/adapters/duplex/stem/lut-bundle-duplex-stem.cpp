#include "lut-bundle-duplex-stem.h"

#include "adapter-design/adapter-component-indices.h"
#include "adapter-design/adapter-design-bundle.h"
#include "adapters/duplex/lut-bundle-duplex.h"
#include "lut-bundle/lut-bundle.h"

namespace xoos::demux {

// Positional constructor; adjacent FixedSequence params encode fixed read-structure positions
// documented in the class header and are constructed positionally by the assembler.
LutBundleDuplexStem::LutBundleDuplexStem(FixedSequence runway_1,  // NOSONAR(S5419)
                                         SeqMatcher sid_2,
                                         FixedSequence start_3,
                                         FixedSequence stem_4,
                                         FixedSequence loop_5,
                                         SeqMatcher sid_6,
                                         [[maybe_unused]] const FixedSequence& stem_8)
    : LutBundleDuplex(std::move(runway_1), std::move(sid_2), std::move(start_3), std::move(loop_5), std::move(sid_6)),
      _stem_4{std::move(stem_4)} {
}

const FixedSequence& LutBundleDuplexStem::Stem4() const {
  return _stem_4;
}

template <>
LutBundleDuplexStem CreateLutBundle<LutBundleDuplexStem>(const AdapterDesignBundle& designs) {
  // Read-order component layout: [runway, sid, start, stem, loop, sid, anchor, stem]. Only the
  // pre-loop stem (index 3) is consumed by this adapter; the stop anchor and post-loop stem are not read.
  using Index = adapter_component_index::DuplexStem;
  return LutBundleDuplexStem{designs.FixedAt(Index::kRunway),
                             *designs.MatcherAt(Index::kSidFwd),
                             designs.FixedAt(Index::kStart),
                             designs.FixedAt(Index::kStemPre),
                             designs.FixedAt(Index::kLoop),
                             *designs.MatcherAt(Index::kSidRev),
                             designs.FixedAt(Index::kStemPost)};
}

}  // namespace xoos::demux
