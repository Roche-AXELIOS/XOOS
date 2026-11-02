#include "lut-bundle-simplex-10x.h"

#include "adapter-design/adapter-component-indices.h"
#include "adapter-design/adapter-design-bundle.h"
#include "lut-bundle/lut-bundle.h"

namespace xoos::demux {

LutBundleSimplex10x::LutBundleSimplex10x(FixedSequence fixed_1, SeqMatcher sid_2, FixedSequence fixed_3)
    : fixed_1{std::move(fixed_1)}, sid_2_matcher{std::move(sid_2)}, fixed_3{std::move(fixed_3)} {
}

template <>
LutBundleSimplex10x CreateLutBundle<LutBundleSimplex10x>(
    const AdapterDesignBundle& designs) {  // NOSONAR(S5536) -- called from flow-manager
  // Read-order component layout: [anchor, sid, stem].
  using Index = adapter_component_index::Simplex10x;
  return LutBundleSimplex10x{
      designs.FixedAt(Index::kAnchor), *designs.MatcherAt(Index::kSid), designs.FixedAt(Index::kStem)};
}

}  // namespace xoos::demux
