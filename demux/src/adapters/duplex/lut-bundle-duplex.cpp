#include "lut-bundle-duplex.h"

#include "adapter-design/adapter-component-indices.h"
#include "adapter-design/adapter-design-bundle.h"
// Clangd's include-cleaner doesn't detect usage from template specializations.
// This include is required for explicit specialization below.
#include "lut-bundle/lut-bundle.h"  // IWYU pragma: keep

namespace xoos::demux {

LutBundleDuplex::LutBundleDuplex(
    FixedSequence runway_1, SeqMatcher sid_2, FixedSequence start_3, FixedSequence loop_4, SeqMatcher sid_6)
    : _runway_1{std::move(runway_1)},
      _sid_2_matcher{std::move(sid_2)},
      _start_3{std::move(start_3)},
      _loop_4{std::move(loop_4)},
      _sid_6_matcher{std::move(sid_6)} {
}

const FixedSequence& LutBundleDuplex::Runway1() const {
  return _runway_1;
}

const SeqMatcher& LutBundleDuplex::Sid2Matcher() const {
  return _sid_2_matcher;
}

const BarcodePool& LutBundleDuplex::Sid2Pool() const {
  return _sid_2_matcher.Pool();
}

const FixedSequence& LutBundleDuplex::StartSequence() const {
  return _start_3;
}

std::string_view LutBundleDuplex::LoopSequence() const {
  return _loop_4.Seq();
}

const SeqMatcher& LutBundleDuplex::Sid6Matcher() const {
  return _sid_6_matcher;
}

template <>
LutBundleDuplex CreateLutBundle<LutBundleDuplex>(const AdapterDesignBundle& designs) {
  // Read-order component layout: [runway, sid, start, loop, anchor, sid]. Components are addressed by
  // position; the two sid matchers differ only in orientation, and the dead stop anchor (index 4) is
  // not read.
  using Index = adapter_component_index::Duplex;
  return LutBundleDuplex{designs.FixedAt(Index::kRunway),
                         *designs.MatcherAt(Index::kSidFwd),
                         designs.FixedAt(Index::kStart),
                         designs.FixedAt(Index::kLoop),
                         *designs.MatcherAt(Index::kSidRev)};
}

}  // namespace xoos::demux
