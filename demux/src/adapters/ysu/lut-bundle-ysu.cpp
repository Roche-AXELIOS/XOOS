#include "lut-bundle-ysu.h"

#include <xoos/error/error.h>

#include "adapter-design/adapter-component-indices.h"
#include "adapter-design/adapter-design-bundle.h"
// Clangd's include-cleaner doesn't detect usage from template specializations.
// This include is required for explicit specialization below.
#include "lut-bundle/lut-bundle.h"  // IWYU pragma: keep

namespace xoos::demux {

// Positional constructor; adjacent SeqMatcher params encode fixed read-structure positions
// documented in the class header and are constructed positionally by the assembler.
LutBundleYsu::LutBundleYsu(SeqMatcher sid_1,  // NOSONAR(S5419)
                           SeqMatcher umi_2,
                           SeqMatcher umi_3,
                           SeqMatcher sid_4,
                           std::optional<std::string> bait_1,
                           std::optional<std::string> bait_4)
    : _sid_1_matcher{std::move(sid_1)},
      _umi_2_matcher{std::move(umi_2)},
      _umi_3_matcher{std::move(umi_3)},
      _sid_4_matcher{std::move(sid_4)},
      _bait_1{std::move(bait_1)},
      _bait_4{std::move(bait_4)} {
}

const SeqMatcher& LutBundleYsu::Sid1Matcher() const {
  return _sid_1_matcher;
}

const BarcodePool& LutBundleYsu::Sid1Pool() const {
  return _sid_1_matcher.Pool();
}

const SeqMatcher& LutBundleYsu::Umi2Matcher() const {
  return _umi_2_matcher;
}

const SeqMatcher& LutBundleYsu::Umi3Matcher() const {
  return _umi_3_matcher;
}

const SeqMatcher& LutBundleYsu::Sid4Matcher() const {
  return _sid_4_matcher;
}

const std::optional<std::string>& LutBundleYsu::Bait1() const {
  return _bait_1;
}

const std::optional<std::string>& LutBundleYsu::Bait4() const {
  return _bait_4;
}

namespace {
/// @brief A component's bait as an optional: empty string -> nullopt so the trimmer skips the search.
std::optional<std::string> BaitOf(const AdapterComponent& component) {
  if (component.bait.empty()) {
    return std::nullopt;
  }
  return component.bait;
}
}  // namespace

template <>
LutBundleYsu CreateLutBundle(const AdapterDesignBundle& designs) {
  // Read-order component layout: [sid_1, umi_2, umi_3, sid_4]. Each sid carries its own bait seed.
  using Index = adapter_component_index::Ysu;
  return LutBundleYsu{
      *designs.MatcherAt(Index::kSidFwd),
      *designs.MatcherAt(Index::kUmiFwd),
      *designs.MatcherAt(Index::kUmiRev),
      *designs.MatcherAt(Index::kSidRev),
      BaitOf(designs.components.at(Index::kSidFwd)),
      BaitOf(designs.components.at(Index::kSidRev)),
  };
}

}  // namespace xoos::demux
