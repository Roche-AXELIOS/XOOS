#include "lut-bundle-duplex.h"

#include <fmt/format.h>

#include <utility>

#include "adapter-design/adapter-design-bundle.h"
// Clangd's include-cleaner doesn't detect usage from template specializations.
// This include is required for explicit specialization below.
#include <xoos/error/error.h>

#include "lut-bundle/lut-bundle.h"  // IWYU pragma: keep

namespace xoos::demux {
// TODO: rewrite to use a struct
LutBundleDuplex::LutBundleDuplex(SeqMatcher sid_5p, SeqMatcher start_sequence, SeqMatcher sid_3p,
                                 SeqMatcher stop_sequence, SeqMatcher loop_sequence)
    : _sid_5p_matcher{std::move(sid_5p)},
      _start_sequence_matcher{std::move(start_sequence)},
      _sid_3p_matcher{std::move(sid_3p)},
      _stop_sequence_matcher{std::move(stop_sequence)},
      _loop_sequence_matcher{std::move(loop_sequence)} {}

const SeqMatcher& LutBundleDuplex::Sid5pMatcher() const { return _sid_5p_matcher; }

const BarcodePool& LutBundleDuplex::Sid5pPool() const { return _sid_5p_matcher.Pool(); }

const SeqMatcher& LutBundleDuplex::StartSequenceMatcher() const { return _start_sequence_matcher; }

const SeqMatcher& LutBundleDuplex::Sid3pMatcher() const { return _sid_3p_matcher; }

const SeqMatcher& LutBundleDuplex::StopSequenceMatcher() const { return _stop_sequence_matcher; }

std::string_view LutBundleDuplex::LoopSequence() const {
  const auto& loop_pool = _loop_sequence_matcher.Pool();
  if (loop_pool.empty()) {
    throw error::Error("Loop pool is empty; expected exactly one loop sequence for duplex demux/trim.");
  }
  if (loop_pool.size() != 1U) {
    throw error::Error(
        fmt::format("Loop pool contains {} sequences; expected exactly one loop sequence for duplex demux/trim.",
                    loop_pool.size()));
  }
  const auto& loop_sequence = loop_pool.front().sequence;
  return loop_sequence;
}

template <>
LutBundleDuplex CreateLutBundle<LutBundleDuplex>(const AdapterDesignBundle& designs) {
  return LutBundleDuplex{*designs.GetMatcher5p(BarcodeType::kSid), *designs.GetMatcher5p(BarcodeType::kAnchor),
                         *designs.GetMatcher3p(BarcodeType::kSid), *designs.GetMatcher3p(BarcodeType::kAnchor),
                         *designs.GetMatcher5p(BarcodeType::kLoop)};
}

}  // namespace xoos::demux
