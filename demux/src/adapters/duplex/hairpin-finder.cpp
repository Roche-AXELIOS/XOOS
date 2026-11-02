#include "hairpin-finder.h"

#include <map>  // IWYU pragma: keep

#include <xoos/enum/enum-util.h>
#include <xoos/types/int.h>

#include "duplex-match.h"
#include "sequence/encoding-constants.h"

namespace xoos::demux {

constexpr u8 kBitapMaxEditDistance = 4;
/// @brief Minimum loop search window length, preserving the legacy 32-base floor for short loops.
constexpr s32 kMinLoopSearchWindow = 32;
/// @brief The loop search window is sized as this multiple of the loop length.
constexpr s32 kLoopSearchWindowMultiplier = 3;

HairpinFinder::HairpinFinder(SeqMatcher sid_5p_matcher, SeqMatcher sid_3p_matcher, const std::string_view loop_sequence)
    : _sid_5p_matcher(std::move(sid_5p_matcher)),
      _sid_3p_matcher{std::move(sid_3p_matcher)},
      _loop_sequence{loop_sequence},
      _cascaded_luts(
          _sid_5p_matcher, _sid_3p_matcher, DuplexMatch::BarcodeType::kSID5p, DuplexMatch::BarcodeType::kSID3p),
      _mask_sid5p(
          static_cast<int64_t>((1UL << (kBitsPerEncodedBase * _cascaded_luts.MaxLength(DuplexMatch::kSID5p))) - 1UL)),
      _mask_sid3p(
          static_cast<int64_t>((1UL << (kBitsPerEncodedBase * _cascaded_luts.MaxLength(DuplexMatch::kSID3p))) - 1UL)),
      _mask_plan_a_full{_mask_sid5p, _mask_sid3p, _mask_sid5p, _mask_sid3p, 0, 0, 0, 0},
      _loop_fw(_loop_sequence, SearchDirection::kForward),
      _loop_bw(_loop_sequence, SearchDirection::kBackward),
      _loop_length(static_cast<s32>(_loop_sequence.length())),
      _sid_5p_length(static_cast<s32>(_sid_5p_matcher.Pool().front().sequence.length())),
      _sid_3p_length(static_cast<s32>(_sid_3p_matcher.Pool().front().sequence.length())),
      _loop_search_window(std::min(std::max(kLoopSearchWindowMultiplier * _loop_length, kMinLoopSearchWindow),
                                   Bitap<2>::kQueryWindowSize)) {
  // Initialize the bitap's using a std::map, then convert to a std::vector. Note that for performance reasons, most
  // member variables in the bitap are const, so we need to initialize them upon construction which makes the use of
  // a vector more difficult.
  std::map<u32, Bitap<kBitapMaxEditDistance>> bitap_3p;
  std::map<u32, Bitap<kBitapMaxEditDistance>> bitap_3p_bw;
  for (const auto& sid : _sid_3p_matcher.Pool()) {
    bitap_3p.try_emplace(sid.id, sid.sequence, SearchDirection::kForward);
    // Backward counterpart, paired with `bitap_3p` for `FindStartEnd` calls.
    bitap_3p_bw.try_emplace(sid.id, sid.sequence, SearchDirection::kBackward);
  }
  // The map will sort the ids, so figure out the max id and create a vector.
  // We should expect that SID list is not empty.
  const auto max_id = bitap_3p.rbegin()->first;

  _sid_3p_bitap.reserve(max_id + 1);
  _sid_3p_bw.reserve(max_id + 1);

  // SID IDs are contiguous starting from 0 — use at() to assert this invariant.
  for (u32 i = 0; i <= max_id; ++i) {
    _sid_3p_bitap.push_back(bitap_3p.at(i));
    _sid_3p_bw.push_back(bitap_3p_bw.at(i));
  }
}

HairpinFinder::HairpinFinder(const LutBundleDuplex& lut_bundle)
    : HairpinFinder(lut_bundle.Sid2Matcher(), lut_bundle.Sid6Matcher(), lut_bundle.LoopSequence()) {
}

}  // namespace xoos::demux
