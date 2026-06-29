#include "match-info.h"

#include <fmt/format.h>

#include <utility>

namespace xoos::demux {
Barcode::Barcode() : Barcode(0, "", "") {}

Barcode::Barcode(uint id, std::string sequence, std::string name)
    : id{id}, sequence{std::move(sequence)}, name{std::move(name)} {}

std::optional<uint> MatchInfo::BarcodeId(const std::optional<MatchInfo>& match_info) {
  return !match_info ? std::nullopt : std::make_optional(match_info->BarcodeId());
}

std::optional<uint> MatchInfo::EDist(const std::optional<MatchInfo>& match_info) {
  return !match_info ? std::nullopt : std::make_optional(match_info->EDist());
}

std::optional<LociRange> MatchInfo::Position(const std::optional<MatchInfo>& match_info) {
  return !match_info ? std::nullopt : std::make_optional(match_info->Position());
}

MatchInfo::MatchInfo(const BarcodeMatch& match, const LociRange& position, MatchType match_type)
    : _match{match}, _position{position}, _match_type{match_type} {}

void MatchInfo::Update(const BarcodeMatch& new_match, MatchType new_match_type, uint spos, uint epos) {
  // This function is called if we have found a match - but we need to determine how good that match is. We
  // keep track of the best match.
  if (new_match.edist > EDist()) {
    return;  // match is worse, ignore.
  }
  if (new_match.edist < EDist()) {
    // This is the best match we have seen so far, update.
    _match = new_match;
    _position = LociRange{spos, epos};
    _match_type = new_match_type;
    return;
  }

  // If we get here, we have a match with the same edit distance as the current best match. This usually
  // doesn't bode well as these are often ambiguous matches. There is one exception: if we find a match with
  // the same barcode_id with a longer length, we prefer that match.
  if (new_match.barcode_id == BarcodeId()) {  // so we might find multiple matches at similar positions (as in:
                                              // variations in position with respect to
    // the center is maximum edit distance.
    const auto match_len{Length()};
    const auto new_match_len{epos - spos};

    if (new_match_len > match_len) {  // this is good news, we have a longer match with the same barcode_id.
      _position = LociRange{spos, epos};
      _match_type = new_match_type;
    } else if (new_match_len == match_len) {  // this is bad news, we have a match with the same length and barcode_id.
      // Before declaring an ambiguous match, please note that with an edit distance of 2, we can have two matches
      // that are not ambiguous. I discovered this by printing out the differences, e.g.
      // "Ambiguous match found: id 597 at 4 to 18 conflicts with position 5 to 19 (edit distance 2)"
      if (std::abs(static_cast<int>(spos) - static_cast<int>(_position.spos)) < EDist()) {
        return;  // not ambiguous! The match we already found is OK.
      }

      _match_type = MatchType::kAmbiguous;
    }
    // The case of new_match_len < match_len is ignored, as we already have a longer match.
    return;
  }

  // If we get here, we have a match with the same edit distance as the current best match, but with a different
  // barcode_id. This is bad news, as we have an ambiguous match. This case was NOT CAUGHT in demuxing code
  // until April 2024 (discovered by KZ) - basically, we were not updating the match_type to ambiguous...
  _match_type = MatchType::kAmbiguous;  // This statement not present before April 2024.
}
}  // namespace xoos::demux
