#pragma once

#include <cstddef>
#include <optional>
#include <span>
#include <string_view>
#include <unordered_map>

#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::demux {

/// @brief Controls which single-ended matches may assign reads for dual-SID sample rows.
enum class DualSidMode {
  /// Partial matching is determined only by 5p/3p SID collisions in the sample sheet.
  kNone,
  /// 3p-only reads are not assigned; 5p-only reads may still assign if their 5p SID is unique.
  kNeed5p,
  /// 5p-only reads are not assigned; 3p-only reads may still assign if their 3p SID is unique.
  kNeed3p,
  /// Single-ended reads are not assigned.
  kNeedBoth,
};

/// @brief Per-sample dual-SID partial assignment policy.
struct DualSidPartialMatch {
  bool allow_5p_partial_match{false};
  bool allow_3p_partial_match{false};
};

/// @brief One sample-sheet pairing reachable from a 3' SID matcher ID.
struct DualSidPairing {
  /// @brief Matching 5' SID matcher ID.
  u32 sid_5p_id;
  /// @brief Sample-sheet row represented by this SID pair.
  u32 sample_id;
};

inline u64 DualSidPairKey(const u32 sid_5p_id, const u32 sid_3p_id) {
  // Pack the 5p side ID into the high 32 bits and the 3p side ID into the low 32 bits.
  constexpr u32 kSid5pShiftBits = 32U;
  return (static_cast<u64>(sid_5p_id) << kSid5pShiftBits) | static_cast<u64>(sid_3p_id);
}

/// @brief Maps unique per-side SID matcher IDs back to sample row IDs.
struct DualSidMatchingPolicy {
  /// @brief Per-sample partial assignment policy indexed by sample row ID.
  vec<DualSidPartialMatch> partial_matches;
  /// @brief Unique 5p side ID -> sample row ID when that side is unambiguous.
  vec<std::optional<u32>> sample_by_5p_sid;
  /// @brief Unique 3p side ID -> sample row ID when that side is unambiguous.
  vec<std::optional<u32>> sample_by_3p_sid;
  /// @brief Packed (5p side ID, 3p side ID) -> sample row ID.
  std::unordered_map<u64, u32> sample_by_sid_pair;
  /// @brief CSR offsets into @ref pairings_by_3p_sid, indexed by 3' SID matcher ID.
  vec<u32> pairing_offsets_by_3p_sid;
  /// @brief Sample-sheet pairings grouped by 3' SID matcher ID.
  vec<DualSidPairing> pairings_by_3p_sid;

  bool Enabled() const {
    return !sample_by_sid_pair.empty();
  }

  /// @brief Returns every sample-sheet pair compatible with a 3' SID matcher ID.
  std::span<const DualSidPairing> PairingsFor3pSid(const u32 sid_3p_id) const {
    const auto offset_index = static_cast<std::size_t>(sid_3p_id);
    if (offset_index + 1 >= pairing_offsets_by_3p_sid.size()) {
      return {};
    }
    const auto begin = pairing_offsets_by_3p_sid.at(offset_index);
    const auto end = pairing_offsets_by_3p_sid.at(offset_index + 1);
    return std::span<const DualSidPairing>{pairings_by_3p_sid}.subspan(begin, end - begin);
  }
};

constexpr std::string_view ToString(const DualSidMode mode) {
  using enum DualSidMode;
  switch (mode) {
    case kNone:
      return "none";
    case kNeed5p:
      return "need-5p";
    case kNeed3p:
      return "need-3p";
    case kNeedBoth:
      return "need-both";
    default:
      return "none";
  }
}

inline std::optional<DualSidMode> ParseDualSidMode(const std::string_view value) {
  using enum DualSidMode;
  if (value == ToString(kNone)) {
    return kNone;
  }
  if (value == ToString(kNeed5p)) {
    return kNeed5p;
  }
  if (value == ToString(kNeed3p)) {
    return kNeed3p;
  }
  if (value == ToString(kNeedBoth)) {
    return kNeedBoth;
  }
  return std::nullopt;
}

}  // namespace xoos::demux
