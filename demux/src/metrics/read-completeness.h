#pragma once

#include <optional>

#include <xoos/error/error.h>
#include <xoos/types/int.h>

namespace xoos::demux {

enum class DiscordantSidMode {
  // Discard reads with discordant SIDs that are equally likely (default)
  kDiscardTied,
  // Discard all reads with discordant SIDs (both significantly detected)
  kDiscardAll,
  // Keep all discordant reads, assign to winner
  kKeep,
};

/**
 * @brief Determines whether a read is full: a landmark (SID or UMI) detected on both ends.
 *
 * For UMI adapters: (SID or UMI) on each end. For non-UMI adapters: a SID on each end.
 *
 * @tparam T Any trim-info type with sid_5p/sid_3p fields and optionally umi_5p/umi_3p.
 * @param[in] record The trim-info record to evaluate.
 * @return True if the read qualifies as a full read.
 */
template <class T>
bool IsFullRead(const T& record) {
  if constexpr (requires { record.umi_5p && record.umi_3p; }) {
    return (record.sid_5p || record.umi_5p) && (record.sid_3p || record.umi_3p);
  } else {
    return (record.sid_5p && record.sid_3p);
  }
}

/**
 * @brief Determines whether a read is partial: the complement of IsFullRead.
 *
 * @tparam T Any trim-info type with sid_5p/sid_3p fields and optionally umi_5p/umi_3p.
 * @param[in] record The trim-info record to evaluate.
 * @return True if the read qualifies as a partial read.
 */
template <class T>
bool IsPartialRead(const T& record) {
  return !IsFullRead(record);
}

/**
 * @brief Whether a YSU read has both SID and UMI detected on at least one end (the UMI gate).
 *
 * @tparam T A UMI trim-info type with sid_5p/sid_3p and umi_5p/umi_3p fields.
 * @param[in] record The trim-info record to evaluate.
 * @return True if at least one end has both SID and UMI detected.
 */
template <class T>
bool HasSidUmiCombination(const T& record) {
  static_assert(
      requires { record.umi_5p && record.umi_3p; },
      "HasSidUmiCombination is only defined for UMI adapter trim-info types");
  return (record.sid_5p && record.umi_5p) || (record.sid_3p && record.umi_3p);
}

/**
 * @brief Determines whether a read is a perfect index match (zero edit distance on both SIDs, concordant).
 *
 * @tparam T Any trim-info type with sid_5p/sid_3p and sid_5p_edist/sid_3p_edist fields.
 * @param[in] record The trim-info record to evaluate.
 * @return True if both SIDs match with zero edits and agree.
 */
template <class T>
bool IsPerfectMatch(const T& record) {
  const bool both_sids_concordant = record.sid_5p && record.sid_3p && record.sid_5p == record.sid_3p;
  const bool zero_edit_distance = record.sid_5p_edist == 0 && record.sid_3p_edist == 0;
  return both_sids_concordant && zero_edit_distance;
}

/**
 * @brief Determines whether both the 5' and 3' SIDs were detected.
 *
 * @tparam T Any trim-info type with sid_5p/sid_3p fields.
 * @param[in] record The trim-info record to evaluate.
 * @return True if both SIDs are present.
 */
template <class T>
bool IsBothSidDetected(const T& record) {
  return record.sid_5p && record.sid_3p;
}

/**
 * @brief Determines whether both SIDs were detected but assigned to different indices (discordant).
 *
 * @tparam T Any trim-info type with sid_5p/sid_3p fields.
 * @param[in] record The trim-info record to evaluate.
 * @return True if both SIDs are present but disagree.
 */
template <class T>
bool IsIndexDiscordant(const T& record) {
  return record.sid_5p && record.sid_3p && record.sid_5p != record.sid_3p;
}

/**
 * @brief Determines whether a discordant read should be discarded based on the configured mode.
 *
 * For kKeep, always returns false. For kDiscardAll, returns true for any discordant read.
 * For kDiscardTied, checks whether the two sides are equally likely.
 *
 * @tparam TrimInfo Any trim-info type with sid_5p/sid_3p and sid_5p_edist/sid_3p_edist fields.
 * @param[in] record The trim-info record to evaluate.
 * @param[in] mode   The configured discordant SID handling mode.
 * @return True if the read should be discarded.
 */
template <typename TrimInfo>
bool ShouldDiscardDiscordant(const TrimInfo& record, const DiscordantSidMode mode) {
  const std::optional<u32>& sid_5p = record.sid_5p;
  const std::optional<u32>& sid_3p = record.sid_3p;

  // Both SIDs must be present to evaluate discordance
  if (!sid_5p.has_value() && !sid_3p.has_value()) [[unlikely]] {
    throw error::Error("SIDs found unset before checking discordance");
  }
  if (!IsIndexDiscordant(record)) {
    return false;
  }
  switch (mode) {
    using enum DiscordantSidMode;
    case kDiscardAll:
      return true;
    case kDiscardTied:
      // Discard only when equally likely. Score-bearing trim-info types (e.g. simplex) tie-break on
      // score; others fall back to edit distance.
      if constexpr (requires { record.score_5p == record.score_3p; }) {
        return (record.score_5p == record.score_3p);
      } else {
        return (record.sid_5p_edist == record.sid_3p_edist);
      }
    case kKeep:
    default:
      return false;
  }
}

}  // namespace xoos::demux
