#pragma once

#include <optional>
#include <string>
#include <vector>

#include <xoos/types/int.h>

#include "sequence/alignment/log-likelihood-scoring.h"
#include "sequence/loci-range.h"

namespace xoos::demux {
enum class MatchType {
  kExact,
  kUnknown,
  kAmbiguous,
  kInExact,
};

// To minimize memory usage, we'd like to limit the amount of memory used by this structure to 16 bits.
// This would then allow us to put all the barcode matches in a very large 8 GB (2^32) lookup table.
static constexpr u16 kUninitialized{3};

struct BarcodeMatch {
  BarcodeMatch() : barcode_id(0), edist(kUninitialized) {
  }

  BarcodeMatch(const u16 id, const u16 dist) : barcode_id(id), edist(dist) {
  }

  u16 barcode_id : 11;  // NOSONAR cpp:S2806 - u16 packs BarcodeMatch into 16 bits; max 2048 IDs
  u16 edist : 2;        // NOSONAR cpp:S2806 - u16 packs BarcodeMatch into 16 bits; 0..2, 3 = uninitialized
};

static_assert(sizeof(BarcodeMatch) == 2, "BarcodeMatch must be exactly 16 bits");

/**
 * @brief A named sequence component within an adapter design.
 *
 * In this codebase "barcode" is used broadly to mean any recognizable short
 * sub-sequence that the LUT matches against — SIDs, UMIs, spacers, runways,
 * stems, etc. (see BarcodeType). The struct simply pairs an ordinal ID with
 * the raw nucleotide sequence and a human-readable name.
 *
 * TODO: Another name for this might make more sense as "adapter component".
 */
struct Barcode {
  u32 id;
  std::string sequence;
  std::string name;

  Barcode();
  Barcode(u32 id, std::string sequence, std::string name);

  auto operator<=>(const Barcode&) const = default;
};

/// @brief Generic collection of adapter-component sequences used by the LUT/matcher layer (SIDs, UMIs, spacers, etc.).
using BarcodePool = std::vector<Barcode>;

/// @brief Domain alias marking a BarcodePool that specifically holds sample-identifying barcodes for demux and metrics.
using SidPool = BarcodePool;

/**
 * @class MatchInfo
 * @brief Represents barcode match information and properties.
 *
 * The MatchInfo class encapsulates information about barcode matches, including the matched
 * barcodes, their positions, match type, and related properties. It allows for updating matches
 * and querying properties of the match information.
 */
class MatchInfo {
 public:
  /**
   * @brief Static method to retrieve the barcode ID from a MatchInfo instance.
   *
   * This static method extracts the barcode ID from the provided MatchInfo instance.
   * If the MatchInfo instance is empty (std::nullopt), returns an empty optional.
   *
   * @param match_info The MatchInfo instance to retrieve the barcode ID from.
   * @return An optional containing the barcode ID if available, otherwise an empty optional.
   */
  static std::optional<u32> BarcodeId(const std::optional<MatchInfo>& match_info);

  u32 BarcodeId() const {
    return _match.barcode_id;
  }

  /**
   * @brief Static method to retrieve the edit distance from a MatchInfo instance.
   *
   * This static method extracts the edit distance from the provided MatchInfo instance.
   * If the MatchInfo instance is empty (std::nullopt), returns an empty optional.
   *
   * @param match_info The MatchInfo instance to retrieve the edit distance from.
   * @return An optional containing the edit distance if available, otherwise an empty optional.
   */
  static std::optional<u32> EDist(const std::optional<MatchInfo>& match_info);

  u32 EDist() const {
    return _match.edist;
  }

  /**
   * @brief Static method to retrieve the position from a MatchInfo instance.
   *
   * This static method extracts the position from the provided MatchInfo instance.
   * If the MatchInfo instance is empty (std::nullopt), returns an empty optional.
   *
   * @param match_info The MatchInfo instance to retrieve the position from.
   * @return An optional containing the position if available, otherwise an empty optional.
   */
  static std::optional<LociRange> Position(const std::optional<MatchInfo>& match_info);

  /**
   * @brief Constructs a MatchInfo instance with provided match information.
   *
   * Initializes a MatchInfo object with the given barcode matches, positions, and match type.
   *
   * @param match The BarcodeMatch instance representing barcode match.
   * @param position The Position instances representing match position.
   * @param match_type The MatchType enum representing the type of the barcode match.
   */
  MatchInfo(const BarcodeMatch& match, const LociRange& position, MatchType match_type);

  /**
   * @brief Initialize a MatchInfo instance with default ("unknown yet") information.
   *
   */
  MatchInfo() = default;

  /**
   * @brief Updates the MatchInfo instance with new barcode match information.
   *
   * Updates the MatchInfo instance with new barcode match information, considering the provided
   * new match, new match type, start position (spos), and end position (epos). It compares the
   * new match to existing ones and decides whether to update or add matches based on various
   * match characteristics.
   *
   * @param new_match The new barcode match information to be considered for update or addition.
   * @param new_match_type The new match type associated with the new matches.
   * @param spos The start position of the new matches.
   * @param epos The end position of the new matches.
   */
  void Update(const BarcodeMatch& new_match, MatchType new_match_type, u32 spos, u32 epos);

  u32 SPos() const {
    return _position.spos;
  }

  u32 EPos() const {
    return _position.epos;
  }

  u32 Length() const {
    return _position.Length();
  }

  bool IsUnknown() const {
    return _match_type == MatchType::kUnknown;
  }

  bool IsAmbiguous() const {
    return _match_type == MatchType::kAmbiguous;
  }

  bool IsExact() const {
    return _match_type == MatchType::kExact;
  }

  bool HasMatches() const {
    return _match.edist != kUninitialized;
  }

  const LociRange& Position() const {
    return _position;
  }

  /**
   * @brief Compute a log-likelihood alignment score from the LUT match metadata.
   *
   * Uses the observed match window length (Length()) and edit distance (EDist()) to decompose
   * the alignment into insertions, deletions, substitutions, and matches. The net indel count is
   * computed as Length() - gt_len: a positive value represents insertions and a negative value
   * represents deletions. After accounting for the absolute net indel count, the remaining edits
   * are treated as substitutions, and the remaining ground-truth positions are treated as matches.
   * The final score is then computed by applying the scoring constants from the scoring namespace.
   *
   * @param gt_len Ground-truth barcode length (e.g. from BarcodePool[BarcodeId()].sequence.size())
   * @return Log-likelihood alignment score
   */
  s32 LogLikelihoodScore(const s32 gt_len) const {
    // positive = net insertions, negative = net deletions
    const s32 net_indels = static_cast<s32>(Length()) - gt_len;
    const s32 abs_net = net_indels < 0 ? -net_indels : net_indels;

    const s32 num_insertions = net_indels > 0 ? net_indels : 0;
    const s32 num_deletions = net_indels < 0 ? -net_indels : 0;
    // Remaining edits after accounting for length-changing operations must be substitutions.
    const s32 num_substitutions = static_cast<s32>(EDist()) - abs_net;
    // Matched bases = reference positions that were neither substituted nor deleted.
    const s32 num_matches = gt_len - num_substitutions - num_deletions;

    return (num_matches * scoring::kMatch) + (num_substitutions * scoring::kSubstitution) +
           (num_insertions * scoring::kInsertion) + (num_deletions * scoring::kDeletion);
  }

 private:
  BarcodeMatch _match;
  LociRange _position;
  MatchType _match_type{MatchType::kUnknown};
};
}  // namespace xoos::demux
