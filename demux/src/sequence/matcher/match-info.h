#pragma once

#include <optional>
#include <string>
#include <vector>

#include "sequence/loci-range.h"

namespace xoos::demux {
enum class MatchType {
  kExact,
  kUnknown,
  kAmbiguous,
  kInExact,
};

// KZ2024 - to minimize memory usage, we'd like to limit the amount of memory used by this structure to 16 bits.
// This would then allow us to put all the barcode matches in a very large 8 GB (2^32) lookup table.
constexpr static uint kUninitialized{3};

struct BarcodeMatch {
  enum BarcodeType { kSID5p = 0, kUMI5p, kSpacer5p, kRunway5p, kSID3p, kUMI3p, kSpacer3p };

  BarcodeMatch() : barcode_id(0), edist(kUninitialized), barcode_type(0) {}

  BarcodeMatch(uint id, uint dist, uint type) : barcode_id(id), edist(dist), barcode_type(type) {}

  BarcodeMatch(uint id, uint dist) : barcode_id(id), edist(dist), barcode_type(0) {}

  uint16_t barcode_id : 11;   // Max 2048 IDs
  uint16_t edist : 2;         // 0..2, the value 3 indicates uninitialized
  uint16_t barcode_type : 3;  // future extension for umi3p, umi5p, sid3p, sid5p, ... (up to 8 types)
};

struct Barcode {
  uint id;
  std::string sequence;
  std::string name;

  Barcode();
  Barcode(uint id, std::string sequence, std::string name);

  auto operator<=>(const Barcode&) const = default;
};

using BarcodePool = std::vector<Barcode>;

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
  static std::optional<uint> BarcodeId(const std::optional<MatchInfo>& match_info);

  /**
   * @brief Static method to retrieve the edit distance from a MatchInfo instance.
   *
   * This static method extracts the edit distance from the provided MatchInfo instance.
   * If the MatchInfo instance is empty (std::nullopt), returns an empty optional.
   *
   * @param match_info The MatchInfo instance to retrieve the edit distance from.
   * @return An optional containing the edit distance if available, otherwise an empty optional.
   */
  static std::optional<uint> EDist(const std::optional<MatchInfo>& match_info);

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

 public:
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
  void Update(const BarcodeMatch& new_match, MatchType new_match_type, uint spos, uint epos);

  uint SPos() const { return _position.spos; }

  uint EPos() const { return _position.epos; }

  uint Length() const { return _position.Length(); }

  uint EDist() const { return _match.edist; }

  uint BarcodeId() const { return _match.barcode_id; }

  bool IsUnknown() const { return _match_type == MatchType::kUnknown; }

  bool IsAmbiguous() const { return _match_type == MatchType::kAmbiguous; }

  bool IsExact() const { return _match_type == MatchType::kExact; }

  bool HasMatches() const { return _match.edist != kUninitialized; }

  const LociRange& Position() const { return _position; }

 private:
  BarcodeMatch _match;
  LociRange _position;
  MatchType _match_type{MatchType::kUnknown};
};
}  // namespace xoos::demux
