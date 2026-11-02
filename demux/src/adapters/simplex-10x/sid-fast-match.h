#pragma once

#include <optional>
#include <string>
#include <string_view>

#include <gtl/phmap.hpp>

#include <xoos/types/int.h>

#include "sequence/matcher/match-info.h"

namespace xoos::demux {

/**
 * @brief Fast-path SID assignment using 2-bit encoded hash table.
 *
 * Checks for an exact runway match at read index 0 or 1, then encodes the
 * subsequent (SID + stem_prefix) bases as a 2-bit-per-base uint32_t key and
 * looks it up in a hash table that tolerates up to 1 substitution in the SID
 * portion. The stem prefix must match exactly.
 *
 * With stem_prefix_len=4: key = SID(12bp) + stem(4bp) = 16bp = 32 bits.
 */
class SidFastMatch {
 public:
  struct Result {
    u32 sid_id;
    u32 edist;
    // Position right after the SID (stem is left in the insert).
    u32 insert_start;
  };

  /**
   * @brief Build the fast-path matcher.
   *
   * Returns nullopt if the SID pool is unsuitable (variable-length SIDs,
   * SID > 16bp, non-ACGT bases, or hash-key collisions between distinct
   * SIDs at edit distance 0). The stem prefix length (kStemPrefixLen) is
   * silently truncated when it would exceed the encodable key length or
   * the actual stem length.
   *
   * @param runway   Full runway sequence (e.g. "CAACAA"). TryMatch checks
   *                 for the truncated runway (runway[1:]) at offset 0
   *                 first, then the full runway at offset 0.
   * @param sid_pool SID barcode pool.
   * @param stem     Full stem sequence.
   */
  static std::optional<SidFastMatch> Build(std::string_view runway, const BarcodePool& sid_pool, std::string_view stem);

  /**
   * @brief Try the fast-path on a raw read sequence.
   *
   * Checks for the truncated runway (runway[1:]) at offset 0 first,
   * then the full runway at offset 0.
   */
  std::optional<Result> TryMatch(const char* seq, size_t length) const;

  /**
   * @brief Encode a DNA sequence to a 2-bit integer.
   * @return The encoded value, or nullopt if any base is not in {A, C, G, T}.
   */
  static std::optional<u32> Encode(const char* seq, u32 len);

 private:
  // Sentinel SID ID for hash entries where two different SIDs both have a
  // 1-mismatch variant producing the same key. These reads fall through
  // to the slow path.
  static constexpr u32 kAmbiguous = std::numeric_limits<u32>::max();

  /** @brief 2-bit encoding: A=0, C=1, G=2, T=3. Bits per DNA base. */
  static constexpr u32 kBitsPerBase = 2;
  /** @brief Max bases encodable in a u32 key. */
  static constexpr u32 kMaxEncodableLen = 32 / kBitsPerBase;
  /** @brief Highest valid 2-bit base value (3 = T). */
  static constexpr u8 kMaxValidBase = (1U << kBitsPerBase) - 1;
  /** @brief 4bp stem prefix + 12bp SID = 16bp = fills the u32 key. */
  static constexpr u32 kStemPrefixLen = 4;

  struct SidEntry {
    u32 sid_id;
    u32 edist;
  };

  /**
   * @brief Validate that the SID pool is suitable for the fast-path.
   * @return The uniform SID length, or nullopt if validation fails.
   */
  static std::optional<u32> ValidatePool(const BarcodePool& sid_pool);

  /**
   * @brief Insert all 1-substitution variants for each SID into the hash table.
   *
   * Ambiguous collisions are marked kAmbiguous; exact-match entries are preserved.
   */
  static void InsertMismatchVariants(const BarcodePool& sid_pool,
                                     u32 sid_len,
                                     u32 stem_bits,
                                     gtl::flat_hash_map<u32, SidEntry>& sid_map);

  SidFastMatch(std::string_view runway, u32 sid_len, u32 key_len, gtl::flat_hash_map<u32, SidEntry> sid_map);

  std::string _runway;
  u32 _sid_len;
  // _sid_len + effective stem_prefix_len.
  u32 _key_len;
  gtl::flat_hash_map<u32, SidEntry> _sid_map;
};

}  // namespace xoos::demux
