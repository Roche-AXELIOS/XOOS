#pragma once

#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <gtl/phmap.hpp>

#include "match-info.h"

namespace xoos::demux {
/**
 * @class IntLut
 * @brief A class that represents a sequence lookup table for barcode matching.
 *
 * This class provides functionality to manage a sequence lookup table (LUT) for barcode matching.
 * Extensive performance analysis showed that unordered_map<uint64, int> can be much faster compared to
 * unordered_map<string, int>  (as we can optimize the hash function ourselves), so we replaced the original
 * approach of one unordered_map<string,int> with vector<unordered_map<uint64,int> (we need multiple maps to be
 * able to support different lengths).
 *
 * In additional to the unordered_maps (which are relatively slow compared to a simple memory lookup), this class
 * maintains an additional lookup table that entirely reside in memory.
 * The lookup table has a size of 24 bits  (that's 2M) and therefore will fit (almost) in processor
 * caches; this size was optimized based on benchmarks. This LUT is referred to as "the pre-filter" and tests
 * whether up to 12 bases in the sequence might match a barcode; most incoming sequences are immediately rejected,
 * which is very fast as we usually will have cache hits. If we find a potential candidate,
 *
 */

const uint8_t kDefaultMaxLength = 4;

class IntLut {
 public:
  // Helper array: Unique encoding of three bit value into non-overlapping bits
  constexpr static int kSetBits[] = {1, 2, 4, 8, 16, 32, 64, 128};

  // Constructor - requires maximum length expected for barcode
  //  explicit IntLut(uint max_length = kDefaultMaxLength);
  explicit IntLut(uint max_length);

  // Used only in case pre-filter and filter LUTs indicate that there might be a match.
  inline std::tuple<BarcodeMatch, MatchType> operator()(uint length, uint64_t hash_val) const {
    const auto it = _hash_tables[length].find(hash_val);  // look for match in hash table with appropriate length

    if (it == _hash_tables[length].end()) {
      return {BarcodeMatch(), MatchType::kUnknown};  // Nada
    }
    // Return the match. Note: we do not support ambiguous matches anymore (they did not occur), we assume only one
    // match
    const auto& match = it->second;
    return {match, match.edist == 0 ? MatchType::kExact : MatchType::kInExact};
  }

  // This function adds a barcode with specified length, hash value and barcode to our 'dictionary'.
  void Add(uint length, uint64_t hash_val, const BarcodeMatch& barcode);

  // Return const reference to our internal representation of the pre-filter LUT - required for speed optimization.
  const auto& PrefilterValues() const { return _prefilter; }

  int PrefilterMask() const { return _prefilter_mask; }

  const auto& HashTables() const { return _hash_tables; }

 private:
  std::vector<gtl::flat_hash_map<uint64_t, BarcodeMatch>> _hash_tables;
  std::vector<uint8_t> _prefilter; /**< Contains set bit for every valid barcode entry (up to 11 bases). */
  int _prefilter_mask{0};          /**< A mask that identifies the relevant bits to be used for the prefilter. */
};

/**
 * @class SeqLut
 * @brief A class that represents a sequence lookup table for barcode matching.
 *
 * This class provides functionality to manage a sequence lookup table (LUT) for barcode matching.
 * It allows finding matches in the LUT based on input sequences and provides information about the
 * matches and their types.
 */
class SeqLut {
 public:
  using Lut = gtl::flat_hash_map<std::string, std::vector<BarcodeMatch>>;
  using FindResult = std::tuple<BarcodeMatch, MatchType>;
  /**
   * @brief Constructs a SeqLut object using provided LUT and BarcodePool.
   *
   * Initializes a SeqLut object with the given lookup table (LUT) and BarcodePool.
   *
   * @param lut A lookup table containing sequence-barcode matches.
   * @param pool A BarcodePool containing the associated barcode information.
   */
  SeqLut(const IntLut& lut, const BarcodePool& pool);
  ~SeqLut() = default;

  /**
   * @brief Finds matches for the given sequence in the LUT.
   *
   * Searches the LUT for matches to the provided sequence and determines the type of match
   * based on edit distance and match count. Returns a FindResult containing match information
   * and the corresponding match type.
   *
   * @param seq The input sequence to find matches for.
   * @return A FindResult containing information about matches and match type.
   */
  FindResult SimpleFind(const std::string_view& seq) const;
  /**
   * @brief Returns a reference to the associated BarcodePool.
   *
   * Provides access to the BarcodePool associated with this SeqLut object.
   *
   * @return A reference to the associated BarcodePool.
   */

  FindResult SimpleFind(const std::string_view& seq, uint spos, uint epos) const;

  /**
   * @brief Returns a reference to the associated BarcodePool.
   *
   * Provides access to the BarcodePool associated with this SeqLut object. This method
   * is using an integer-based hash table, which should be faster than using strings.
   *
   * @return A reference to the associated BarcodePool.
   */

  inline FindResult SimpleFind(uint length, uint64_t hash_value) const { return _lut(length, hash_value); }

  /**
   * @brief Returns a match in the hash table for the specified hash value and length.
   *
   * @return The match value.
   */

  const BarcodePool& Pool() const;

  /**
   * @brief Calculates a 64-bit hash value using the specified string, position, and length.
   *
   * This function calculates hash values for substrings up to a length of 32 (32 x 2 bits = 64 bits).
   * No check is done on the value of the specified length (for performance reasons).
   *
   * @return A 64-bit hash value.
   */
  static uint64_t HashValue(const std::string_view& seq, uint start, uint length);

  /**
   * @brief converts a 2-bit representation to a string with specified start and length position.
   *
   * This function is a performance-optimized (no error checking) converter that converts a 2-bit representation
   * to a string (8-bits).
   *
   * Make sure that the p_input destination memory is large enough to receive the string representation;
   * no error checking is done on that.
   */
  static std::string ConvertFrom2Bit(const uint8_t* p_input, uint length);

  const auto& PrefilterValues() const { return _lut.PrefilterValues(); }

  auto PrefilterMask() const { return _lut.PrefilterMask(); }

  const auto& HashTables() const { return _lut.HashTables(); }

 private:
  const IntLut _lut;       /**< Lookup table using integer hash values, one per length. */
  const BarcodePool _pool; /**< The associated BarcodePool containing barcode information. */
};

using SeqLutPtr = std::shared_ptr<SeqLut>;

}  // namespace xoos::demux
