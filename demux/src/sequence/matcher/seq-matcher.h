#pragma once

#include <memory>
#include <vector>

#include "seq-lut.h"

// KZ2024 - speed optimization by converting the character sequence to a 2-bit representation.
const int kMaxMemorySequence{4096};  // this is a sequence of 16K bases, much more than we encounter

namespace xoos::demux {
/**
 * @class SequenceTwoBit
 * @brief A class for 2-bit representation of sequences
 *
 * This class converts the incoming sequence to a 2-bit representation, which can be converted to hash values
 * much more efficiently. To avoid costly memory allocations, data is stored on the stack; we reserve 4K of memory
 * which should be more than plenty for short reads.
 */
class SequenceTwoBit {
  static constexpr int kOffset{64};  // add some data before the actual 2-bit representation to allow code optimization
 public:
  explicit SequenceTwoBit(const std::string_view& seq);

  size_t Length() const { return _length; }

  const uint8_t* Data() const { return _two_bit_data + kOffset; }

 private:
  size_t _length;
  uint8_t _two_bit_data[kMaxMemorySequence];
};

enum class ReadEnd {
  k5p,
  k3p,
};

/**
 * @class SeqMatcher
 * @brief A class for sequence matching and excision locus generation.
 *
 * This class provides functionality for sequence matching and generating excision loci
 * based on specified parameters. It assists in identifying and processing sequence excisions
 * relative to a ground truth sequence, considering edit distance and wiggle room.
 */

struct Loci {
  uint64_t mask64{0ul};
  int spos{0};
  int epos{0};
  int length{0};
  int skip{0};

  Loci(uint64_t mask, int start, int end, int l) {
    mask64 = mask;
    spos = start;
    epos = end;
    length = l;
  }

  Loci(int start, int end) {
    spos = start;
    epos = end;
  }
};

class SeqMatcher {
 public:
  using ExcisionLoci = std::vector<Loci>;
  static ExcisionLoci CreateRelativeExcisionLoci(int gt_seq_len, int max_edist, int max_wiggle_left,
                                                 int max_wiggle_right);
  /**
   * @brief Constructs a SeqMatcher object with provided parameters and LUT.
   *
   * Initializes a SeqMatcher object with the given sequence length, maximum edit distance,
   * maximum wiggle distances, and the provided sequence lookup table (LUT). Generates relative
   * excision loci based on the provided sequence length and parameters.
   *
   * @param seq_len The length of the sequences being processed.
   * @param max_edist The maximum allowed edit distance for excisions.
   * @param max_wiggle_left The maximum allowed left wiggle distance for excisions.
   * @param max_wiggle_right The maximum allowed right wiggle distance for excisions.
   * @param lut The sequence lookup table for barcode matching.
   */
  SeqMatcher(uint seq_len, int max_edist, int max_wiggle_left, int max_wiggle_right, SeqLutPtr lut);

  /**
   * @brief Finds barcode matches for a given read sequence.
   *
   * Searches for barcode matches in the sequence lookup table (LUT) based on the provided
   * read end, start position, and read sequence. Considers various potential start and end
   * positions for barcode matching, updating match information based on found matches.
   *
   * @param read_end The end of the read to match barcode against (k3p or k5p).
   * @param start_pos The start position of the barcode matching.
   * @param read_seq The read sequence to match against barcodes.
   * @return A MatchInfo structure containing match details and type.
   */
  MatchInfo FindBarcode(ReadEnd read_end, uint start_pos, const uint8_t* two_bit, size_t length) const;
  const BarcodePool& Pool() const;

  // Allow caller to peek into the LUT
  const auto& Lut() const { return _lut; }

 private:
  uint _seq_len; /**< The length of the ground truth sequence. */
  ExcisionLoci _relative_excision_loci;
  SeqLutPtr _lut; /**< The sequence lookup table for barcode matching. */
  const size_t _nr_loci;
};

using SeqMatcherPtr = std::shared_ptr<SeqMatcher>;
}  // namespace xoos::demux
