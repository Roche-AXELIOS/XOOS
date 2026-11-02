#pragma once

#include <memory>
#include <string_view>
#include <vector>

#include <xoos/types/int.h>

#include "seq-lut.h"
#include "sequence/sequence-enums.h"

namespace xoos::demux {

// optimization by converting the character sequence to a 2-bit representation.
// this is a sequence of 16K bases, much more than we encounter
constexpr u16 kMaxMemorySequence{4096};

/**
 * @class SequenceTwoBit
 * @brief A class for 2-bit representation of sequences
 *
 * This class converts the incoming sequence to a 2-bit representation, which can be converted to hash values
 * much more efficiently. To avoid costly memory allocations, data is stored on the stack; we reserve 4K of memory
 * which should be more than plenty for short reads.
 */
class SequenceTwoBit {
  static constexpr u8 kOffset{64};  // add some data before the actual 2-bit representation to allow code optimization
 public:
  explicit SequenceTwoBit(const std::string_view& seq);

  size_t Length() const {
    return _length;
  }

  const u8* Data() const {
    return _two_bit_data + kOffset;
  }

 private:
  size_t _length;
  u8 _two_bit_data[kMaxMemorySequence];
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
  u64 mask64{0ULL};
  s32 spos{0};
  s32 epos{0};
  s32 length{0};
  s32 skip{0};

  Loci(const u64 mask, const s32 start, const s32 end, const s32 l) : mask64(mask), spos(start), epos(end), length(l) {
  }

  Loci(const s32 start, const s32 end) : spos(start), epos(end) {
  }
};

class SeqMatcher {
 public:
  using ExcisionLoci = std::vector<Loci>;
  static ExcisionLoci CreateRelativeExcisionLoci(s32 gt_seq_len,
                                                 s32 max_edist,
                                                 s32 max_wiggle_left,
                                                 s32 max_wiggle_right);
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
  SeqMatcher(u32 seq_len, s32 max_edist, s32 max_wiggle_left, s32 max_wiggle_right, SeqLutPtr lut);

  /**
   * @brief Finds barcode matches for a given read sequence.
   *
   * Searches for barcode matches in the sequence lookup table (LUT) based on the provided
   * read end, start position, and read sequence. Considers various potential start and end
   * positions for barcode matching, updating match information based on found matches.
   *
   * @param read_end The end of the read to match barcode against (k3p or k5p).
   * @param start_pos The start position of the barcode matching.
   * @param two_bit The encoded read sequence to match against barcodes.
   * @param seq_length The length of the read sequence.
   * @return A MatchInfo structure containing match details and type.
   */
  MatchInfo FindBarcode(ReadEnd read_end, u32 start_pos, const u8* two_bit, size_t seq_length) const;

  /**
   * @brief Greedy scan for the next barcode match(es) from a scan position.
   *
   * Slides a barcode-length window across the remaining sequence (from `pos` to the end for k5p,
   * from `pos` backward for k3p) and populates `results` with the first match(es) found. Unlike
   * FindBarcode, this does NOT use the precomputed excision loci or wiggle parameters — it scans
   * the entire remaining sequence. The prefilter LUT is still used for efficiency.
   *
   * For k5p the return position is EPos (right edge); for k3p it is SPos (left edge).
   * All matches sharing the same return position are returned (rather than collapsing to "ambiguous").
   *
   * @param read_end     Search direction (k5p = left-to-right, k3p = right-to-left).
   * @param pos          Start scan position.
   * @param two_bit      2-bit encoded read sequence.
   * @param seq_length   Length of the read in bases.
   * @param[out] results Cleared and filled with all matches at the first return position (empty if no match).
   */
  void FindNextBarcode(
      ReadEnd read_end, s32 pos, const u8* two_bit, size_t seq_length, std::vector<MatchInfo>& results) const;

  const BarcodePool& Pool() const;

  // Allow caller to peek into the LUT
  const auto& Lut() const {
    return _lut;
  }

  // Common operation for easy access of front sequence
  const auto& GetFrontSeq() const {
    return _lut->Pool().front().sequence;
  }

  u32 GetMaxDist() const {
    return _max_edist;
  }

 private:
  /**< The length of the ground truth sequence. */
  u32 _seq_len;
  /**< Maximum edit distance for matching. */
  u32 _max_edist;
  ExcisionLoci _relative_excision_loci;
  /**< The sequence lookup table for barcode matching. */
  SeqLutPtr _lut;
  const size_t _nr_loci;
};

using SeqMatcherPtr = std::shared_ptr<SeqMatcher>;
}  // namespace xoos::demux
