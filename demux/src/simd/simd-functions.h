/**
 *  Functions that are internally using SIMD optimizations.
 */

#pragma once

#include <array>
#include <bit>

#include <xoos/types/int.h>

#include "adapters/duplex/duplex-match.h"

namespace simd {

/**
 * *brief Calculate the 2-bit encoding of a DNA sequence.
 * *param p_seq Pointer to the input sequence.
 * *param length Length of the input sequence.
 * *param p_output Pointer to the output sequence.
 * This function calculates a 2-bit representation for a specified part of the input sequence.
 * Every base is assigned 2 bits (A=00, C=01, G=10, T=11); this assignment is advantageous as
 * A/T and C/G are complementary (find the complement by flipping the bits).
 */
void ConvertTo2Bit(const xoos::u8* p_seq, xoos::u32 start, xoos::u32 length, xoos::u8* p_output);

/**
 * @brief Find the index of the first IsSpace() character in a string.
 * When parsing fastq files, the following loop was computationally expensive:
 *      for (i = this->begin; i < this->end; ++i) {
 *              if (IsSpace(this->buf[i])) break;
 *      }
 * This function uses SIMD instructions to find the first space character in a string.
 */
xoos::u64 FindFirstSpace(const xoos::u8* buf, xoos::u64 begin, xoos::u64 end);

/**
 * @brief Find the index of the first IsSpace() character in a string which is not ' '.
 * When parsing fastq files, the following loop was computationally expensive:
 *      for (i = this->begin; i < this->end; ++i) {
 *              if (IsSpace(this->buf[i])) break;
 *      }
 * This function uses SIMD instructions to find the first space character except ' ' in a string.
 */
xoos::u64 FindFirstNonTrivialSpace(const xoos::u8* buf, xoos::u64 begin, xoos::u64 end);

// old union approach used UB trick as mentioned in https://news.ycombinator.com/item?id=38472174
inline xoos::u64 Load64(const xoos::u8* const b) {
  static_assert(std::endian::native == std::endian::little, "This function assumes a little-endian architecture.");
  const std::array arr = {b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]};
  return std::bit_cast<xoos::u64>(arr);
}

/** @brief An optimized version of functions used by edlib (pairwise alignment)
 * This function converts a sequence into an alternate representation where each base is assigned an 8-bit integer
 * number (this functionality is performed by edlib's transformSequences()) and, additionally, create a representation
 * that consists of four binary masks (performed by edlib's buildPeq()).
 * @param buf Pointer to the input sequence
 * @param begin Start position in the string (note: start is assumed to be smaller than end)
 * @param end End position in the string (inclusive to avoid confusion with reversing strings).
 * @param digit In/Out digit representation of the sequence (sequences using values 0..3)
 * @param reverse Reverse the sequence (true/false)
 * @param invert Reverse complement the sequence (true/false)
 */

void TransformSequence(
    const char* buf, xoos::s32 begin, xoos::s32 end, xoos::u8* digit, bool reverse = false, bool invert = false);

/** @brief BuildPeq, an optimized version of functions used by edlib (pairwise alignment)
 * This function creates a representation
 * that consists of four binary masks (optimized version of edlib's buildPeq()).
 * @param buf Pointer to the input sequence
 * @param digit In digit representation of the sequence (sequences using values 0..3)
 * @param bit   In/Out bit representation of the sequence (four binary masks)
 */

void BuildPeq(const xoos::u8* digit, xoos::s32 length, xoos::u64* bit);

/**
 * @brief Copies matching bases from input to consensus sequence during alignment.
 *
 * This function identifies consecutive matching positions in an edlib alignment
 * (encoded as zero values) and copies the corresponding bases from the input
 * sequence to the consensus. The function processes up to 64 bases at a time
 * using SIMD instructions.
 *
 * @param alignment Pointer to the edlib alignment array (0 = match, non-zero = mismatch/gap)
 * @param bases_in Pointer to the input base sequence
 * @param consensus_out Pointer to the output consensus sequence buffer
 * @return Number of consecutive matching bases found at the start of the alignment
 *
 * @note The function may copy more bases than the returned count for efficiency;
 *       subsequent processing will overwrite any incorrect values.
 */
xoos::s32 CopyMatchingSeq(const xoos::u8* alignment, char* bases_in, char* consensus_out);

/**
 * @brief Sets quality values for matching positions during alignment.
 *
 * This function identifies consecutive matching positions in an edlib alignment
 * (encoded as zero values) and sets the corresponding quality values to a
 * specified base quality. The function processes up to 64 positions at a time
 * using SIMD instructions.
 *
 * @param alignment Pointer to the edlib alignment array (0 = match, non-zero = mismatch/gap)
 * @param quality_out Pointer to the output quality sequence buffer
 * @param base The quality character to assign to matching positions
 * @return Number of consecutive matching bases found at the start of the alignment
 *
 * @note The function may set more quality values than the returned count for efficiency;
 *       subsequent processing will overwrite any incorrect values.
 */
xoos::s32 CopyMatchingQual(const xoos::u8* alignment, char* quality_out, char base);

/* @brief Replacement for std::reverse_copy of 8-bit values, which turned out to be time consuming.
 * Raw pointers are required by the highway HWY_DYNAMIC_DISPATCH boundary; the [start, end) range
 * already carries the size, and callers pass fixed-size hot-path buffers. */
void ReverseCopy(const xoos::u8* const p_start, const xoos::u8* const p_end, xoos::u8* const p_out);  // NOSONAR(S945)

}  // namespace simd

namespace xoos::demux {
struct FixedReadRecord;
/**
 * This function is used in Duplex HD demuxing as the first step of processing; we pick a position close to the
 * 3' end of the read and search for the exact "reverse complementary" sequence at the 5p end of the read.
 * Assumptions: input data is in 2-bit format and input parameters are valid (no error checking within this function).
 * @param start_ptr Start address of the 2-bit sequence (base 0 at 5p end).
 * @param length Length of the 2-bit sequence (number of bases).
 * @return Position of the exact match (0-based) or -1 if no exact match was found.
 *
 * Note that we optimized this function to find either an exact match or a match with a single substitution; we do
 * pairwise comparison of 16 bases at a time and use a threshold of 0 (exact) or 1 (1 substitution OK). No attempts
 * at handling deletes/insertions were made.
 *
 */
s32 FindSymmetryPosition(FixedReadRecord& record);

bool FindMarker(
    const CascadedLUTs& luts, FixedReadRecord& record, const s64* offsets, const s64* masks, const s64* types);

// This function looks for the hairpin sequence in the read using an approach that resembles Igor Mandric's algorithm.
// It evaluates the popcnt() ("matching score") for 32 positions starting at the specified offset; it returns the
// maximum score that was found.
s64 FindHairpinSliding(FixedReadRecord& record, s64 offset);

}  // namespace xoos::demux
