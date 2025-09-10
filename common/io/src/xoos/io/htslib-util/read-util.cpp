#include "read-util.h"

#include <string>

#include <htslib/sam.h>

#include <xoos/types/int.h>

namespace xoos::io {

/**
 * Get the base at a given position in a read.
 *
 * @param seq The read sequence.
 * @param qpos The position in the read.
 *
 * @return The base at position `qpos`.
 *
 * @warning This function does not perform bounds checking. The behavior is undefined if `qpos` is greater than or equal
 * to the read length.
 */
char GetBase(const u8* seq, const u32 qpos) {
  // Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G, 8 for T and 15 for N.
  // Two bases are packed in one byte with the base at the higher 4 bits having smaller coordinate on the read.
  const u8 base = bam_seqi(seq, qpos);
  return seq_nt16_str[base];
}

/**
 * Get the sequence from a read.
 * This method is provided to as a wrapper around bam_seqi to get the sequence from a read.
 * It is used primarily for a limited range, such as a portion of reference sequence.
 *
 * @warning This function does not perform bounds checking. The behavior is undefined if start + len is greater than the
 * read length.
 *
 * @param seq The read sequence.
 * @param start The start position of the sequence.
 * @param len The length of the sequence to fetch.
 *
 * @return The sequence containing bases from position start (inclusive) to start + len (exclusive).
 */
std::string GetSequence(const u8* seq, const u32 start, const u32 len) {
  std::string sequence;
  sequence.reserve(len);
  for (u32 i = start; i < start + len; ++i) {
    sequence.push_back(GetBase(seq, i));
  }
  return sequence;
}

/**
 *  This method removes the soft clipped bases from sequence length from the first or last cigar operation.
 */
u32 ApproximateGenomicInsertLength(const u32 l_seq, const u32* cigar, const u32 n_cigar_op) {
  auto read_length = l_seq;
  // test and remove the first soft-clipped base
  if (n_cigar_op > 0 && bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP) {
    read_length -= bam_cigar_oplen(cigar[0]);
  }
  // test and remove the last soft-clipped base
  if (n_cigar_op > 1 && bam_cigar_op(cigar[n_cigar_op - 1]) == BAM_CSOFT_CLIP) {
    read_length -= bam_cigar_oplen(cigar[n_cigar_op - 1]);
  }
  return read_length;
}

}  // namespace xoos::io
