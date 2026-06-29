#include "util/read-util.h"

#include <string>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <xoos/error/error.h>
#include <xoos/types/int.h>
#include <xoos/util/math.h>

namespace xoos::read_collapser {

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
char GetBase(const u8* seq, const u64 qpos) {
  // Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G, 8 for T and 15 for N.
  // Two bases are packed in one byte with the base at the higher 4 bits having smaller coordinate on the read.
  const u8 base = bam_seqi(seq, qpos);
  return seq_nt16_str[base];
}

/**
 * Get the sequence from a read.
 *
 * @param seq The read sequence.
 * @param start The start position of the sequence.
 * @param len The length of the sequence to fetch.
 *
 * @return The sequence containing bases from position start (inclusive) to start + len (exclusive).
 *
 * @warning This function does not perform bounds checking. The behavior is undefined if start + len is greater than the
 * read length.
 */
std::string GetSequence(const u8* seq, const u64 start, const u64 len) {
  std::string sequence;
  sequence.reserve(len);
  for (u64 i = start; i < start + len; ++i) {
    sequence.push_back(GetBase(seq, i));
  }
  return sequence;
}

/**
 * Get the base qualities from a read.
 *
 * @param quals
 * @param start The start position of the qualities.
 * @param len The length of the qualities to fetch.
 *
 * @return The base qualities of a sequence from position start (inclusive) to start + len (exclusive).
 *
 * @warning This function does not perform bounds checking. The behavior is undefined if start + len is greater than the
 * read length.
 */
vec<u8> GetQualities(const u8* const quals, const u32 start, const u32 len) {
  return {quals + start, quals + start + len};
}

}  // namespace xoos::read_collapser
