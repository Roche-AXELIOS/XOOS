#pragma once

#include <string>

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/yc-decode/yc-decoder.h>

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
char GetBase(const u8* seq, u64 qpos);

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
std::string GetSequence(const u8* seq, u64 start, u64 len);

/**
 * Get the base qualities from a read.
 *
 * @param quals The base qualities of a read sequence.
 * @param start The start position of the qualities.
 * @param len The length of the qualities to fetch.
 *
 * @return The base qualities of a sequence from position start (inclusive) to start + len (exclusive).
 *
 * @warning This function does not perform bounds checking. The behavior is undefined if start + len is greater than the
 * read length.
 */
vec<u8> GetQualities(const u8* quals, u32 start, u32 len);

}  // namespace xoos::read_collapser
