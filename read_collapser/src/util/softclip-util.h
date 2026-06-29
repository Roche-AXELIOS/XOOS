#pragma once

#include <string>

#include <xoos/types/vec.h>

#include "io/alignment.h"

namespace xoos::read_collapser {

/**
 * @brief Given a list of alignments in a cluster, find the leftmost reference position.
 * If trim_overhangs is true, the leftmost position is adjusted to account for any soft-clipped
 * bases at the start of reads in the cluster so that all reads have "blunt" left ends.
 */
u32 LeftmostRefPos(const vec<AlignmentPtr>& reads_in_cluster, bool trim_overhangs);

/**
 * @brief Given a list of alignments in a cluster, find the rightmost reference position.
 * If trim_overhangs is true, the rightmost position is adjusted to account for any soft-clipped
 * bases at the end of reads in the cluster so that all reads have "blunt" right ends.
 */
u32 RightmostRefPos(const vec<AlignmentPtr>& reads_in_cluster, bool trim_overhangs);

/**
 * @brief Given a vector of alignments in a cluster, extract the left soft-clipped sequences.
 *
 * NOTE: If reads in a cluster have different start positions, bases from the actual start of each read
 * up to (but not including) the rightmost boundary between soft-clips and aligned sequences will be
 * extracted as soft-clipped sequences.
 *
 * Example:
 * Read 1: 5'---ACGTACGT---3' (start=100, left soft-clip length=3)
 * Read 2: 5'------GACGT---3' (start=101, left soft-clip length=1)
 * Left soft-clipped sequences: ["ACGT", "G"]
 * Note that for Read 1, although the first aligned base 'T' at position 100 is not part of the soft-clip,
 * it is still handled as part of the left soft-clip sequence to ensure that the right end of the soft-clip
 * sequence and the left end of the aligned sequence can be stitched together correctly later.
 */
vec<std::string> GetLeftSoftclipSequences(const vec<AlignmentPtr>& reads_in_cluster);

/**
 * @brief Given a vector of alignments in a cluster, extract the right soft-clipped sequences.
 *
 * NOTE: If reads in a cluster have different end positions, bases from the actual end of each read
 * back to (but not including) the leftmost boundary between aligned sequences and soft-clips will
 * be extracted as soft-clipped sequences.
 *
 * Example:
 * Read 1: 5'---ACGTACGT---3' (end=108, right soft-clip length=3)
 * Read 2: 5'---ACGTC------3' (end=107, right soft-clip length=1)
 * Right soft-clipped sequences: ["ACGT", "C"]
 * Note that for Read 1, although the last aligned base 'A' at position 107 is not part of the soft-clip,
 * it is still handled as part of the right soft-clip sequence to ensure the left end of the soft-clip
 * sequence and the right end of the aligned sequence can be stitched together correctly later.
 */
vec<std::string> GetRightSoftclipSequences(const vec<AlignmentPtr>& reads_in_cluster);

}  // namespace xoos::read_collapser
