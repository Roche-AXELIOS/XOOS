#pragma once

#include <cstddef>

#include <xoos/types/int.h>

namespace xoos::alignment_metrics {

class Alignment;

/**
 * Result of trimming query bases from one end of an alignment. Reports the index of the first (for leading trim)
 * or last (for trailing trim) CIGAR operation that remains after trimming, how many of its bases were trimmed away,
 * and the resulting position/length fields after trimming.
 */
struct TrimmingResult {
  // Index of the first (for leading trim) or
  // last (for trailing trim) CIGAR operation that remains after trimming.
  size_t boundary_cigar_index{};
  // Number of bases trimmed from the first (for leading trim) or
  // last (for trailing trim) CIGAR operation that remains after trimming.
  u32 bases_trimmed_from_boundary_cigar{};
  s64 new_rpos{};
  u32 new_qpos{};
  u32 new_read_length{};
  u32 new_read_length_without_softclips{};
};

/**
 * Trim `bases_to_trim` query bases from the start of the read, walking the CIGAR operations forward.
 * Returns the index of the first CIGAR operation that remains after trimming, how many of its bases were
 * trimmed, and the resulting rpos, qpos, read_length, and read_length_without_softclips.
 */
TrimmingResult TrimLeading(const Alignment& alignment, u32 bases_to_trim);

/**
 * Trim `bases_to_trim` query bases from the end of the read, walking the CIGAR operations backward.
 * Returns the index of the last CIGAR operation that remains after trimming, how many of its bases were
 * trimmed, and the resulting read_length and read_length_without_softclips (rpos and qpos are passed through
 * unchanged).
 */
TrimmingResult TrimTrailing(const Alignment& alignment, u32 bases_to_trim);

}  // namespace xoos::alignment_metrics
