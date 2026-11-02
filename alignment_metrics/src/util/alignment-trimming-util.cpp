#include "util/alignment-trimming-util.h"

#include <cstddef>

#include <htslib/sam.h>

#include <xoos/types/int.h>
#include <xoos/util/math.h>

#include "core/alignment.h"

namespace xoos::alignment_metrics {

TrimmingResult TrimLeading(const Alignment& alignment, const u32 bases_to_trim_arg) {
  const bam1_t* const alignment_ptr = alignment.alignment_ptr;
  size_t first_cigar_index = 0;
  u32 bases_trimmed_from_start_of_first_cigar = 0;
  u32 bases_to_trim = bases_to_trim_arg;
  // Local copies seeded from the alignment; mutated as we walk the CIGAR and returned for the caller to apply.
  s64 rpos = alignment.rpos;
  u32 qpos = alignment.qpos;
  u32 read_length = alignment.read_length;
  u32 read_length_without_softclips = alignment.read_length_without_softclips;
  for (size_t i = 0; i < alignment_ptr->core.n_cigar; ++i) {
    const u32 cigar = bam_get_cigar(alignment_ptr)[i];
    const u32 cigar_op = bam_cigar_op(cigar);
    const u32 cigar_len = bam_cigar_oplen(cigar);
    const bool consumes_query = bam_cigar_type(cigar_op) & 1;
    const bool consumes_reference = bam_cigar_type(cigar_op) & 2;
    if (consumes_query && cigar_len <= bases_to_trim) {
      // If the CIGAR operation consumes query bases and is shorter than to the remaining bases to trim,
      // we skip the entire operation
      if (consumes_reference) {
        rpos += cigar_len;
      }
      if (cigar_op != BAM_CSOFT_CLIP) {
        read_length_without_softclips -= cigar_len;
      }
      read_length -= cigar_len;
      qpos += cigar_len;
      first_cigar_index = i + 1;
      bases_trimmed_from_start_of_first_cigar = 0;
      bases_to_trim -= cigar_len;
      if (bases_to_trim == 0) {
        break;
      }
    } else if (consumes_query && cigar_len > bases_to_trim) {
      // If the CIGAR operation consumes query bases and is longer than the remaining bases to trim,
      // it means we have reached the first CIGAR operation that will remain after trimming
      if (consumes_reference) {
        rpos += bases_to_trim;
      }
      if (cigar_op != BAM_CSOFT_CLIP) {
        read_length_without_softclips -= bases_to_trim;
      }
      read_length -= bases_to_trim;
      qpos += bases_to_trim;
      first_cigar_index = i;
      bases_trimmed_from_start_of_first_cigar = bases_to_trim;
      bases_to_trim = 0;
      break;
    } else {
      // If the CIGAR operation does not consume query bases, we do not trim anything and just move to the next
      // operation
      if (consumes_reference) {
        rpos += cigar_len;
      }
      first_cigar_index = i;
      bases_trimmed_from_start_of_first_cigar = 0;
    }
  }
  return {.boundary_cigar_index = first_cigar_index,
          .bases_trimmed_from_boundary_cigar = bases_trimmed_from_start_of_first_cigar,
          .new_rpos = rpos,
          .new_qpos = qpos,
          .new_read_length = read_length,
          .new_read_length_without_softclips = read_length_without_softclips};
}

TrimmingResult TrimTrailing(const Alignment& alignment, const u32 bases_to_trim_arg) {
  const bam1_t* const alignment_ptr = alignment.alignment_ptr;
  size_t last_cigar_index = math::SatSub(alignment_ptr->core.n_cigar, u32{1});
  u32 bases_trimmed_from_end_of_last_cigar = 0;
  u32 bases_to_trim = bases_to_trim_arg;
  // Local copies seeded from the alignment; trailing trim leaves rpos/qpos unchanged.
  u32 read_length = alignment.read_length;
  u32 read_length_without_softclips = alignment.read_length_without_softclips;
  for (s32 i = ToSigned(alignment_ptr->core.n_cigar) - 1; i >= 0; --i) {
    const u32 cigar = bam_get_cigar(alignment_ptr)[i];
    const u32 cigar_op = bam_cigar_op(cigar);
    const u32 cigar_len = bam_cigar_oplen(cigar);
    const bool consumes_query = bam_cigar_type(cigar_op) & 1;
    if (consumes_query && cigar_len <= bases_to_trim) {
      // If the CIGAR operation consumes query bases and is shorter than to the remaining bases to trim,
      // we skip the entire operation
      if (cigar_op != BAM_CSOFT_CLIP) {
        read_length_without_softclips -= cigar_len;
      }
      read_length -= cigar_len;
      last_cigar_index = math::SatSub(static_cast<size_t>(i), size_t{1});
      bases_trimmed_from_end_of_last_cigar = 0;
      bases_to_trim -= cigar_len;
      if (bases_to_trim == 0) {
        break;
      }
    } else if (consumes_query && cigar_len > bases_to_trim) {
      // If the CIGAR operation consumes query bases and is longer than the remaining bases to trim,
      // it means we have reached the last CIGAR operation that will remain after trimming
      if (cigar_op != BAM_CSOFT_CLIP) {
        read_length_without_softclips -= bases_to_trim;
      }
      read_length -= bases_to_trim;
      last_cigar_index = i;
      bases_trimmed_from_end_of_last_cigar = bases_to_trim;
      bases_to_trim = 0;
      break;
    } else {
      // If the CIGAR operation does not consume query bases, we do not trim anything and just move to the next
      // operation
      last_cigar_index = i;
      bases_trimmed_from_end_of_last_cigar = 0;
    }
  }
  return {.boundary_cigar_index = last_cigar_index,
          .bases_trimmed_from_boundary_cigar = bases_trimmed_from_end_of_last_cigar,
          .new_rpos = alignment.rpos,
          .new_qpos = alignment.qpos,
          .new_read_length = read_length,
          .new_read_length_without_softclips = read_length_without_softclips};
}

}  // namespace xoos::alignment_metrics
