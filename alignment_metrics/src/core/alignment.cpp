#include "core/alignment.h"

#include <span>

#include <htslib/sam.h>

#include <xoos/error/error.h>
#include <xoos/io/htslib-util/read-util.h>
#include <xoos/types/int.h>
#include <xoos/util/math.h>

#include "alignment-metrics-options.h"
#include "metadata/alignment-metadata.h"
#include "util/homopolymer-util.h"

namespace xoos::alignment_metrics {

Alignment::Alignment(bam1_t* const ptr, const AlignmentMetricsOptions& options)
    : alignment_ptr(ptr),
      alignment_metadata(CreateAlignmentMetadata(alignment_ptr, options.ignore_family)),
      // Skip HP masking and base type decoding if specified or if we are only calculating read metrics
      // since base types are not used for read metrics
      qualities(MaskReadHomopolymers(alignment_ptr,
                                     !NeedHpMasking(options),
                                     !NeedBaseTypeDecoding(options),
                                     options.base_quality_threshold_for_hp_masking)),
      cigars(bam_get_cigar(alignment_ptr), alignment_ptr->core.n_cigar),
      rpos(alignment_ptr->core.pos),
      qpos(0),
      reference_length(bam_cigar2rlen(ToSigned(alignment_ptr->core.n_cigar), bam_get_cigar(alignment_ptr))),
      read_length(alignment_ptr->core.l_qseq),
      read_length_without_softclips(io::GetReadLengthWithoutSoftClips(
          alignment_ptr->core.l_qseq, bam_get_cigar(alignment_ptr), alignment_ptr->core.n_cigar)),
      first_cigar_len(alignment_ptr->core.n_cigar > 0 ? bam_cigar_oplen(bam_get_cigar(alignment_ptr)[0]) : 0),
      last_cigar_len(
          alignment_ptr->core.n_cigar > 0
              ? bam_cigar_oplen(bam_get_cigar(alignment_ptr)[math::SatSub(alignment_ptr->core.n_cigar, u32{1})])
              : 0) {
}

/**
 * @brief Returns the length of the CIGAR operation at the specified index.
 *
 * This method accounts for partial trimming of the first and last CIGAR operations.
 * For untrimmed alignments, this returns the same value as bam_cigar_oplen.
 * For trimmed alignments, the first and last CIGAR operations may have reduced lengths.
 *
 * @param cigar_index The index of the CIGAR operation
 * @return The length of the CIGAR operation at the given index
 */
u32 Alignment::GetCigarLen(const size_t cigar_index) const {
  if (cigar_index == 0) {
    return first_cigar_len;
  }
  if (cigar_index == cigars.size() - 1) {
    return last_cigar_len;
  }
  return bam_cigar_oplen(cigars[cigar_index]);
}

/**
 * @brief Returns the CIGAR operation and its length at the specified index as a Cigar struct.
 *
 * This is a convenience wrapper around GetCigarLen that also extracts the CIGAR operation type.
 *
 * @param cigar_index The index of the CIGAR operation
 * @return A Cigar struct containing the operation type and length
 */
Cigar Alignment::GetCigarAt(const size_t cigar_index) const {
  const u32 cigar = cigars[cigar_index];
  return {bam_cigar_op(cigar), GetCigarLen(cigar_index)};
}

/**
 * @brief Trims the alignment by removing the specified number of bases from the start and end of the read.
 *
 * This method adjusts the alignment by removing bases from the query sequence (read) without modifying
 * the underlying BAM record. It updates the internal representation of the alignment by:
 * - Adjusting the reference position (rpos) to account for leading bases trimmed
 * - Adjusting the query position (qpos) to the first non-trimmed base
 * - Updating the CIGAR span to exclude fully or partially trimmed operations
 * - Recalculating read_length, read_length_without_softclips, and reference_length
 * - Updating first_cigar_len and last_cigar_len for partial CIGAR operation trimming
 *
 * This is primarily used to remove low-quality ends from analysis.
 *
 * The method preserves the original BAM record while allowing the Alignment object to represent
 * a trimmed view of the alignment. This is more efficient than creating a new BAM record and
 * allows multiple metrics to be calculated on different trimmed versions of the same alignment.
 *
 * EDGE CASES:
 * - If trimming more bases than the read length, the alignment is set to zero length
 * - CIGAR operations that don't consume query bases (deletions, ref skips) are handled correctly
 * - Partial trimming of CIGAR operations at boundaries is supported
 * - Trailing deletions/ref skips after trimming are excluded from the final CIGAR span
 * - An alignment can only be trimmed once (throws error on subsequent attempts)
 *
 * @param trim_leading_bases Number of bases to trim from the start of the read
 * @param trim_trailing_bases Number of bases to trim from the end of the read
 */
void Alignment::Trim(const u32 trim_leading_bases, const u32 trim_trailing_bases) {
  if (_trimmed) {
    throw error::Error("Alignment has already been trimmed and cannot be trimmed again.");
  }
  if (trim_leading_bases + trim_trailing_bases == 0) {
    return;
  }
  if (trim_leading_bases + trim_trailing_bases >= read_length) {
    // If we are trimming more bases than the read length, we set the read length to 0 and the reference span to 0
    qpos = 0;
    read_length = 0;
    read_length_without_softclips = 0;
    reference_length = 0;
    cigars = std::span<u32>();
    first_cigar_len = 0;
    last_cigar_len = 0;
    return;
  }
  if (alignment_ptr->core.n_cigar == 0) {
    // The read has no CIGAR operations, we just update the read length and qpos
    qpos = trim_leading_bases;
    read_length = read_length - trim_leading_bases - trim_trailing_bases;
    read_length_without_softclips = read_length;
    return;
  }
  // If we reach here, we have a non-zero number of CIGAR operations and we need to adjust them
  size_t first_cigar_index = 0;
  u32 bases_trimmed_from_start_of_first_cigar = 0;
  u32 bases_to_trim = trim_leading_bases;
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
  bases_to_trim = trim_trailing_bases;
  size_t last_cigar_index = math::SatSub(alignment_ptr->core.n_cigar, u32{1});
  u32 bases_trimmed_from_end_of_last_cigar = 0;
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
  // TODO: I think we also need to handle the edge case where after trimming the first or last CIGAR operation becomes
  // a deletion/ref skip.

  // Edge case: if the last CIGAR after trimming is a deletion/ref skip, we skip past it
  if (last_cigar_index > first_cigar_index) {
    const u32 last_cigar = bam_get_cigar(alignment_ptr)[last_cigar_index];
    const bool consumes_query = bam_cigar_type(bam_cigar_op(last_cigar)) & 1;
    const bool consumes_reference = bam_cigar_type(bam_cigar_op(last_cigar)) & 2;
    if (!consumes_query && consumes_reference) {
      last_cigar_index = math::SatSub(last_cigar_index, size_t{1});
    }
  }
  // Update the CIGAR pointer stored in the Alignment object to point to the trimmed CIGAR operations
  cigars = std::span<u32>(bam_get_cigar(alignment_ptr) + first_cigar_index, last_cigar_index - first_cigar_index + 1);
  if (first_cigar_index == last_cigar_index) {
    first_cigar_len = bam_cigar_oplen(bam_get_cigar(alignment_ptr)[first_cigar_index]) -
                      bases_trimmed_from_start_of_first_cigar - bases_trimmed_from_end_of_last_cigar;
    last_cigar_len = first_cigar_len;
  } else {
    first_cigar_len =
        bam_cigar_oplen(bam_get_cigar(alignment_ptr)[first_cigar_index]) - bases_trimmed_from_start_of_first_cigar;
    last_cigar_len =
        bam_cigar_oplen(bam_get_cigar(alignment_ptr)[last_cigar_index]) - bases_trimmed_from_end_of_last_cigar;
  }
  // Update reference_length
  reference_length = 0;
  for (size_t i = 0; i < cigars.size(); ++i) {
    const u32 cigar = cigars[i];
    const u32 cigar_op = bam_cigar_op(cigar);
    const u32 cigar_len = GetCigarLen(i);
    const bool consumes_reference = bam_cigar_type(cigar_op) & 2;
    if (consumes_reference) {
      reference_length += cigar_len;
    }
  }
  _trimmed = true;
}

u32 Alignment::CountAlignedBases() const {
  u32 count = 0;
  for (size_t cigar_index = 0; cigar_index < cigars.size(); ++cigar_index) {
    const u32 cigar = cigars[cigar_index];
    const u32 cigar_op = bam_cigar_op(cigar);
    if (cigar_op == BAM_CMATCH || cigar_op == BAM_CDIFF || cigar_op == BAM_CEQUAL) {
      count += GetCigarLen(cigar_index);
    }
  }
  return count;
}

}  // namespace xoos::alignment_metrics
