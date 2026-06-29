#include "util/softclip-util.h"

#include <algorithm>
#include <string>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/util/math.h>

#include "io/alignment.h"
#include "util/read-util.h"

namespace xoos::read_collapser {

/**
 * @brief Checks if the alignment has leading soft-clip or insertion CIGAR operations.
 * Insertions immediately following a leading soft-clip are also considered as leading insertions.
 *
 * @param alignment The alignment to check.
 * @return True if the alignment has leading soft-clip or insertion, false otherwise.
 */
static bool HasLeadingSoftclipOrInsertion(const AlignmentPtr& alignment) {
  const u32 n_cigar = alignment->record->core.n_cigar;
  if (n_cigar == 0) {
    return false;
  }
  const bool has_leading_softclip = bam_cigar_op(bam_get_cigar(alignment->record.get())[0]) == BAM_CSOFT_CLIP;
  const bool has_leading_insertion = bam_cigar_op(bam_get_cigar(alignment->record.get())[0]) == BAM_CINS;
  const bool has_leading_insertion_after_softclip =
      alignment->record->core.n_cigar > 1 && bam_cigar_op(bam_get_cigar(alignment->record.get())[1]) == BAM_CINS &&
      has_leading_softclip;
  return has_leading_softclip || has_leading_insertion || has_leading_insertion_after_softclip;
}

/**
 * @brief Checks if the alignment has trailing soft-clip or insertion CIGAR operations.
 * Insertions immediately preceding a trailing soft-clip are also considered as trailing insertions.
 *
 * @param alignment The alignment to check.
 * @return True if the alignment has trailing soft-clip or insertion, false otherwise.
 */
static bool HasTrailingSoftclipOrInsertion(const AlignmentPtr& alignment) {
  const u32 n_cigar = alignment->record->core.n_cigar;
  if (n_cigar == 0) {
    return false;
  }
  const bool has_trailing_softclip =
      bam_cigar_op(bam_get_cigar(alignment->record.get())[n_cigar - 1]) == BAM_CSOFT_CLIP;
  const bool has_trailing_insertion = bam_cigar_op(bam_get_cigar(alignment->record.get())[n_cigar - 1]) == BAM_CINS;
  const bool has_trailing_insertion_before_softclip =
      n_cigar > 1 && bam_cigar_op(bam_get_cigar(alignment->record.get())[n_cigar - 2]) == BAM_CINS &&
      has_trailing_softclip;
  return has_trailing_softclip || has_trailing_insertion || has_trailing_insertion_before_softclip;
}

/**
 * @brief Gets the length of leading soft-clip or insertion CIGAR operations.
 * Insertions immediately following a leading soft-clip are also counted.
 *
 * @param alignment The alignment to check.
 * @return The total length of leading soft-clip and insertion operations; 0 if none exist.
 */
static u32 GetLengthOfLeadingSoftclipOrInsertion(const AlignmentPtr& alignment) {
  u32 length = 0;
  const u32 n_cigar = alignment->record->core.n_cigar;
  const u32* cigar = bam_get_cigar(alignment->record.get());
  if (n_cigar == 0) {
    return length;
  }
  if (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP || bam_cigar_op(cigar[0]) == BAM_CINS) {
    length += bam_cigar_oplen(cigar[0]);
    if (n_cigar > 1 && (bam_cigar_op(cigar[1]) == BAM_CINS) && (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP)) {
      length += bam_cigar_oplen(cigar[1]);
    }
  }
  return length;
}

/**
 * @brief Gets the length of trailing soft-clip or insertion CIGAR operations.
 * Insertions immediately preceding a trailing soft-clip are also counted.
 *
 * @param alignment The alignment to check.
 * @return The total length of trailing soft-clip and insertion operations; 0 if none exist.
 */
static u32 GetLengthOfTrailingSoftclipOrInsertion(const AlignmentPtr& alignment) {
  u32 length = 0;
  const u32 n_cigar = alignment->record->core.n_cigar;
  const u32* cigar = bam_get_cigar(alignment->record.get());
  if (n_cigar == 0) {
    return length;
  }
  if (bam_cigar_op(cigar[n_cigar - 1]) == BAM_CSOFT_CLIP || bam_cigar_op(cigar[n_cigar - 1]) == BAM_CINS) {
    length += bam_cigar_oplen(cigar[n_cigar - 1]);
    if (n_cigar > 1 && (bam_cigar_op(cigar[n_cigar - 2]) == BAM_CINS) &&
        (bam_cigar_op(cigar[n_cigar - 1]) == BAM_CSOFT_CLIP)) {
      length += bam_cigar_oplen(cigar[n_cigar - 2]);
    }
  }
  return length;
}

u32 LeftmostRefPos(const vec<AlignmentPtr>& reads_in_cluster, const bool trim_overhangs) {
  u32 leftmost_rpos = std::numeric_limits<u32>::max();
  if (reads_in_cluster.empty()) {
    return 0;
  }
  for (const auto& alignment : reads_in_cluster) {
    leftmost_rpos = std::min(leftmost_rpos, alignment->StartPos());
  }
  // Adjust the leftmost reference positions if we are trimming overhangs
  // due to soft-clips
  //
  // Example:
  // Without trimming overhangs:
  // Read1 (start=100, end=114):   AAAAACTCCTTTTT
  // Read2 (start=105, end=115):        CTCCTTTTTT
  // Consensus matrix range: [100, 115)
  //
  // With trimming overhangs:
  // Read1 (start=100, end=114):   AAAAACTCCTTTTT
  // Read2 (start=105, end=115):        CTCCTTTTTT
  // Consensus matrix range: [105, 114)
  if (trim_overhangs) {
    for (const auto& alignment : reads_in_cluster) {
      if (alignment->umi5p.has_value() && HasLeadingSoftclipOrInsertion(alignment)) {
        leftmost_rpos = std::max(leftmost_rpos, alignment->StartPos());
      }
    }
  }
  return leftmost_rpos;
}

u32 RightmostRefPos(const vec<AlignmentPtr>& reads_in_cluster, const bool trim_overhangs) {
  u32 rightmost_rpos = std::numeric_limits<u32>::min();
  if (reads_in_cluster.empty()) {
    return 0;
  }
  for (const auto& alignment : reads_in_cluster) {
    rightmost_rpos = std::max(rightmost_rpos, alignment->EndPos());
  }
  // Adjust the rightmost reference positions if we are trimming overhangs
  // due to soft-clips
  //
  // Example:
  // Without trimming overhangs:
  // Read1 (start=100, end=114):   AAAAACTCCTTTTT
  // Read2 (start=105, end=115):        CTCCTTTTTT
  // Consensus matrix range: [100, 115)
  //
  // With trimming overhangs:
  // Read1 (start=100, end=114):   AAAAACTCCTTTTT
  // Read2 (start=105, end=115):        CTCCTTTTTT
  // Consensus matrix range: [105, 114)
  if (trim_overhangs) {
    for (const auto& alignment : reads_in_cluster) {
      if (alignment->umi3p.has_value() && HasTrailingSoftclipOrInsertion(alignment)) {
        rightmost_rpos = std::min(rightmost_rpos, alignment->EndPos());
      }
    }
  }
  return rightmost_rpos;
}

vec<std::string> GetLeftSoftclipSequences(const vec<AlignmentPtr>& reads_in_cluster) {
  vec<std::string> left_softclips;
  if (reads_in_cluster.empty()) {
    return left_softclips;
  }
  // Find the boundary between soft-clips and aligned sequences.
  // This boundary should be consistent across all reads in the cluster that have soft clips and
  // may result in extraction of some aligned bases as part of the soft-clip sequences
  const u32 boundary = LeftmostRefPos(reads_in_cluster, true);
  for (const auto& read : reads_in_cluster) {
    std::string softclip;
    u32 softclip_length_with_overhang = math::SatSub(boundary, read->StartPos());
    u32 read_start = 0;
    if (HasLeadingSoftclipOrInsertion(read)) {
      if (read->umi5p.has_value()) {
        softclip_length_with_overhang += GetLengthOfLeadingSoftclipOrInsertion(read);
      } else {
        read_start = GetLengthOfLeadingSoftclipOrInsertion(read);
      }
    }
    softclip = GetSequence(bam_get_seq(read->record.get()), read_start, softclip_length_with_overhang);
    left_softclips.push_back(softclip);
  }

  return left_softclips;
}

vec<std::string> GetRightSoftclipSequences(const vec<AlignmentPtr>& reads_in_cluster) {
  vec<std::string> right_softclips;
  if (reads_in_cluster.empty()) {
    return right_softclips;
  }
  // Find the boundary between aligned sequences and soft-clips.
  // This boundary should be consistent across all reads in the cluster that have soft clips and
  // may result in extraction of some aligned bases as part of the soft-clip sequences
  const u32 boundary = RightmostRefPos(reads_in_cluster, true);
  for (const auto& read : reads_in_cluster) {
    std::string softclip;
    u32 softclip_length_with_overhang = math::SatSub(read->EndPos(), boundary);
    u32 right_offset = math::SatSub(read->EndPos(), boundary);
    const u32 read_length = ToUnsigned(read->record->core.l_qseq);
    if (HasTrailingSoftclipOrInsertion(read)) {
      if (read->umi3p.has_value()) {
        softclip_length_with_overhang += GetLengthOfTrailingSoftclipOrInsertion(read);
      }
      right_offset += GetLengthOfTrailingSoftclipOrInsertion(read);
    }
    softclip = GetSequence(bam_get_seq(read->record.get()), read_length - right_offset, softclip_length_with_overhang);
    right_softclips.push_back(softclip);
  }

  return right_softclips;
}

}  // namespace xoos::read_collapser
