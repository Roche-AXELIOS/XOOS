#include "core/alignment.h"

#include <algorithm>
#include <optional>
#include <span>
#include <string>

#include <htslib/sam.h>

#include <xoos/error/error.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/io/htslib-util/read-util.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/util/math.h>

#include "alignment-metrics-options.h"
#include "core/interval.h"
#include "core/region-lookup.h"
#include "core/supplementary-alignment.h"
#include "metadata/alignment-metadata.h"
#include "util/alignment-trimming-util.h"
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
  return {.type = bam_cigar_op(cigar), .length = GetCigarLen(cigar_index)};
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
  // If we reach here, we have a non-zero number of CIGAR operations and we need to adjust them.
  // Apply the leading trim's results before computing the trailing trim, since the trailing trim
  // walks from the read_length/read_length_without_softclips left by the leading trim.
  const TrimmingResult leading_result = TrimLeading(*this, trim_leading_bases);
  rpos = leading_result.new_rpos;
  qpos = leading_result.new_qpos;
  read_length = leading_result.new_read_length;
  read_length_without_softclips = leading_result.new_read_length_without_softclips;
  const size_t first_cigar_index = leading_result.boundary_cigar_index;
  const u32 bases_trimmed_from_start_of_first_cigar = leading_result.bases_trimmed_from_boundary_cigar;
  const TrimmingResult trailing_result = TrimTrailing(*this, trim_trailing_bases);
  // Trimming trailing bases does not change rpos/qpos.
  read_length = trailing_result.new_read_length;
  read_length_without_softclips = trailing_result.new_read_length_without_softclips;
  size_t last_cigar_index = trailing_result.boundary_cigar_index;
  const u32 bases_trimmed_from_end_of_last_cigar = trailing_result.bases_trimmed_from_boundary_cigar;
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
  // Update the CIGAR pointer stored in the Alignment object to point to the trimmed CIGAR operations.
  // bam_get_cigar returns a raw C array from htslib, so pointer arithmetic is required to build the span.
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
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
  // Recompute reference_length by summing the reference-consuming lengths of the (possibly trimmed) CIGAR operations.
  reference_length = 0;
  for (size_t i = 0; i < cigars.size(); ++i) {
    const u32 cigar = cigars[i];
    const bool consumes_reference = bam_cigar_type(bam_cigar_op(cigar)) & 2;
    if (consumes_reference) {
      reference_length += GetCigarLen(i);
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

namespace {

// Per-base coverage flags for a single soft-clipped base of the primary read, recording whether it
// is mapped/aligned by some supplementary alignment.
struct SoftClipBaseMappingStatus {
  bool mapped_by_supplementary{};
  bool aligned_by_supplementary{};
  bool aligned_by_supplementary_passing_filter{};
};

// Tracks, for one contiguous soft-clipped region of the trimmed primary read, whether each base is
// mapped (match/mismatch/insertion) or aligned (match/mismatch) by some supplementary alignment.
//
// Positions are expressed in the primary read's stored-query coordinates (the orientation in which
// the read is stored in the BAM, i.e. the primary alignment's orientation). Only soft-clipped
// regions are tracked, so the flag array stays small even when the read is mostly mapped.
struct SoftClipRegion {
  // Stored-query coordinate of the first base in the region.
  s64 start{};
  // Per-base flags, indexed by (stored position - start).
  vec<SoftClipBaseMappingStatus> bases;

  explicit SoftClipRegion(const s64 region_start, const s64 length)
      : start{region_start}, bases(static_cast<size_t>(length)) {
  }

  s64 End() const {
    return start + static_cast<s64>(bases.size());
  }
};

// Walk one supplementary alignment's CIGAR (already oriented to the primary's stored-query frame),
// marking which bases in the primary's soft-clipped regions it maps or aligns.
void RecordSoftClipRegionCoverageBySA(const vec<u32>& supplementary_cigar,
                                      const bool passes_filter,
                                      vec<SoftClipRegion>& soft_clip_regions) {
  s64 supplementary_read_position = 0;
  for (const u32 cigar : supplementary_cigar) {
    const u32 cigar_op = bam_cigar_op(cigar);
    const u32 cigar_len = bam_cigar_oplen(cigar);
    const bool consumes_query = (bam_cigar_type(cigar_op) & u32{1}) != 0;
    // Hard clips do not consume the query but still advance the read cursor in the SA-tag frame.
    if (!consumes_query && cigar_op != BAM_CHARD_CLIP) {
      continue;
    }
    const bool is_aligned = cigar_op == BAM_CMATCH || cigar_op == BAM_CDIFF || cigar_op == BAM_CEQUAL;
    const bool is_mapped = is_aligned || cigar_op == BAM_CINS;
    if (is_mapped) {
      const s64 op_start = supplementary_read_position;
      const s64 op_end = supplementary_read_position + cigar_len;
      // Record coverage over the (small) soft-clipped regions only.
      for (SoftClipRegion& region : soft_clip_regions) {
        const s64 overlap_start = std::max(op_start, region.start);
        const s64 overlap_end = std::min(op_end, region.End());
        for (s64 pos = overlap_start; pos < overlap_end; ++pos) {
          SoftClipBaseMappingStatus& flags = region.bases[static_cast<size_t>(pos - region.start)];
          flags.mapped_by_supplementary = true;
          flags.aligned_by_supplementary = flags.aligned_by_supplementary || is_aligned;
          flags.aligned_by_supplementary_passing_filter =
              flags.aligned_by_supplementary_passing_filter || (is_aligned && passes_filter);
        }
      }
    }
    supplementary_read_position += cigar_len;
  }
}

// Marks the primary's soft-clipped bases mapped/aligned by a single supplementary alignment and adds
// its aligned bases to `counts`. Supplementary alignments whose reference span does not overlap a
// region of interest (when a BED file was supplied) are skipped. `primary_flag` and
// `primary_is_reverse` describe the primary alignment; `soft_clip_regions` is updated in place.
void AccumulateSupplementaryCoverage(const SaTag& supplementary,
                                     const SaTagFilter& filter,
                                     const RegionLookupTable& regions_of_interest,
                                     const u32 primary_flag,
                                     const bool primary_is_reverse,
                                     vec<SoftClipRegion>& soft_clip_regions,
                                     SoftClipBaseCounts& counts) {
  // A supplementary alignment only contributes when its own reference span overlaps a region of
  // interest. The reference span is orientation-independent, so the SA CIGAR is used as parsed (no
  // reversal). SA positions are 0-based, half-open on the end. When no BED file was supplied, every
  // supplementary alignment contributes and the overlap check is skipped entirely.
  const s64 supplementary_reference_length =
      bam_cigar2rlen(static_cast<int>(supplementary.cigar.size()), supplementary.cigar.data());
  const Interval supplementary_span{supplementary.pos, supplementary.pos + supplementary_reference_length};
  if (regions_of_interest.HasBedFile() &&
      !regions_of_interest.OverlapsAnyRegionOfInterest(supplementary.rname, supplementary_span)) {
    return;
  }

  // Supplementary alignment flag inferred from the primary alignment's flags and the SA tag's strand field.
  // Inheriting the primary's BAM_FDUP is important: many tools do not duplicate-mark supplementary
  // records whose primary is a duplicate, so the primary's duplicate status must drive the
  // passing-filter counts here rather than the supplementary record's own (possibly missing) flag.
  const u32 inferred_supplementary_flag = u32{BAM_FSUPPLEMENTARY} | (primary_flag & u32{BAM_FDUP}) |
                                          (supplementary.is_reverse ? u32{BAM_FREVERSE} : u32{0});
  const bool passes_filter =
      (supplementary.mapq >= filter.min_mapq) && ((inferred_supplementary_flag & filter.exclude_flags) == 0);
  const vec<u32> supplementary_cigar = (supplementary.is_reverse != primary_is_reverse)
                                           ? vec<u32>(supplementary.cigar.rbegin(), supplementary.cigar.rend())
                                           : supplementary.cigar;
  // Record which of the primary's soft-clipped bases are mapped or aligned by this supplementary.
  RecordSoftClipRegionCoverageBySA(supplementary_cigar, passes_filter, soft_clip_regions);

  // Tally this supplementary's aligned bases toward aligned_bases. The aligned-base count is
  // orientation-independent, so the SA CIGAR is used as parsed (no reversal).
  u32 supplementary_aligned_bases = 0;
  for (const u32 op : supplementary.cigar) {
    const u32 cigar_op = bam_cigar_op(op);
    if (cigar_op == BAM_CMATCH || cigar_op == BAM_CDIFF || cigar_op == BAM_CEQUAL) {
      supplementary_aligned_bases += bam_cigar_oplen(op);
    }
  }
  counts.supplementary_aligned_bases += supplementary_aligned_bases;
  if (passes_filter) {
    counts.supplementary_aligned_bases_passing_filter += supplementary_aligned_bases;
  }
}

}  // namespace

SoftClipBaseCounts Alignment::CountSoftClippedBasesWithSupplementary(
    const SaTagFilter& filter, const RegionLookupTable& regions_of_interest) const {
  SoftClipBaseCounts counts{};
  if (read_length == 0) {
    return counts;
  }
  const auto sa_tag = io::BamAuxGet<std::string>(alignment_ptr, "SA");
  if (!sa_tag.has_value()) {
    // No SA tag, so all soft-clipped bases are unaligned and unmapped
    return counts;
  }

  // Identify the soft-clipped regions of the trimmed primary alignment, if any. Soft clips can only
  // appear at the ends of the CIGAR, so there are at most two regions. Tracking is allocated only
  // over these (typically small) regions rather than the whole, mostly-mapped read. Positions are
  // expressed in the primary read's query coordinates: the trimmed read occupies
  // [qpos, qpos + read_length), so the leading soft clip starts at qpos and the trailing soft clip
  // ends at qpos + read_length.
  vec<SoftClipRegion> soft_clip_regions;
  if (!cigars.empty() && bam_cigar_op(cigars.front()) == BAM_CSOFT_CLIP) {
    soft_clip_regions.emplace_back(qpos, GetCigarLen(0));
  }
  if (cigars.size() > 1 && bam_cigar_op(cigars.back()) == BAM_CSOFT_CLIP) {
    const u32 trailing_soft_clip_len = GetCigarLen(cigars.size() - 1);
    soft_clip_regions.emplace_back(qpos + read_length - trailing_soft_clip_len, trailing_soft_clip_len);
  }

  const vec<SaTag> supplementary_alignments = ParseSaTag(sa_tag.value());
  const u32 primary_flag = alignment_ptr->core.flag;
  const bool primary_is_reverse = (primary_flag & u32{BAM_FREVERSE}) != 0;
  // Walk each supplementary alignment's SA-tag CIGAR, marking which bases in the primary's
  // soft-clipped regions are mapped or aligned by that supplementary alignment. SA-tag CIGARs are
  // expressed in the supplementary's own orientation, so a CIGAR is reversed when the supplementary
  // is on the opposite strand of the primary, aligning its read cursor with the primary's
  // stored-query coordinate frame.
  for (const SaTag& supplementary : supplementary_alignments) {
    AccumulateSupplementaryCoverage(
        supplementary, filter, regions_of_interest, primary_flag, primary_is_reverse, soft_clip_regions, counts);
  }

  for (const SoftClipRegion& region : soft_clip_regions) {
    for (const SoftClipBaseMappingStatus& flags : region.bases) {
      if (flags.mapped_by_supplementary) {
        ++counts.mapped_soft_clipped_bases;
      }
      if (flags.aligned_by_supplementary) {
        ++counts.aligned_soft_clipped_bases;
      }
      if (flags.aligned_by_supplementary_passing_filter) {
        ++counts.aligned_soft_clipped_bases_passing_filter;
      }
    }
  }
  return counts;
}

}  // namespace xoos::alignment_metrics
