#include "util/duplex-util.h"

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <xoos/error/error.h>
#include <xoos/types/int.h>
#include <xoos/util/math.h>

namespace xoos::read_collapser {

std::optional<f64> GetDiscordantDuplexErrorRate(const bam1_t* record) {
  yc_decode::YcTag yc_tag = yc_decode::DeserializeYcTag(record);
  // Check if the YC tag indicates a duplex read.
  if (!yc_tag.IsDuplex()) {
    return std::nullopt;
  }

  // NOTE: The denominator cannot be zero because we
  // know that the duplex segment is not empty if we reach this point.
  return static_cast<f64>(yc_tag.CountDiscordantBase()) / static_cast<f64>(yc_tag.CountDuplexBase());
}

/**
 * Update the base qualities of an unmapped, duplex read by replacing any concordant base qualities with the deconvolved
 * quality score.
 *
 * @param read The BAM record.
 */
void UpdateUnmappedDeconvolvedQualities(const bam1_t* read) {
  const auto yc_tag = yc_decode::DeserializeYcTag(read);
  const auto base_types = yc_tag.GetBaseTypes();
  const auto quals = bam_get_qual(read);
  // Ensure the base types length matches the quality scores and read sequence length
  yc_decode::RequireConsistentReadLength(read, base_types);
  for (size_t i = 0; i < base_types.size(); ++i) {
    if (base_types[i] == yc_decode::BaseType::kConcordant) {
      quals[i] = kConcordantDeconvolvedBaseQual;
    }
  }
}

/**
 * Get the average cluster size for a deconvolved read based on its base types.
 *
 * @param read The BAM record.
 *
 * @return The average cluster size for the deconvolved read.
 */
u32 GetUnmappedDeconvolvedAvgClusterSize(const bam1_t* read) {
  using enum yc_decode::BaseType;
  const auto yc_tag = yc_decode::DeserializeYcTag(read);
  const auto base_types = yc_tag.GetBaseTypes();
  // Ensure the base types length matches the quality scores and read sequence length
  yc_decode::RequireConsistentReadLength(read, base_types);
  u32 concordant_base_count = 0;
  u32 discordant_base_count = 0;
  u32 simplex_base_count = 0;
  for (const auto& base_type : base_types) {
    switch (base_type) {
      case kDiscordant:
        ++discordant_base_count;
        break;
      case kSimplex:
        ++simplex_base_count;
        break;
      case kConcordant:
        ++concordant_base_count;
        break;
      default:
        throw yc_decode::BaseTypeError(bam_get_qname(read), base_type);
    }
  }
  return 1 + math::DivideAndRound(concordant_base_count,
                                  concordant_base_count + simplex_base_count + discordant_base_count);
}

vec<AlignmentPtr> DeconvolveDuplexReads(const vec<AlignmentPtr>& original_alignments, const bool is_parent_parent) {
  vec<AlignmentPtr> deconvolved_alignments;
  for (const auto& original_alignment : original_alignments) {
    const auto yc_tag = yc_decode::DeserializeYcTag(original_alignment->record.get());
    auto deconvolved_bam_records = yc_decode::DecodeToBamRecords(original_alignment->record.get());
    if (std::holds_alternative<std::pair<io::Bam1Ptr, io::Bam1Ptr>>(deconvolved_bam_records)) {
      // Two records deconvolved, add both to the alignments
      auto [r1, r2] = std::move(std::get<std::pair<io::Bam1Ptr, io::Bam1Ptr>>(deconvolved_bam_records));
      auto alignment_r1 = std::make_shared<Alignment>(std::move(r1));
      auto alignment_r2 = std::make_shared<Alignment>(std::move(r2));

      alignment_r1->cluster = original_alignment->cluster;
      // NOTE: `mate_has_left_overhang` and `mate_has_right_overhang` are set for
      // the read not containing the overhangs (in other words, the mate of the read with the overhangs).
      alignment_r1->mate_has_left_overhang = !yc_tag.left_overhang_is_r1 && yc_tag.left_overhang > 0;
      alignment_r1->mate_has_right_overhang = !yc_tag.right_overhang_is_r1 && yc_tag.right_overhang > 0;
      alignment_r1->duplex_strand = DuplexStrand::kR1;
      if (is_parent_parent) {
        // Set R1 to the forward strand for parent-parent workflow by clearing the reverse flag if it was set
        alignment_r1->record->core.flag &= ~BAM_FREVERSE;
      }
      alignment_r2->cluster = original_alignment->cluster;
      alignment_r2->mate_has_left_overhang = yc_tag.left_overhang_is_r1 && yc_tag.left_overhang > 0;
      alignment_r2->mate_has_right_overhang = yc_tag.right_overhang_is_r1 && yc_tag.right_overhang > 0;
      alignment_r2->duplex_strand = DuplexStrand::kR2;
      if (is_parent_parent) {
        // Set R2 to the reverse strand for parent-parent workflow
        alignment_r2->record->core.flag |= BAM_FREVERSE;
      }
      deconvolved_alignments.emplace_back(std::move(alignment_r1));
      deconvolved_alignments.emplace_back(std::move(alignment_r2));
    } else if (std::holds_alternative<io::Bam1Ptr>(deconvolved_bam_records)) {
      // One record deconvolved, use it as is
      auto alignment_r1 = std::make_shared<Alignment>(std::move(std::get<io::Bam1Ptr>(deconvolved_bam_records)));
      alignment_r1->cluster = original_alignment->cluster;
      deconvolved_alignments.emplace_back(std::move(alignment_r1));
    } else {
      // Simplex read, no deconvolution needed
      deconvolved_alignments.emplace_back(original_alignment);
    }
  }
  return deconvolved_alignments;
}

}  // namespace xoos::read_collapser
