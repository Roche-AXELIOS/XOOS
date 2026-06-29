#include "consensus/consensus-matrix-builder.h"

#include <ranges>
#include <unordered_map>

#include <htslib/sam.h>

#include <xoos/error/error.h>
#include <xoos/types/int.h>
#include <xoos/util/math.h>
#include <xoos/yc-decode/yc-decoder.h>

#include "consensus/base-encoder.h"
#include "consensus/consensus.h"
#include "consensus/poa-msa.h"
#include "io/alignment.h"
#include "util/read-util.h"
#include "util/softclip-util.h"

namespace xoos::read_collapser {

struct Positions {
  size_t index{};
  s64 qpos{};
  s64 rpos{};
  // This is the position in the new gapped coordinate
  size_t new_pos{};
  // This is the position in the old coordinate
  size_t old_pos{};
  char previous_base{};
  u64 hp_start{};
  u64 hp_end{};

  u32 leftmost_rpos{};
  u32 rightmost_rpos{};

  Positions(const size_t read_index,
            const size_t read_start_pos,
            const InsertionMap& aligned_insertions,
            const u32 leftmost_pos,
            const u32 rightmost_pos)
      : index(read_index),
        qpos(0),
        rpos(static_cast<s64>(read_start_pos)),
        new_pos(static_cast<size_t>((rpos >= leftmost_pos) ? rpos - leftmost_pos : 0)),
        old_pos(0),
        previous_base(kBaseP),
        hp_start(0),
        hp_end(0),
        leftmost_rpos(leftmost_pos),
        rightmost_rpos(rightmost_pos) {
    // Adjust the start position if there is any read that starts
    // before the current read and has an insertion before the current position.
    for (s64 pos = leftmost_rpos; pos < rpos; ++pos) {
      // Check if an insertion exists at this position and if this position
      // if within the valid range of the consensus matrix
      if (aligned_insertions.contains(pos) && pos < rightmost_rpos) {
        new_pos += aligned_insertions.at(pos).front().size();
      }
    }
  }
};

// Calculate the maximum width of the consensus sequence after accounting for insertions
// in all alignments in the cluster. This will be the width of the consensus matrix.
static inline size_t ConsensusMatrixWidth(const InsertionMap& aligned_insertions,
                                          const u32 leftmost_rpos,
                                          const u32 rightmost_rpos) {
  size_t total_aligned_insertion_length = 0;
  for (const auto& [position, aligned_insertion] : aligned_insertions) {
    // Only consider insertions within the leftmost and rightmost reference positions.
    // The leftmost and rightmost positions are usually just the smallest start and largest end
    // positions across all reads. However, if overhang trimming is enabled, these positions are
    // adjusted so we need to check that the insertion falls within this range.
    if (position >= leftmost_rpos && position < rightmost_rpos) {
      // After MSA, all sequences have the same length for the same position
      // so we only need to look at the first one.
      total_aligned_insertion_length += aligned_insertion.front().length();
    }
  }
  return total_aligned_insertion_length + (rightmost_rpos - leftmost_rpos);
}

/**
 * @brief Preprocesses the consensus matrix for a specific read.
 *
 * This method prepares the consensus matrix for the given `read_index` by analyzing
 * the alignment data provided in the `read` parameter. It ensures that the matrix
 * is updated with the necessary information to reflect the alignment of the read
 * to the consensus sequence.
 *
 * @param read_index The index of the read in the cluster to preprocess.
 * @param read A pointer to the alignment data for the read.
 */
static inline void PreprocessConsensusMatrix(ConsensusMatrix& matrix,
                                             const size_t read_index,
                                             const AlignmentPtr& read) {
  const bool is_reverse = read->record->core.flag & BAM_FREVERSE;
  const bool is_partial = read->IsPartial();
  const bool has_duplex_overhang = read->mate_has_left_overhang || read->mate_has_right_overhang;

  matrix.SetStrand(read_index, is_reverse ? ReadStrand::kRev : ReadStrand::kFwd);
  // If the read is a duplex read, this will assign R1/R2 to the read
  // If the read is simplex, it will set the duplex strand to kSimplex, indicating that it is a simplex read.
  matrix.SetDuplexStrand(read_index, read->duplex_strand);

  // Update metadata
  matrix.UpdateMetadata(is_reverse, is_partial);

  // Fill the matrix with 'P' for partial reads or reads with duplex overhangs
  if (is_partial || has_duplex_overhang) {
    // For partial read, leading and trailing gaps should be set to 'P'; 'P' is used to indicate the lack of
    // coverage from the current read at a given position.
    // For read with duplex overhangs, we have a left or right overhang in the other strand of the duplex read,
    // so we need to fill the leading and trailing gaps with 'P' similar to a partial read
    // Example:
    // Duplex read AAAAACTCCC with YC tag 5-5-0 -> deconvolves into AAAAACTCCC and CTCCC
    // R1 has a left overhang of 5 bases (AAAAA) so we need to pad R2 with 'P' for the leading gaps
    // to indicate that the first 5 bases are only covered by R1, giving us the following consensus matrix:
    // R1: AAAAACTCCC
    // R2: PPPPPCTCCC
    matrix.FillRow(read_index, kBaseP);
  }
}

static inline void UpdateHomopolymerRange(ConsensusMatrix& matrix, Positions& positions, const char current_base) {
  if (current_base != kBaseGap) {
    if (current_base == positions.previous_base) {
      // Extend the homopolymer range
      positions.hp_end = positions.new_pos;
    } else {
      // We encountered a different base, so we close the previous homopolymer range,
      // record it, and start a new one.
      if (positions.hp_start < positions.hp_end) {
        matrix.AddHomopolymerRange(positions.hp_start, positions.hp_end);
      }
      // Start a new homopolymer range
      positions.hp_start = positions.new_pos;
      positions.hp_end = positions.new_pos;
      positions.previous_base = current_base;
    }
  }
}

/**
 * @brief Adds a base to the consensus matrix at the specified position.
 * This function also updates the homopolymer range if HD deconvolution is enabled.
 *
 * @param matrix The consensus matrix to update.
 * @param base_types The vector of base types corresponding to the read. Can be empty if no base type information is
 * available.
 * @param positions The current positions in the read and consensus matrix.
 * @param current_base The base to add to the consensus matrix.
 * @param hd_deconvolve_enabled A boolean indicating whether HD deconvolution is enabled.
 */
static inline void AddBase(ConsensusMatrix& matrix,
                           const vec<yc_decode::BaseType>& base_types,
                           Positions& positions,
                           const char current_base,
                           const bool hd_deconvolve_enabled) {
  if (hd_deconvolve_enabled) {
    // Update the homopolymer range based on the new character if HD deconvolution is enabled
    UpdateHomopolymerRange(matrix, positions, current_base);
  }

  // Set the base in the consensus matrix
  matrix.SetBase(positions.index, positions.new_pos, current_base);

  // Early exit if there is no base type information
  if (!matrix.HasBaseTypeInfo() || base_types.empty()) {
    ++positions.new_pos;
    return;
  }

  // Update base type if it's not a gap
  if (current_base != kBaseGap && positions.old_pos < base_types.size()) {
    matrix.SetBaseType(positions.index, positions.new_pos, base_types[positions.old_pos]);
    ++positions.old_pos;
  }

  ++positions.new_pos;
}

/**
 * @brief Handles the processing of match CIGAR operations in the consensus matrix.
 * This function updates the consensus matrix with bases from the alignment,
 * taking into account any insertions at the current reference position.
 *
 * @param matrix The consensus matrix to update.
 * @param base_types The vector of base types corresponding to the read. Can be empty if no base type information is
 * available.
 * @param alignment The alignment containing the read data.
 * @param positions The current positions in the read and consensus matrix.
 * @param cigar_len The length of the match CIGAR operation.
 * @param aligned_insertions A map of aligned insertions for the reads in the cluster.
 * @param hd_deconvolve_enabled A boolean indicating whether HD deconvolution is enabled.
 */
static inline void HandleMatchCigar(ConsensusMatrix& matrix,
                                    const AlignmentPtr& alignment,
                                    const vec<yc_decode::BaseType>& base_types,
                                    Positions& positions,
                                    const size_t cigar_len,
                                    const InsertionMap& aligned_insertions,
                                    const bool hd_deconvolve_enabled) {
  // For matches, we need to copy the base from the read
  // and also check for insertions at the current reference position.
  // If at least one read has an insertion at the current reference position,
  // we need to add gaps/insertions to the consensus matrix prior to the match.
  for (size_t j = 0; j < cigar_len; ++j) {
    if (positions.rpos >= positions.leftmost_rpos && positions.rpos < positions.rightmost_rpos) {
      // Handle aligned insertions at the current position
      if (aligned_insertions.contains(positions.rpos)) {
        const bool is_at_read_start = (positions.rpos == alignment->StartPos()) && (j == 0);
        const auto& insertion_bases = aligned_insertions.at(positions.rpos)[positions.index];
        for (const char insertion_base : insertion_bases) {
          if (alignment->IsPartial() && is_at_read_start) {
            // For partial reads, the matrix was already filled with 'P'.
            // We just need to advance new_pos for the length of the gapped insertion
            // to leave the 'P' bases that are already there.
            ++positions.new_pos;
          } else {
            AddBase(matrix, base_types, positions, insertion_base, hd_deconvolve_enabled);
          }
        }
      }

      // Handle the matched or mismatched base
      const char current_base = GetBase(bam_get_seq(alignment->record.get()), static_cast<size_t>(positions.qpos));
      AddBase(matrix, base_types, positions, current_base, hd_deconvolve_enabled);
    }

    ++positions.rpos;
    ++positions.qpos;
  }
}

/**
 * @brief Handles the processing of deletion CIGAR operations in the consensus matrix.
 * This function updates the consensus matrix with gaps for deletions,
 * taking into account any insertions at the current reference position.
 *
 * @param matrix The consensus matrix to update.
 * @param alignment The alignment containing the read data.
 * @param base_types The vector of base types corresponding to the read. Can be empty if no base type information is
 * available.
 * @param positions The current positions in the read and consensus matrix.
 * @param cigar_len The length of the deletion CIGAR operation.
 * @param aligned_insertions A map of aligned insertions for the reads in the cluster.
 * @param hd_deconvolve_enabled A boolean indicating whether HD deconvolution is enabled.
 */
static inline void HandleDeletionCigar(ConsensusMatrix& matrix,
                                       const AlignmentPtr& alignment,
                                       const vec<yc_decode::BaseType>& base_types,
                                       Positions& positions,
                                       const size_t cigar_len,
                                       const InsertionMap& aligned_insertions,
                                       const bool hd_deconvolve_enabled) {
  // Deletions break the homopolymer range, so we need to add any unclosed range
  if (positions.hp_start < positions.hp_end && hd_deconvolve_enabled) {
    matrix.AddHomopolymerRange(positions.hp_start, positions.hp_end);
  }
  // Start a new homopolymer range
  positions.hp_start = positions.new_pos;
  positions.hp_end = positions.new_pos;
  positions.previous_base = kBaseGap;

  // For deletions, we skip the reference position but still check for insertions
  // at the current reference position.
  for (size_t j = 0; j < cigar_len; ++j) {
    if (positions.rpos >= positions.leftmost_rpos && positions.rpos < positions.rightmost_rpos) {
      if (aligned_insertions.contains(positions.rpos)) {
        const auto& insertion_bases = aligned_insertions.at(positions.rpos)[positions.index];
        for (const char insertion_base : insertion_bases) {
          AddBase(matrix, base_types, positions, insertion_base, hd_deconvolve_enabled);
        }
      }
      // Always add a gap for the deletion itself
      AddBase(matrix, base_types, positions, kBaseGap, hd_deconvolve_enabled);
    }
    // Move to the next reference position
    ++positions.rpos;
  }
}

/**
 * @brief Fills the consensus matrix with bases from the given alignment.
 *
 * This function processes the CIGAR string of the provided alignment and updates
 * the consensus matrix accordingly. It handles matches, deletions, insertions,
 * and soft clips, while also recording homopolymer ranges if HD deconvolution is enabled.
 *
 * @param matrix The consensus matrix to be filled.
 * @param alignment A shared pointer to the alignment containing the read data.
 * @param positions The current positions in the read and consensus matrix.
 * @param aligned_insertions A map of aligned insertions for the reads in the cluster.
 * @param hd_deconvolve_enabled A boolean indicating whether HD deconvolution is enabled.
 */
static inline void FillConsensusMatrix(ConsensusMatrix& matrix,
                                       const AlignmentPtr& alignment,
                                       Positions& positions,
                                       const InsertionMap& aligned_insertions,
                                       const bool hd_deconvolve_enabled) {
  // Keep track of homopolymer ranges
  positions.previous_base = kBaseGap;
  positions.hp_start = positions.new_pos;
  positions.hp_end = positions.new_pos;

  vec<yc_decode::BaseType> base_types;
  if (matrix.HasBaseTypeInfo()) {
    base_types = alignment->GetBaseTypes();
  }

  // Iterate through the CIGAR operations to fill the consensus matrix
  const u32* cigar = bam_get_cigar(alignment->record.get());
  for (u32 cigar_index = 0; cigar_index < alignment->record->core.n_cigar; ++cigar_index) {
    const u32 cigar_op = bam_cigar_op(cigar[cigar_index]);
    const u32 cigar_len = bam_cigar_oplen(cigar[cigar_index]);
    switch (cigar_op) {
      case BAM_CMATCH:
      case BAM_CDIFF:
      case BAM_CEQUAL:  // Match, Diff, Equal
        HandleMatchCigar(
            matrix, alignment, base_types, positions, cigar_len, aligned_insertions, hd_deconvolve_enabled);
        break;
      case BAM_CDEL:
        HandleDeletionCigar(
            matrix, alignment, base_types, positions, cigar_len, aligned_insertions, hd_deconvolve_enabled);
        break;
      case BAM_CINS:        // Insertions
      case BAM_CSOFT_CLIP:  // Soft Clip
        positions.qpos += cigar_len;
        break;
      case BAM_CREF_SKIP:  // Reference Skip
        positions.rpos += cigar_len;
        break;
      default:
        break;
    }
  }
  // If we have an unclosed homopolymer range at the end of the read,
  // we need to add it to the consensus matrix after the main loop.
  // Note that we only add the homopolymer range if it is not empty.
  if (hd_deconvolve_enabled && positions.hp_start < positions.hp_end) {
    matrix.AddHomopolymerRange(positions.hp_start, positions.hp_end);
  }
}

/**
 * @brief Collects and processes insertion events for a given read alignment.
 *
 * This function iterates through the CIGAR string of the provided alignment to identify
 * insertion operations (CIGAR operation 'I'). It ensures that only valid insertions
 * (not at the beginning or end of a read) are processed. Malformed insertions, such as
 * those adjacent to soft-clips, are ignored.
 *
 * @param alignment A pointer to the alignment object containing the read and its associated data.
 * @param read_index The index of the current read being processed.
 * @param insertions A map to store processed insertion data for all reads.
 * @param total_reads The total number of reads being processed.
 */
static void CollectInsertionsForRead(const AlignmentPtr& alignment,
                                     const size_t read_index,
                                     InsertionMap& insertions,
                                     const size_t total_reads) {
  s64 rpos = alignment->StartPos();
  s64 qpos = 0;
  const u32* cigar = bam_get_cigar(alignment->record.get());
  for (u32 cigar_index = 0; cigar_index < alignment->record->core.n_cigar; ++cigar_index) {
    const u32 cigar_op = bam_cigar_op(cigar[cigar_index]);
    const u32 cigar_len = bam_cigar_oplen(cigar[cigar_index]);
    if (cigar_op == BAM_CINS) {
      // Only process the insertions that are not at the beginning or the end of a read
      // Example: 1I2M, 2S1I2M, 2M1I2S, 2M1I are malformed and will be treated as soft-clips
      const bool is_leading_insertion =
          cigar_index == 0 || (cigar_index > 0 && bam_cigar_op(cigar[cigar_index - 1]) == BAM_CSOFT_CLIP);
      const bool is_trailing_insertion =
          cigar_index == alignment->record->core.n_cigar - 1 ||
          (cigar_index < alignment->record->core.n_cigar - 1 && bam_cigar_op(cigar[cigar_index + 1]) == BAM_CSOFT_CLIP);
      if (cigar_index > 0 && !is_leading_insertion && !is_trailing_insertion) {
        if (!insertions.contains(rpos)) {
          insertions[rpos] = std::vector<std::string>(total_reads);
        }
        insertions[rpos][read_index] =
            GetSequence(bam_get_seq(alignment->record.get()), static_cast<u64>(qpos), cigar_len);
      }
    }
    // CIGAR consumes query
    if (bam_cigar_type(cigar_op) & 1) {
      qpos += cigar_len;
    }
    // CIGAR consumes reference
    if (bam_cigar_type(cigar_op) & 2) {
      rpos += cigar_len;
    }
  }
}

/**
 * @brief Aligns insertion sequences using multiple sequence alignment (MSA).
 *
 * This function takes a map of insertion positions and their corresponding sequences,
 * and aligns the sequences at each position using a PoaAligner. The aligned sequences
 * are then updated in the input map.
 *
 * @param insertions A map where the key is the position of the insertion (int) and the value
 *                   is a vector of strings representing the sequences to be aligned.
 */
static void AlignInsertions(InsertionMap& insertions) {
  PoaAligner msa_aligner(kMatchScore, kMismatchPenalty, kGapOpenPenalty, kGapExtPenalty);
  for (const auto& [position, insertionSequences] : insertions) {
    insertions[position] = msa_aligner.ComputeMSA(insertionSequences);
  }
}

/**
 * Calculate the maximum possible length of the consensus sequence.
 * For us to call consensus, there must be at least some overlaps between the reads in the cluster.
 * Therefore, we can calculate the maximum consensus length
 * as the sum of the lengths of all gapped aligned insertions
 * plus the sum of the reference span covered by the reads in the cluster.
 * This is equivalent to the consensus matrix width assuming all reads in the cluster
 * have no overlaps.
 */
static size_t MaxConsensusLength(const vec<AlignmentPtr>& reads_in_cluster, const InsertionMap& aligned_insertions) {
  size_t max_length = 0;
  for (const auto& read : reads_in_cluster) {
    if (read->EndPos() < read->StartPos()) {
      throw error::Error("Read end position is smaller than start position for read {}. Invalid alignment.",
                         bam_get_qname(read->record.get()));
    }
    max_length += read->EndPos() - read->StartPos();
  }
  for (const auto& [position, aligned_insertion] : aligned_insertions) {
    // After MSA, all sequences have the same length for the same position
    // so we only need to look at the first one.
    max_length += aligned_insertion.front().length();
  }
  return max_length;
}

InsertionMap GetAlignedInsertions(const vec<AlignmentPtr>& reads_in_cluster) {
  InsertionMap insertions;
  for (size_t i = 0; i < reads_in_cluster.size(); ++i) {
    const AlignmentPtr& alignment = reads_in_cluster[i];
    CollectInsertionsForRead(alignment, i, insertions, reads_in_cluster.size());
  }
  AlignInsertions(insertions);
  return insertions;
}

ConsensusMatrix BuildConsensusMatrix(const InsertionMap& aligned_insertions,
                                     const vec<AlignmentPtr>& reads_in_cluster,
                                     const bool hd_deconvolve_enabled,
                                     const bool trim_overhangs) {
  // Step 1: Pre-calculate the consensus length and initialize the matrix
  // Find the leftmost and rightmost reference positions across all alignments
  // in the original reference coordinate
  const u32 leftmost_rpos = LeftmostRefPos(reads_in_cluster, trim_overhangs);
  const u32 rightmost_rpos = RightmostRefPos(reads_in_cluster, trim_overhangs);

  // consensus_matrix_width is the width of the consensus matrix,
  // which is the sum of the lengths of all gapped aligned insertions
  // plus the length of the reference region covered by the longest alignment
  const size_t consensus_matrix_width = ConsensusMatrixWidth(aligned_insertions, leftmost_rpos, rightmost_rpos);
  const size_t max_consensus_length = MaxConsensusLength(reads_in_cluster, aligned_insertions);
  if (consensus_matrix_width > max_consensus_length) {
    throw error::Error(
        "Calculated consensus matrix width {} exceeds maximum possible width {}. "
        "This usually indicates an issue with clustering or duplex read decoding. "
        "Cannot proceed with consensus calling.",
        consensus_matrix_width,
        max_consensus_length);
  }

  // Initialize the consensus matrix
  ConsensusMatrix consensus_matrix(reads_in_cluster.size(), consensus_matrix_width);

  // Step 2: Initialize base types if needed
  // If the base type vector of the first read in the cluster is not empty, we can initialize the base type vector for
  // the consensus matrix. Otherwise, we will just keep the base type vector empty to imply that all bases are to be
  // treated as simplex.
  if (!reads_in_cluster.empty() && hd_deconvolve_enabled) {
    consensus_matrix.SetBaseTypes(
        consensus_matrix.GetReadCount(), consensus_matrix.GetConsensusLength(), yc_decode::BaseType::kSimplex);
  }

  // Step 3: Fill the consensus matrix
  size_t read_index = 0;
  for (const auto& read : reads_in_cluster) {
    Positions positions(read_index, read->StartPos(), aligned_insertions, leftmost_rpos, rightmost_rpos);
    // Set the strand for the read, preprocess the consensus matrix, and collect metadata that will be used for
    // read name generation
    PreprocessConsensusMatrix(consensus_matrix, read_index, read);
    // Fill the consensus matrix for the current read
    FillConsensusMatrix(consensus_matrix, read, positions, aligned_insertions, hd_deconvolve_enabled);
    ++read_index;
  }
  return consensus_matrix;
}

ConsensusMatrix BuildConsensusMatrix(const std::vector<AlignmentPtr>& reads_in_cluster,
                                     const bool hd_deconvolve_enabled,
                                     const bool trim_overhangs) {
  const auto aligned_insertions = GetAlignedInsertions(reads_in_cluster);

  // Step 2: Delegate to the other BuildConsensusMatrix function
  return BuildConsensusMatrix(aligned_insertions, reads_in_cluster, hd_deconvolve_enabled, trim_overhangs);
}

ConsensusMatrix BuildConsensusMatrix(const vec<std::string>& sequences, const vec<AlignmentPtr>& reads_in_cluster) {
  if (sequences.empty()) {
    return ConsensusMatrix{};
  }
  if (reads_in_cluster.size() != sequences.size()) {
    throw error::Error("Number of sequences does not match number of alignments in the cluster");
  }
  PoaAligner msa_aligner(kMatchScore, kMismatchPenalty, kGapOpenPenalty, kGapExtPenalty);
  auto msa_result = msa_aligner.ComputeMSA(sequences);
  ConsensusMatrix consensus_matrix(sequences.size(), msa_result.front().length());
  for (size_t read_index = 0; read_index < sequences.size(); ++read_index) {
    PreprocessConsensusMatrix(consensus_matrix, read_index, reads_in_cluster[read_index]);
    const auto& seq = msa_result[read_index];
    size_t pos = 0;
    size_t last_pos = seq.size();
    // If the read is partial with missing UMI, then it is likely truncated at the start
    // or the end, so we need to skip leading and trailing gaps and instead fill them with 'P'
    if (reads_in_cluster[read_index]->IsPartial()) {
      const size_t first_non_gap_base = seq.find_first_not_of(kBaseGap);
      const size_t last_non_gap_base = seq.find_last_not_of(kBaseGap);
      pos = first_non_gap_base != std::string::npos ? first_non_gap_base : 0;
      last_pos = last_non_gap_base != std::string::npos ? last_non_gap_base + 1 : 0;
    }
    for (; pos < last_pos; ++pos) {
      consensus_matrix.SetBase(read_index, pos, seq[pos]);
    }
  }
  return consensus_matrix;
}

}  // namespace xoos::read_collapser
