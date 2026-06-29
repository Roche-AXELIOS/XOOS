#include "alignment-util.h"

#include <xoos/error/error.h>

namespace xoos::demux {

/**
 * @brief Tries to left-align the indels in the alignment.
 *
 * This algorithm does not change the consensus sequence. It works by iterating over the alignment backwards
 * and moving the indel to the left if it does not cause a mismatch. This algorithm assumes Levenshtein distance
 * (edit distance mismatches, insertions, and deletions are all equal). It would not work properly if we use
 * something that has gap penalties. This works over the entire consensus region, including the adapter trimmed region.
 * Thus, it is possible to move indels into the end adapter trimmed region (effectively removing indels from the read).
 *
 * This function does not move indels in a way that will allow for gap openings to occur, which means that for
 * groups of colliding indels it is possible to collapse gaps together. This is probably desirable behavior as errors
 * are often clustered together in the same region.
 *
 * @param consensus The consensus sequence. Only the duplex region.
 * @param alignment The alignment data to be modified. Must be at least as long as the consensus.
 */
void LeftAlignConvert(std::string_view consensus, u8* alignment) {
  u32 indel_length = 0;

  // Iterate over the alignment backwards and reverse the indel positions, Basic idea:
  // Iterate backwards until we find an indel
  // Once found we set the indel length to 1
  // If we encounter another indel increment the length
  // If we encounter a non-indel we at this we can begin swapping alignment indel positions
  // If the swapped is invalid (i.e. causes a mismatch) we reset the indel length to 0
  for (u32 consensus_pos = consensus.length(); consensus_pos-- > 0;) {
    bool is_indel = false;
    switch (alignment[consensus_pos]) {
      case kEdopInsert:
      case kEdopDelete:
        // we found an indel so we increment the length
        ++indel_length;
        is_indel = true;
      default:
        break;
    }
    // non indel state
    if (!is_indel && indel_length > 0) {
      // our last position was an indel so enter into the inner loop
      u32 indel_pos = consensus_pos + indel_length;  // furthest right position of the indel
      if (consensus.at(consensus_pos) == consensus.at(indel_pos)) {
        // if the current base is the same as the indel base we swap them
        std::swap(alignment[indel_pos], alignment[consensus_pos]);
      } else {
        indel_length = 0;
      }
    }
  }
}

/**
 * This function groups the indels by performing left alignment followed by right alignment so it is symmetric to left
 * alignment. It works by first left aligning indels in the alignment and then reversing the grouped indels.
 *
 * To be clear, this function does not change the consensus sequence it just modifies the alignment data.
 *
 * This the idea is to make both strands to be flipped the same way and this seems to improve the positioning of quality
 * bases locations, empirically suggesting an affine gap penalty based alignment would work better than Levenshtein
 * distance
 *
 * @param consensus The consensus sequence. Only the duplex region.
 * @param alignment The alignment data to be modified. Must be at least as long as the consensus.
 */
void GroupRightAlignConvert(std::string_view consensus, u8* alignment) {
  LeftAlignConvert(consensus, alignment);
  // reverse the grouped indels
  u32 indel_length = 0;
  for (u32 consensus_pos = 0; consensus_pos < consensus.length(); ++consensus_pos) {
    bool is_indel = false;
    switch (alignment[consensus_pos]) {
      case kEdopInsert:
      case kEdopDelete:
        // we found an indel so we increment the length
        ++indel_length;
        is_indel = true;
      default:
        break;
    }
    // non-indel state
    if (!is_indel && indel_length > 0) {
      // our last position was an indel so enter into the inner loop
      u32 indel_pos = consensus_pos - indel_length;  // furthest left position of the indel
      if (consensus.at(consensus_pos) == consensus.at(indel_pos)) {
        // if the current base is the same as the indel base we swap them
        std::swap(alignment[indel_pos], alignment[consensus_pos]);
      } else {
        indel_length = 0;
      }
    }
  }
}

/**
 * @brief Performs methylation conversion on the consensus and alignment based on the R2 read.
 *
 * This function scans through the alignment 2 bases at a time looking for specific patterns that indicate
 * methylation status. When such patterns are found, it converts the mismatches back to non-methylated matches
 * (C->T or G->A) in the consensus and updates the alignment accordingly. It also updates the xm tag to reflect
 * the methylation status.
 *
 * @param consensus The consensus sequence to be modified (duplex region only, post endadapter duplex region only).
 * @param alignment The alignment data to be modified. (post endadapter trimming start)
 * @param alignment_length The length of the alignment. (post endadapter trimming shortened)
 * @param xm_tag_start The start of the xm tag to be updated.
 * @param r2 The R2 read sequence.
 * @param swapped A boolean indicating if R2 is larger than R1.
 * @return The number of converted methylation sites.
 */
void MethylConvert(char* consensus, u8* alignment, s32 alignment_length, char* xm_tag_start, const u8* r2,
                   bool swapped) {
  // xm tag starts after the simplex region and is assumed to be initialized to all '.'
  // perform methylation conversion on the alignment by scanning through the alignment 2 bases at a time looking for:
  // 1. TG on the forward strand and CA on the reverse strand causing 2 mismatches
  // 2. or TG on the forward strand and CG on the reverse strand causing a mismatch (partial methylation)
  // 3. or CG on the forward strand and CA on the reverse strand causing a mismatch (partial methylation)
  // When we find the above patterns we convert the mismatches back to non-methylated matches (C->T or G->A)
  // e.g. if old consensus is ATGT (methylated C->T) we convert back to ACGT (non-methylated C)
  // Then we update the xm tag to
  // Z (R1 and R2 methylated)
  // U (R1 methylated only)
  // u (R2 methylated only) or
  // z (unmethylated)
  // e.g. R1: TATGTTCGTCGGTG
  //      R2: TACATTCGTCAGCG
  // Consens: TACGTTCGTCGGAG
  // XM_Tag : ..Z...z..u..U.
  // XM tag will be simplex_length + consensus.length() long

  if (alignment_length < 2) {
    // the consensus is too short to have CpG sites
    return;
  }

  s32 r2_pos = 0;

  for (s32 consensus_pos = 0; consensus_pos < alignment_length - 1; ++consensus_pos) {
    u8& op1 = alignment[consensus_pos];
    u8& op2 = alignment[consensus_pos + 1];

    // if either operation is a deletion or insertion we skip this position
    if (op1 != kEdopDelete && op1 != kEdopInsert && op2 != kEdopDelete && op2 != kEdopInsert) {
      char& base1 = consensus[consensus_pos];
      // shouldn't need to be changed since we currently use R1 as the base if mismatches occur
      // may need to change if we switch to using the higher quality read as the base
      const char base2 = consensus[consensus_pos + 1];
      const char r2_base1 = kAlphabet[r2[r2_pos]];
      const char r2_base2 = kAlphabet[r2[r2_pos + 1]];

      if (op1 == kEdopMatch && op2 == kEdopMatch && base1 == 'C' && base2 == 'G') {
        // CG site but unmethylated (z)
        xm_tag_start[consensus_pos] = 'z';
      } else if (op1 == kEdopMismatch && op2 == kEdopMismatch && base1 == 'T' && base2 == 'G' && r2_base1 == 'C' &&
                 r2_base2 == 'A') {
        // CpG site and fully methylated (Z)
        // Correct the consensus
        base1 = 'C';
        // update xm tag
        xm_tag_start[consensus_pos] = 'Z';
        // update the alignment to match (but methylated)
        op1 = kEdopMatchMethyl;
        op2 = kEdopMatchMethyl;
      } else if (op1 == kEdopMismatch && op2 == kEdopMatch && base1 == 'T' && base2 == 'G' && r2_base1 == 'C' &&
                 r2_base2 == 'G') {
        // CpG site and hemi-methylated on R1 (U)
        // Correct the consensus
        base1 = 'C';
        // update xm tag
        xm_tag_start[consensus_pos] = 'U';
        // update the alignment to match (but methylated)
        op1 = kEdopMatchMethyl;
      } else if (op1 == kEdopMatch && op2 == kEdopMismatch && base1 == 'C' && base2 == 'G' && r2_base1 == 'C' &&
                 r2_base2 == 'A') {
        // CpG site and hemi-methylated on R2 (u)
        // Consensus uses R1 so does not need to change base
        // update xm tag
        xm_tag_start[consensus_pos] = 'u';
        // update the alignment to match (but methylated)
        op2 = kEdopMatchMethyl;
      }
    }
    if (swapped) {
      // R2 larger than R1 so it becomes the reference, so we need to increment when we see deletions
      switch (op1) {
        case kEdopMatchMethyl:
        case kEdopMismatch:
        case kEdopMatch:
        case kEdopDelete:
          ++r2_pos;
          break;
        case kEdopInsert:
          // deletion in consensus means we do not move r2
          break;
        default:
          throw error::Error("Invalid alignment operation code encountered during methylation conversion.");
      }
    } else {
      // R2 is smaller than R1 so it becomes the query, so we need to increment when we see insertions
      switch (op1) {
        case kEdopMatchMethyl:
        case kEdopMismatch:
        case kEdopMatch:
        case kEdopInsert:
          ++r2_pos;
          break;
        case kEdopDelete:
          // deletion in consensus means we do not move r2
          break;
        default:
          throw error::Error("Invalid alignment operation code encountered during methylation conversion.");
      }
    }
  }
}
}  // namespace xoos::demux
