#include "task/alignment.h"

#include <fmt/format.h>
#include <xoos/error/error.h>
#include <xoos/log/logging.h>

#include <algorithm>
#include <utility>

#include "metrics/duplex-metrics.h"
#include "task/flow-context.h"
#include "task/flow-manager.h"
#include "utility/alignment-util.h"
#include "utility/math-util.h"
#include "utility/stop-watch.h"

/*
 * For alignment, we are using a heavily modified version of Martin Sosic's Edlib library
 (https://github.com/Martinsos/edlib).
 * That software is licensed under the MIT License, which is included below per license terms:

    The MIT License (MIT)

    Copyright (c) 2014 Martin Šošić

    Permission is hereby granted, free of charge, to any person obtaining a copy of
    this software and associated documentation files (the "Software"), to deal in
    the Software without restriction, including without limitation the rights to
    use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
    the Software, and to permit persons to whom the Software is furnished to do so,
    subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
    FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
    COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
    IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 * We initially modified various parts of the library to better suit Roche's
 * needs, including removing unneeded code (notably Hirschberg support as our reads are short), removing support for
 * alignment modes that we don't need (MODE_HW and MODE_SW) as that simplifies the code and makes it faster, replacing
 * heap allocation by stack allocation, massive renaming of variables/functions to satisfy Roche's coding standards,
 * converting the various static functions into a class (so we can pass around state), using a fixed alphabet
 * (thus eliminating the overhead of creating one on-the-fly) and other speed optimizations.
 *
 * Code has continued to undergo various modifications to better conform to Roche's coding standards and implement new
 * features as needed.
 *
 * The purpose of this code is to align the R1 and R2 regions of a duplex read to create a consensus sequence.
 * When aligning, we'll determine edit distance and alignment path (and start and end locations of it in reference).
 *
 * We support both prefix and infix alignment modes. Infix penalizes gaps at the beginning and end of the query
 * sequence, but prefix does not. Prefix alignment is approximately 60% faster than infix alignment. In theory using the
 * prefix alignment mode, where the gap at the query end is not penalized, makes sense on partial duplex reads because
 * the "simplex" region should not align. For full duplex reads (expected that R1 and R2 both have endadapters) it
 * may make more sense to use infix alignment, though if the gap penalty at the end of the alignment is not important,
 * then prefix alignment might better due to speed.
 *
 * TODO: We currently default to using infix alignment. We should test prefix alignment.
*/

namespace xoos::demux {

/**
 * Sets the sequences to be aligned based on the read and trimming information.
 * @param [in] read  The read containing the sequences and trimming info.
 * @return true if sequences were set successfully, false if sequence is too short
 */
bool Alignment::SetSequences(const FixedReadRecord& read) {
  // assumes midadapter was found; we will perform alignment using data for which the midadapter has been trimmed out,
  // but endadapters are still present.  The midadapter is the SID region of R1 + loop + SID region of R2.

  // do not trim out the start adapter
  insert_start1 = 0;
  // INCLUSIVE
  insert_end1 = static_cast<s32>(read.trim_info_duplex.midadapter_range->spos - 1);
  insert_start2 = static_cast<s32>(read.trim_info_duplex.midadapter_range->epos + 1);
  // this assumes endadapter is not found.
  // INCLUSIVE
  insert_end2 = static_cast<s32>(read.SeqLen() - 1);

  reference_length = insert_end1 - insert_start1 + 1;
  query_length = insert_end2 - insert_start2 + 1;

  // In edlib, the shorter sequence is the query. If the query is longer, we swap the two sequences.
  is_swapped = query_length > reference_length;

  if (is_swapped) {  // need to swap the sequences to make query the shortest
    std::swap(query, reference);
    std::swap(insert_start1, insert_start2);
    std::swap(insert_end1, insert_end2);
    std::swap(query_length, reference_length);
  }
  // Do not attempt alignment if the insert is too short
  if (reference_length + query_length <= ToSigned(_params.min_trimmed_read_len)) {
    return false;
  }

  /*------------ TRANSFORM SEQUENCES INTO ALPHABET -----------*/
  // Calculate reverse and complement flag for the first sequence, the second sequence has the inverse values.
  bool do_reverse = is_prefix;

  // default for first read
  bool do_complement = false;

  if (is_swapped) {  // invert the values if we have swapped the strings (i.e. R2 was the longest)
    do_reverse = !do_reverse;
    do_complement = !do_complement;
  }
  simd::TransformSequence(read.Seq(), insert_start1, insert_end1, reference, do_reverse, do_complement);
  simd::TransformSequence(read.Seq(), insert_start2, insert_end2, query, !do_reverse, !do_complement);

  return true;
}

/**
 * Main edlib method.
 * @return true if alignment was successful, false if something went wrong (e.g. alignment length is too long).
 */
bool Alignment::Align() {
  // Do not attempt alignment if the insert is too long
  if (static_cast<uint>(reference_length) > kMaxAlignLength) {
    return false;
  }

  edit_distance = -1;
  num_locations = 0;
  alignment_length = 0;

  /*--------------------- INITIALIZATION ------------------*/
  // bmax in Myers
  const s32 max_num_blocks = CeilDiv(query_length, kWordSize);
  // number of redundant cells in last level blocks
  const s32 num_redundant_cells = max_num_blocks * kWordSize - query_length;

  simd::BuildPeq(query, query_length, peq_table);

  /*------------------ MAIN CALCULATION -------------------*/
  // only allow a specific error rate for the alignment. This means that we can calculate the maximum
  // edit distance that we'd allow - and smaller edit distances will be calculated faster.
  constexpr float kPercentToDecimal = 0.01f;

  s32 k =
      1 + std::max(0, static_cast<s32>((max_error_rate * static_cast<float>(std::min(query_length, reference_length))) *
                                       kPercentToDecimal));

  MyersCalcEditDistanceSemiGlobal(peq_table, num_redundant_cells, max_num_blocks, k, is_prefix, reference,
                                  reference_length, edit_distance, end_locations, num_locations);

  if (edit_distance >= 0) {
    // If there is solution.

    // We need to penalize the edit distance if we have excess bases in the reference sequence after alignment.
    // TODO: This is a temporary solution, we should use a better alignment methodology
    const auto reference_tail_length = reference_length - 1 - end_locations[num_locations - 1];
    if ((edit_distance + reference_tail_length) > k) {
      // alignment is higher error rate than allowed
      edit_distance = -1;
    } else {
      // Initialize starting locations.
      if (!is_prefix) {  // we're using infix (MODE_HW) alignment
        simd::ReverseCopy(reference, reference + reference_length, reverse_reference);
        simd::ReverseCopy(query, query + query_length, reverse_query);
        simd::BuildPeq(reverse_query, query_length, reverse_peq_table);
        for (s32 i = 0; i < num_locations; ++i) {
          const s32 end_location = end_locations[i];
          if (end_location == -1) {
            start_locations[i] = 0;
          } else {
            // score is not used here
            s32 score;
            s32 num_reverse_positions;
            MyersCalcEditDistanceSemiGlobal(reverse_peq_table, num_redundant_cells, max_num_blocks, k, true,
                                            reverse_reference + reference_length - end_location - 1, end_location + 1,
                                            score, reverse_positions, num_reverse_positions);
            start_locations[i] = end_location - reverse_positions[num_reverse_positions - 1];
          }
        }
      } else {
        for (s32 i = 0; i < num_locations; ++i) {
          start_locations[i] = 0;
        }
      }

      // Find alignment
      // For infix alignment we want the alignment closest to the end of the reference if there are multiple alignments.
      const s32 aln_start = start_locations[num_locations - 1];
      const s32 aln_end = end_locations[num_locations - 1];
      const u8* aln_reference = reference + aln_start;
      const s32 aln_reference_length = aln_end - aln_start + 1;
      simd::ReverseCopy(reference, reference + reference_length, reverse_aln_reference);
      simd::ReverseCopy(query, query + query_length, reverse_aln_query);
      ObtainAlignment(aln_reference, aln_reference_length);
      // This shouldn't happen much as we check kMaxAlignLength against reference length before this.
      // The inserted bases after alignment cause this to happen.
      if (static_cast<uint>(start_locations[num_locations - 1] + alignment_length + reference_tail_length) >
          kMaxAlignLength) {
        // don't need to adjust edit_distance, we are stating here that alignment failed due to too long consensus
        return false;
      }
    }
  }
  return true;
}

/**
 * Corresponds to Advance_Block function from Myers.
 * Calculates one word(block), which is part of a column.
 * Highest bit of word (one most to the left) is most bottom cell of block from column.
 * pos_v[i] and minus_v[i] define vin of cell[i]: vin = cell[i] - cell[i-1].
 * @param [in] pos_v  Bitset, pos_v[i] == 1 if vin is +1, otherwise pos_v[i] == 0.
 * @param [in] minus_v  Bitset, minus_v[i] == 1 if vin is -1, otherwise minus_v[i] == 0.
 * @param [in] eq  Bitset, eq[i] == 1 if match, 0 if mismatch.
 * @param [in] pos_hin  Bitset, 00..01 when hin == +1.
 * @param [in] minus_hin  Bitset, 00..01 when hin == -1.
 * @param [out] pos_hout  Bitset, pos_hout[i] == 1 if vout is +1, otherwise pos_hout[i] == 0.
 * @param [out] minus_hout  Bitset, minus_hout[i] == 1 if vout is -1, otherwise minus_hout[i] == 0.
 *
 * Note: this is a modification introduced by Ragnar Groot Koerkamp to allow for the calculation of the hin and hout
 * using +1/-1 indicators
 * (https://github.com/Martinsos/edlib/pull/214/commits/f2a465332c1bd240b528eaf72b7e561b2e809e86). This gave a 5%
 * performance improvement on overall throughput.
 */
inline Delta Alignment::CalculateBlock(const Word pos_v, const Word minus_v, Word eq, const Word pos_hin,
                                       const Word minus_hin, Word& pos_hout, Word& minus_hout) {
  // hin can be 1, -1 or 0.
  // 1  -> 00...01
  // 0  -> 00...00
  // -1 -> 11...11 (2-complement)

  const Word xv = eq | minus_v;
  // This is instruction below written using 'if': if (hin < 0) eq |= (Word)1;
  eq |= minus_hin;
  const Word xn = (((eq & pos_v) + pos_v) ^ pos_v) | eq;

  Word ph = minus_v | ~(xn | pos_v);
  Word mh = pos_v & xn;

  // This is instruction below written using 'if': if (ph & kHighBitMask) hout = 1;
  const auto p_hout = static_cast<s32>(ph >> (kWordSize - 1));
  // This is instruction below written using 'if': if (mh & kHighBitMask) hout = -1;
  const auto m_hout = static_cast<s32>(mh >> (kWordSize - 1));

  // used to be ph <<= 1 -> this is equivalent but usually faster
  ph += ph;
  mh += mh;

  // This is instruction below written using 'if': if (hin < 0) mh |= (Word)1;
  mh |= minus_hin;
  // This is instruction below written using 'if': if (hin > 0) ph |= (Word)1;
  ph |= pos_hin;

  pos_hout = mh | ~(xv | ph);
  minus_hout = ph & xv;

  return {p_hout, m_hout};
}

static inline s32 Min(const s32 x, const s32 y) { return x < y ? x : y; }

static inline s32 Max(const s32 x, const s32 y) { return x > y ? x : y; }

/**
 * @param [in] block
 * @param [in/out] scores array with scores
 * @return Values of cells in block, starting with bottom cell in block.
 */
void Alignment::GetBlockCellValues(const Block& block, s32* scores) {
  s32 score = block.score;
  Word mask = kHighBitMask;
  for (s32 i = 0; i < kWordSize - 1; ++i) {
    scores[i] = score;
    if (block.p & mask) {
      --score;
    }
    if (block.m & mask) {
      ++score;
    }
    mask >>= 1;
  }
  scores[kWordSize - 1] = score;
}

/**
 * Uses Myers' bit-vector algorithm to find edit distance for one of semi-global alignment methods.
 * @param [in] w  Size of padding in last block.
 * @param [in] max_num_blocks  Number of blocks needed to cover the whole query.
 * @param [in] k
 */
void Alignment::MyersCalcEditDistanceSemiGlobal(const Word peq_input[], const s32 w, const s32 max_num_blocks, s32 k,
                                                const bool do_prefix, const u8* target, const s32 target_length,
                                                s32& score, s32* positions_out, s32& num_positions_out) {
  num_positions_out = 0;

  // first_block is 0-based index of first block in Ukkonen band.
  // last_block is 0-based index of last block in Ukkonen band.
  s32 first_block = 0;
  // y in Myers
  s32 last_block = Min(CeilDiv(k + 1, kWordSize), max_num_blocks) - 1;
  Block blocks[kMaxNrBlocks];

  // For HW, solution will never be larger than query_length.
  k = Min(query_length, k);

  // Initialize P, M and score
  for (s32 b = 0; b <= last_block; ++b) {
    blocks[b].score = (b + 1) * kWordSize;
    // All 1s
    blocks[b].p = static_cast<Word>(-1);
    blocks[b].m = static_cast<Word>(0);
  }

  s32 best_score = -1;
  s32 nr_positions = 0;

  // gap before query is penalized in prefix mode
  const Word pos_start_hout = do_prefix ? 1 : 0;
  Block* bl{nullptr};
  const u8* reference_char = target;
  for (s32 c = 0; c < target_length; ++c) {  // for each column
    const Word* peq_c = peq_input + (*reference_char) * max_num_blocks;

    //----------------------- Calculate column -------------------------//
    auto pos_hout = static_cast<s32>(pos_start_hout);
    s32 minus_hout = 0;
    s32 hout = 0;
    for (s32 b = first_block; b <= last_block; ++b) {
      auto& block{blocks[b]};
      const auto delta_hout = CalculateBlock(block.p, block.m, peq_c[b], pos_hout, minus_hout, block.p, block.m);
      pos_hout = delta_hout.pos_hout;
      minus_hout = delta_hout.minus_hout;
      hout = pos_hout - minus_hout;
      block.score += hout;
    }
    bl = blocks + last_block;
    peq_c += last_block;

    //---------- Adjust number of blocks according to Ukkonen ----------//
    // bl is pointing to last block
    if ((last_block < max_num_blocks - 1) && (bl->score - hout <= k) && ((*(peq_c + 1) & kWord1) || (hout < 0))) {
      // peq_c is pointing to last block
      // If score of left block is not too big, calculate one more block
      ++last_block;
      ++bl;
      ++peq_c;
      // All 1s
      bl->p = static_cast<Word>(-1);
      bl->m = static_cast<Word>(0);
      const auto delta_hout = CalculateBlock(bl->p, bl->m, *peq_c, pos_hout, minus_hout, bl->p, bl->m);
      bl->score = (bl - 1)->score - hout + kWordSize + delta_hout.pos_hout - delta_hout.minus_hout;
    } else {
      while ((last_block >= first_block) && (bl->score >= (k + kWordSize))) {
        --last_block;
        --bl;
        --peq_c;
      }
    }

    // For infix mode, even if all cells are > k, there still may be a solution in the next column
    // because starting conditions at upper boundary are 0. This means that first block is always a
    // candidate for solution, and we can never end calculation before the last column.
    if (!do_prefix && (last_block == -1)) {
      ++last_block;
      ++bl;
      ++peq_c;
    }

    if (do_prefix) {
      // Reduce band by increasing first block if possible. Not applicable to HW.
      while ((first_block <= last_block) && (blocks[first_block].score >= (k + kWordSize))) {
        ++first_block;
      }
    }

    // If band stops to exist finish
    if (last_block < first_block) {
      score = best_score;
      if (best_score != -1) {
        num_positions_out = nr_positions;
        std::copy_n(positions, nr_positions, positions_out);
      }
      return;
    }

    //------------------------- Update best score ----------------------//
    if (last_block == max_num_blocks - 1) {
      const s32 col_score = bl->score;
      if (col_score <= k) {  // Scores > k dont have correct values (so we cannot use them), but are certainly > k.
        // NOTE: Score that I find in column c is actually score from column c-w
        if ((best_score == -1) || (col_score <= best_score)) {
          if (col_score != best_score) {
            nr_positions = 0;
            best_score = col_score;
            // Change k so we will look only for equal or better
            // scores then the best found so far.
            k = best_score;
          }
          positions[nr_positions++] = c - w;
        }
      }
    }
    ++reference_char;
  }

  // Obtain results for last w columns from last column.
  if (last_block == max_num_blocks - 1) {
    s32 block_scores[kWordSize];
    GetBlockCellValues(*bl, block_scores);
    for (s32 i = 0; i < w; ++i) {
      if (const s32 col_score = block_scores[i + 1];
          (col_score <= k) && (best_score == -1 || (col_score <= best_score))) {
        if (col_score != best_score) {
          nr_positions = 0;
          k = best_score = col_score;
        }
        positions[nr_positions++] = target_length - w + i;
      }
    }
  }

  score = best_score;
  if (best_score != -1) {
    num_positions_out = nr_positions;
    std::copy_n(positions, nr_positions, positions_out);
  }
}

/**
 * Uses Myers' bit-vector algorithm to find edit distance for global(NW) alignment method.
 * @param [in] w  Size of padding in last block.
 * @param [in] max_num_blocks  Number of blocks needed to cover the whole query.
 * @param [in] align_reference
 * @param [in] align_reference_length
 */
void Alignment::MyersCalcEditDistanceNW(const s32 w, const s32 max_num_blocks, const u8* align_reference,
                                        s32 align_reference_length) {
  s32 k = edit_distance;
  if (k < ::abs(align_reference_length - query_length)) {
    return;
  }

  // Upper bound for k
  k = Min(k, Max(query_length, align_reference_length));

  // first_block is 0-based index of first block in Ukkonen band.
  // last_block is 0-based index of last block in Ukkonen band.
  s32 first_block = 0;
  // This is optimal now, by my formula.
  s32 last_block =
      Min(max_num_blocks, CeilDiv(Min(k, (k + query_length - align_reference_length) / 2) + 1, kWordSize)) - 1;

  Block blocks[kMaxNrBlocks];

  // Initialize P, M and score
  for (s32 b = 0; b <= last_block; ++b) {
    blocks[b].score = (b + 1) * kWordSize;
    // All 1s
    blocks[b].p = static_cast<Word>(-1);
    blocks[b].m = static_cast<Word>(0);
  }

  const u8* reference_char = align_reference;
  for (s32 c = 0; c < align_reference_length; ++c) {  // for each column
    const Word* peq_c = peq_table + *reference_char * max_num_blocks;

    //----------------------- Calculate column -------------------------//
    s32 hout = 1;
    s32 pos_hout = 1;
    s32 minus_hout = 0;
    // Unroll loop to process two blocks at a time
    for (s32 b = first_block; b <= last_block; b += 2) {
      auto& block1 = blocks[b];
      const auto delta_hout1 = CalculateBlock(block1.p, block1.m, peq_c[b], pos_hout, minus_hout, block1.p, block1.m);
      pos_hout = delta_hout1.pos_hout;
      minus_hout = delta_hout1.minus_hout;
      hout = pos_hout - minus_hout;
      block1.score += hout;

      if (b + 1 <= last_block) {
        auto& block2 = blocks[b + 1];
        const auto delta_hout2 =
            CalculateBlock(block2.p, block2.m, peq_c[b + 1], pos_hout, minus_hout, block2.p, block2.m);
        pos_hout = delta_hout2.pos_hout;
        minus_hout = delta_hout2.minus_hout;
        hout = pos_hout - minus_hout;
        block2.score += hout;
      }
    }

    auto* bl = &blocks[last_block];
    //------------------------------------------------------------------//
    // bl now points to last block

    // Update k. I do it only on end of column because it would slow calculation too much otherwise.
    // NOTICE: I add w when in last block because it is actually result from w cells to the left and w cells up.
    k = Min(k, bl->score + Max(align_reference_length - c - 1, query_length - ((1 + last_block) * kWordSize - 1) - 1) +
                   (last_block == max_num_blocks - 1 ? w : 0));

    //---------- Adjust number of blocks according to Ukkonen ----------//
    //--- Adjust last block ---//
    // If block is not beneath band, calculate next block. Only next because others are certainly beneath band.
    if (((last_block + 1) < max_num_blocks) &&
        !(  // score[last_block] >= k + kWordSize ||  // NOTICE: this condition could be satisfied if above block also!
            ((last_block + 1) * kWordSize - 1) >
            (k - bl->score + 2 * kWordSize - 2 - align_reference_length + c + query_length))) {
      ++last_block;
      ++bl;
      // All 1s
      bl->p = static_cast<Word>(-1);
      bl->m = static_cast<Word>(0);
      const auto delta_hout = CalculateBlock(bl->p, bl->m, peq_c[last_block], pos_hout, minus_hout, bl->p, bl->m);
      pos_hout = delta_hout.pos_hout;
      minus_hout = delta_hout.minus_hout;
      const s32 new_hout = pos_hout - minus_hout;
      bl->score = (bl - 1)->score - hout + kWordSize + new_hout;
      hout = new_hout;
    }

    // While block is out of band, move one block up.
    // NOTE: Condition used here is more loose than the one from the article, since I simplified the Max() part of it.
    // I could consider adding that Max part, for optimal performance.
    while ((last_block >= first_block) &&
           ((bl->score >= (k + kWordSize)) ||
            (((last_block + 1) * kWordSize - 1) >
             (k - bl->score + 2 * kWordSize - 2 - align_reference_length + c + query_length + 1)))) {
      --last_block;
      --bl;
    }

    //--- Adjust first block ---//
    // While outside of band, advance block
    while ((first_block <= last_block) &&
           ((blocks[first_block].score >= (k + kWordSize)) ||
            (((first_block + 1) * kWordSize - 1) <
             (blocks[first_block].score - k - align_reference_length + query_length + c)))) {
      ++first_block;
    }

    // If band stops to exist finish
    if (last_block < first_block) {
      return;
    }
    //------------------------------------------------------------------//

    //---- Save column so it can be used for reconstruction ----//
    if (c < align_reference_length) {
      for (s32 b = first_block; b <= last_block; ++b) {
        Ps[max_num_blocks * c + b] = blocks[b].p;
        Ms[max_num_blocks * c + b] = blocks[b].m;
        scores[max_num_blocks * c + b] = blocks[b].score;
      }
      first_blocks[c] = first_block;
      last_blocks[c] = last_block;
    }
    //----------------------------------------------------------//
    //---- If this is stop column, save it and finish ----//
    if (c == -1) {
      for (s32 b = first_block; b <= last_block; ++b) {
        Ps[b] = blocks[b].p;
        Ms[b] = blocks[b].m;
        scores[b] = blocks[b].score;
      }
      first_blocks[0] = first_block;
      last_blocks[0] = last_block;
      return;
    }
    //----------------------------------------------------//

    ++reference_char;
  }

  if (last_block == max_num_blocks - 1) {  // If last block of last column was calculated
    // Obtain best score from block -> it is complicated because query is padded with w cells
    s32 block_scores[kWordSize];
    GetBlockCellValues(blocks[last_block], block_scores);
    const s32 best_score_tmp = block_scores[w];
    if (best_score_tmp <= k) {
      return;
    }
  }
}

/**
 * Finds one possible alignment that gives optimal score by moving back through the dynamic programming matrix,
 * that is stored in  Consumes large amount of memory: O(query_length * aligned_reference_length).
 * @param [in] aligned_reference_length  Normal length, without W.
 * @return Status code.
 */
void Alignment::ObtainAlignmentTraceback(const s32 aligned_reference_length) {
  const auto best_score = edit_distance;
  const s32 max_num_blocks = CeilDiv(query_length, kWordSize);
  const s32 w = max_num_blocks * kWordSize - query_length;
  alignment_length = 0;
  // index of column
  s32 c = aligned_reference_length - 1;
  // index of block in column
  s32 b = max_num_blocks - 1;
  // Score of current cell
  s32 current_score = best_score;
  // Score of left cell
  s32 left_score = -1;
  // Score of upper cell
  s32 upper_score = -1;
  // Score of upper left cell
  s32 upper_left_score = -1;
  // P of current block
  Word curr_p = Ps[c * max_num_blocks + b];
  // M of current block
  Word curr_m = Ms[c * max_num_blocks + b];
  // True if block to left exists and is in band
  bool there_is_left_block = (c > 0) && (b >= first_blocks[c - 1]) && (b <= last_blocks[c - 1]);
  // We set initial values of l_p and l_m to 0 only to avoid compiler warnings, they should not affect the
  // calculation as both l_p and l_m should be initialized at some moment later (but compiler can not
  // detect it since this initialization is guaranteed by "business" logic).
  Word l_p = 0;
  Word l_m = 0;
  if (there_is_left_block) {
    // P of block to the left
    l_p = Ps[(c - 1) * max_num_blocks + b];
    // M of block to the left
    l_m = Ms[(c - 1) * max_num_blocks + b];
  }
  curr_p <<= w;
  curr_m <<= w;
  // 0 based index of current cell in block_pos
  s32 block_pos = kWordSize - w - 1;

  while (true) {
    if (c == 0) {
      there_is_left_block = true;
      left_score = b * kWordSize + block_pos + 1;
      upper_left_score = left_score - 1;
    }

    //---------- Calculate scores ---------//
    if ((left_score == -1) && there_is_left_block) {
      // score of block to the left
      left_score = scores[(c - 1) * max_num_blocks + b];

      // Modified code (Roche)
      const s32 nr_bits{kWordSize - block_pos - 1};
      if (nr_bits) {
        // The above loop is equivalent to a popcount of the highest bits of l_p and l_m.
        const auto mask{~((1UL << (block_pos + 1)) - 1UL)};
        left_score -= std::popcount(l_p & mask);
        left_score += std::popcount(l_m & mask);
        l_p <<= nr_bits;
        l_m <<= nr_bits;
      }
    }
    if (upper_left_score == -1) {
      if (left_score != -1) {
        upper_left_score = left_score;
        if (l_p & kHighBitMask) {
          --upper_left_score;
        }
        if (l_m & kHighBitMask) {
          ++upper_left_score;
        }
      } else if ((c > 0) && ((b - 1) >= first_blocks[c - 1]) && ((b - 1) <= last_blocks[c - 1])) {
        // This is the case when upper left cell is last cell in block,
        // and block to left is not in band so left_score is -1.
        upper_left_score = scores[(c - 1) * max_num_blocks + b - 1];
      }
    }
    if (upper_score == -1) {
      upper_score = current_score;
      if (curr_p & kHighBitMask) {
        --upper_score;
      }
      if (curr_m & kHighBitMask) {
        ++upper_score;
      }
      curr_p <<= 1;
      curr_m <<= 1;
    }

    //-------------- Move --------------//
    // Move up - insertion to reference - deletion from query
    if ((upper_score != -1) && ((upper_score + 1) == current_score)) {
      current_score = upper_score;
      left_score = upper_left_score;
      upper_score = upper_left_score = -1;
      if (block_pos == 0) {  // If entering new (upper) block
        if (b == 0) {        // If there are no cells above (only boundary cells)
          // Move up
          alignment[alignment_length++] = kEdopInsert;
          for (s32 i = 0; i < c + 1; ++i) {  // Move left until end
            alignment[alignment_length++] = kEdopDelete;
          }
          break;
        } else {
          block_pos = kWordSize - 1;
          --b;
          curr_p = Ps[c * max_num_blocks + b];
          curr_m = Ms[c * max_num_blocks + b];
          if ((c > 0) && (b >= first_blocks[c - 1]) && (b <= last_blocks[c - 1])) {
            there_is_left_block = true;
            l_m = Ms[(c - 1) * max_num_blocks + b];
            l_p = Ps[(c - 1) * max_num_blocks + b];
          } else {
            there_is_left_block = false;
          }
        }
      } else {
        --block_pos;
        l_p <<= 1;
        l_m <<= 1;
      }
      // Mark move
      alignment[alignment_length++] = kEdopInsert;
    } else {
      // Move left - deletion from reference - insertion to query
      if ((left_score != -1) && ((left_score + 1) == current_score)) {
        current_score = left_score;
        upper_score = upper_left_score;
        left_score = upper_left_score = -1;
        --c;
        if (c == -1) {  // If there are no cells to the left (only boundary cells)
          // Move left
          alignment[alignment_length++] = kEdopDelete;
          const s32 num_up = b * kWordSize + block_pos + 1;
          for (s32 i = 0; i < num_up; ++i) {  // Move up until end
            alignment[alignment_length++] = kEdopInsert;
          }
          break;
        }
        curr_p = l_p;
        curr_m = l_m;
        if ((c > 0) && (b >= first_blocks[c - 1]) && (b <= last_blocks[c - 1])) {
          there_is_left_block = true;
          l_p = Ps[(c - 1) * max_num_blocks + b];
          l_m = Ms[(c - 1) * max_num_blocks + b];
        } else {
          if (c == 0) {  // If there are no cells to the left (only boundary cells)
            there_is_left_block = true;
            left_score = b * kWordSize + block_pos + 1;
            upper_left_score = left_score - 1;
          } else {
            there_is_left_block = false;
          }
        }
        // Mark move
        alignment[alignment_length++] = kEdopDelete;
      } else {
        // Move up left - (mis)match
        if (upper_left_score != -1) {
          const u8 move_code = upper_left_score == current_score ? kEdopMatch : kEdopMismatch;
          current_score = upper_left_score;
          upper_score = left_score = upper_left_score = -1;
          --c;
          if (c == -1) {  // If there are no cells to the left (only boundary cells)
            // Move left
            alignment[alignment_length++] = move_code;
            const s32 num_up = b * kWordSize + block_pos;
            for (s32 i = 0; i < num_up; ++i) {  // Move up until end
              alignment[alignment_length++] = kEdopInsert;
            }
            break;
          }
          if (block_pos == 0) {  // If entering upper left block
            if (b == 0) {        // If there are no more cells above (only boundary cells)
              // Move up left
              alignment[alignment_length++] = move_code;
              for (s32 i = 0; i < c + 1; ++i) {  // Move left until end
                alignment[alignment_length++] = kEdopDelete;
              }
              break;
            }
            block_pos = kWordSize - 1;
            --b;
            curr_p = Ps[c * max_num_blocks + b];
            curr_m = Ms[c * max_num_blocks + b];
          } else {  // If entering left block
            --block_pos;
            curr_p = l_p;
            curr_m = l_m;
            curr_p <<= 1;
            curr_m <<= 1;
          }
          // Set new left block
          if ((c > 0) && (b >= first_blocks[c - 1]) && (b <= last_blocks[c - 1])) {
            there_is_left_block = true;
            l_p = Ps[(c - 1) * max_num_blocks + b];
            l_m = Ms[(c - 1) * max_num_blocks + b];
          } else {
            if (c == 0) {  // If there are no cells to the left (only boundary cells)
              there_is_left_block = true;
              left_score = b * kWordSize + block_pos + 1;
              upper_left_score = left_score - 1;
            } else {
              there_is_left_block = false;
            }
          }
          // Mark move
          alignment[alignment_length++] = move_code;
        } else {
          // Reached end - finished!
          break;
        }
      }
    }
  }

  // Edlib reverses the alignment here to get it in the correct order, but as I am aligning in
  // reverse order, I can just return the alignment as is when using prefix mode. Otherwise, do the reverse.
  if (!is_prefix) {
    std::reverse(alignment, alignment + alignment_length);
  }
}

/**
 * Finds one possible alignment that gives optimal score.
 * @param [in] align_reference
 * @param [in] aln_reference_length
 */
void Alignment::ObtainAlignment(const u8* align_reference, const s32 aln_reference_length) {
  const s32 max_num_blocks = CeilDiv(query_length, kWordSize);
  const s32 w = max_num_blocks * kWordSize - query_length;

  simd::BuildPeq(query, query_length, peq_table);
  MyersCalcEditDistanceNW(w, max_num_blocks, align_reference, aln_reference_length);
  return ObtainAlignmentTraceback(aln_reference_length);
}

PairwiseAlignment::PairwiseAlignment(FlowContext& exec, size_t batch_nr)
    : Task(exec, batch_nr, fmt::format("Align{}", batch_nr)) {}

void Alignment::ComputeConsensus(FixedReadRecord& read, DuplexMetrics& metrics) {
  is_prefix = alignment_mode == AlignmentMode::kPrefix;

  // prepare sequences for alignment
  if (!SetSequences(read)) {
    read.SetStatus(FixedReadRecord::Status::kTrimmedTooShortFail);
    UpdateMetrics(read, metrics);
    return;
  }

  // perform alignment
  if (!Align()) {
    read.SetStatus(FixedReadRecord::Status::kDuplexTooLongFail);
    UpdateMetrics(read, metrics);
    return;
  }

  // if edit_distance is -1, the alignment failed because it exceeded the maximum edit distance given error rate
  if (edit_distance >= 0) {
    // Successful alignment, this read will be written out if it passes the filters.

    // To get the fastest alignment we align the sequences in the reverse order (i.e. starting at the hairpin loop)
    // To write the consensus and quality sequences efficiently, we use this reverse representation of the alignment
    // to make the new results directly. To do we currently do so in 2 passes of this alignment representation.
    //
    // Main goal of function: create new consensus sequence and quality string as well as IUPAC-like tag
    //
    // Basics of 2 pass algorithm of processing alignment:
    // 1st pass:
    // - We create an untrimmed consensus sequence representation including the endadapter if it exists
    // - We try to find the start adapter from this sequence and keep track its position (usually we succeed)
    // - If the start adapter is not found we currently do not trim anything and still output the read
    // 2nd pass:
    // - Use the end position of the start adapter and trim it off to create the final consensus and quality string
    // - Once we do this we have to recalculate the final alignment.
    // - During this pass we also generate the IUPAC-like tag (i.e. YC tag)
    //
    // Note: For mismatches in consensus we replace the base with the base from the reference sequence (i.e. R1)
    // TODO: Instead of defaulting to reference base, consider using base with best quality base
    //
    // The output of the alignment is written to the FixedReadRecord object:
    // - The consensus sequence is built inside the dedicated consensus_seq field
    // - The new quality string is built the quality string Qual() field
    // - The when the consensus sequence is fully constructed the IUPAC-like tag will be copied into the Seq() field
    //
    // Storing the IUPAC-like tag string inside read.Seq() is a bit of a hack, but it allows us to reuse the space
    // as we don't need the original sequence anymore.
    // Notes on the IUPAC-like tag:
    // - It is based on the alignment (showing directionality of changes of R1 relative to R2) in the consensus
    //   allowing for reconstruction of the original sequences (minus adapters)

    // Assume the last position of the alignment is the start of the simplex region (in the reverse order).
    simplex_length = is_prefix ? (reference_length - end_locations[0] - 1) : start_locations[num_locations - 1];

    // Speed optimization: mark the base at the alignment length as an insert; this will force the loop track the
    // matching sections to stop, but that insert will be ignored as it is past the alignment length. And this
    // then will allow us to count the number of leading zero bytes in the alignment data, which can be implemented
    // efficiently using SIMD.
    // This mark is detected by the CopyMatchingBases() SIMD function, which will stop copying when it detects it
    alignment[alignment_length] = kEdopInsert;

    // Here we create the consensus sequence. We then look for the start adapter in the consensus sequence; if it
    // is present, we trim the consensus sequence. Note the simplex_length is not yet to the final value yet.
    if (!GenerateConsensusSeq(read)) {
      // read failed during consensus generation probably due to missing end adapter
      UpdateMetrics(read, metrics);
      return;
    }

    // Filter out reads less than min_length
    if (static_cast<size_t>(read.consensus_seq_len) < _params.min_trimmed_read_len) {
      read.SetStatus(FixedReadRecord::Status::kTrimmedTooShortFail);
      UpdateMetrics(read, metrics);
      // finish early
      return;
    }

    // perform strand detection
    if (_params.strand_classifier.has_value()) {
      // strand detection is enabled, so we need to check the strand
      _strand = _params.strand_classifier->Classify(std::string_view(read.ConsensusSeq(), read.consensus_seq_len));
      // All consensus alignments are in reverse direction so are naturally left aligned without intervention
      // We need to left align the forward strand reads because they are right aligned relative to the reference
      ProcessStrand(read);
    }

    // Here we set quality string and we create the IUPAC-like tag
    GenerateConsensusQualAndIUPAC(read);

    // Everything is finished and not filtered. What should be true here:
    // - Consensus sequence string is stored in consensus_seq (and offset is set to trim off bases)
    // - Consensus quality string is stored in read.Qual()
    // - IUPAC-like tag string is stored in read.Seq()
    // - read.status should remain as kDemultiplexed
  } else {
    // Alignment failed due to exceeding maximum edit distance
    read.SetStatus(FixedReadRecord::Status::kDuplexEditDistanceFail);
  }
  // update the metrics for this read
  UpdateMetrics(read, metrics);
}

// This function updates the duplex metrics for the read
// Uses the read status to determine what metric to update
void Alignment::UpdateMetrics(FixedReadRecord& read, DuplexMetrics& metrics) const {
  const auto& sid = read.trim_info_duplex.matches[DuplexMatch::kSID3p].match.barcode_id;
  switch (read.GetStatus()) {
    using enum FixedReadRecord::Status;
    case kDemultiplexed:
      // everything is fine, so we can continue
      break;
    case kDuplexEditDistanceFail:
      ++metrics.failed_assigned_counts.too_many_errors[sid];
      return;
    case kDuplexTooLongFail:
      ++metrics.failed_assigned_counts.consensus_too_long[sid];
      return;
    case kTrimmedTooShortFail:
      ++metrics.failed_assigned_counts.trimmed_read_too_short[sid];
      return;
    case kFailedMidadapterTrimFail:
      ++metrics.failed_assigned_counts.failed_hairpin_stem_trim_reads[sid];
      return;
    default:
      // This should never happen so log an error
      throw error::Error("Alignment::UpdateMetrics: unexpected read status {} for read {}",
                         static_cast<s32>(read.GetStatus()), read.Name());
  }

  // update non-failure but strange cases metrics
  metrics.passing_counts.longer_r2[sid] += is_swapped ? 1 : 0;

  const auto trim_size = static_cast<u32>(_endadapter_trim_pos);
  if (trim_size == 0) {
    ++metrics.passing_counts.no_endadapter[sid];
  }

  metrics.endadapter_position_distr[sid].AddCountToHistogram(trim_size, 1);

  // update base count metrics
  metrics.base_counts.concordant[sid] += _concordant_duplex_bases_total;
  metrics.base_counts.discordant[sid] += static_cast<u64>(_discordant_bases_total);
  metrics.base_counts.simplex[sid] += simplex_length;

  const auto report_length = static_cast<u32>(read.consensus_seq_len);

  // update the metrics on the read length distributions
  if (simplex_length == 0) {
    // full duplex read
    metrics.passing_counts.longer_r2_full_duplex[sid] += is_swapped ? 1 : 0;
    metrics.passing_counts.full_duplex[sid] += 1;
    metrics.full_duplex_length_distr[sid].AddCountToHistogram(report_length, 1);
  } else {
    // partial duplex read
    metrics.passing_counts.partial_duplex[sid] += 1;
    metrics.partial_duplex_length_distr[sid].AddCountToHistogram(report_length, 1);
  }
  metrics.passing_length_distr[sid].AddCountToHistogram(report_length, 1);

  // update strand detection counts
  switch (_strand) {
    using enum strand::StrandType;
    case kForward:
      ++metrics.strand_counts.fw[sid];
      break;
    case kForwardSig:
      ++metrics.strand_counts.fw_sig[sid];
      break;
    case kReverse:
      ++metrics.strand_counts.rv[sid];
      break;
    case kReverseSig:
      ++metrics.strand_counts.rv_sig[sid];
      break;
    default:
      break;
  }

  // update UMI metrics if applicable
  if (_adapter_type == AdapterType::kDuplexUMI) {
    if (read.trim_info_duplex.umi_3p.has_value()) {
      ++metrics.passing_counts.only_3p_umi[sid];
      if (read.trim_info_duplex.umi_5p.has_value()) {
        ++metrics.passing_counts.both_umi[sid];
      }
    }
    if (read.trim_info_duplex.umi_5p.has_value()) {
      ++metrics.passing_counts.only_5p_umi[sid];
    }
  }
}

/**
 * Helper function to process extra trimming and UMI determination if relevant.
 * @param read The read to process
 * @return True if the processing was successful, false if it failed (due to missing midadapter trim)
 */
bool Alignment::ProcessExtraTrim(FixedReadRecord& read) {
  // first is 5' and 3' side trimming positions
  auto [trim_start, trim_end] = _trimmer.FindUMIPos(read);

  // trim off the 3' side
  if (trim_end != -1) {
    read.consensus_seq_len = trim_end;

    if (read.consensus_seq_len < simplex_length) {
      // the simplex region cannot be longer than the consensus sequence so set it to the consensus sequence length
      // duplex region is completely trimmed off
      simplex_length = read.consensus_seq_len;
      // we don't have an alignment anymore so set the alignment length to 0
      alignment_length = 0;
    } else {
      // the simplex region is still shorter than the consensus sequence and not complete trimmed off
      // Alter alignment length to reflect the side it was trimmed
      // it is possible that the new alignment length is still shorter than the current alignment length because the
      // alignment isn't flush to the loop side of the read (these are technically indels)
      if (const auto new_alignment_length = read.consensus_seq_len - simplex_length;
          new_alignment_length < alignment_length) {
        alignment_length = new_alignment_length;
        // mark the end of the alignment for the SIMD function used afterward
        alignment[alignment_length] = kEdopInsert;
      }
    }
    // trim off the 5' side
    if (trim_start != -1) {
      _endadapter_trim_pos = trim_start;
    } else if (read.consensus_seq_len == 0) {
      // we can't trim something that doesn't have any bases
      read.SetStatus(FixedReadRecord::Status::kTrimmedTooShortFail);
      return false;
    } else {
      // we can't find so just trim via the end adapter position
      _endadapter_trim_pos = _trimmer.FindStartAdapterInConsensus(read);
    }
    return true;
  }
  // we can't find the extra loop adapter sequence so we don't know where the consensus sequence ends
  read.SetStatus(FixedReadRecord::Status::kFailedMidadapterTrimFail);
  return false;
}

/**
 * @brief Creates the consensus sequence and trims it to the start position.
 *
 * This function generates a consensus sequence from aligned reads and trims it to remove the start adapter.
 * It updates the `reference_location` and `query_location` arrays to reflect the new start position after trimming.
 * The resulting consensus sequence is stored in the `ConsensusSeq()` field of the `FixedReadRecord`,
 * with the `consensus_seq_offset` and `consensus_seq_len` fields adjusted accordingly.
 * @param [in,out] read The read record containing the sequences and alignment information.
 * @return True if the consensus sequence was successfully generated and trimmed; false otherwise.
 */
bool Alignment::GenerateConsensusSeq(FixedReadRecord& read) {
  // we have not cut anything off yet so assume we start in the beginning of the sequence
  s32 consensus_pos = 0;
  // start of duplex data
  s32 reference_pos = simplex_length;
  s32 query_pos = 0;

  // If read record is reused, we need to reset the consensus sequence length and offset
  auto* consensus_seq = read.consensus_buffer;

  // still in alphabet representation
  // original reference sequence IF NOT SWAPPED
  auto* reference_in{read.Seq() + insert_start1};

  // these are in alphabet (0..3) representations
  const auto* p_query{is_prefix ? reverse_aln_query : query};  // input data is now in correct output order,
  const auto* p_reference{is_prefix ? reverse_aln_reference : reference};

  // We first copy non-duplex bases (i.e. unaligned data) into consensus_seq, though likely later assign a shorter
  // consensus_seq_len if we find the endadapter (to effectively trim it off).
  if (is_swapped) {  // The S2 is longer than S1, should only occur with full reads
    for (s32 i = 0; i < simplex_length; ++i) {
      // copy via alphabet representation
      consensus_seq[i] = kAlphabet[p_reference[i]];
    }
  } else {  // fast path as we can directly copy from the reference (reference) sequence
    std::copy_n(reference_in, simplex_length, consensus_seq);
  }
  s32 output_size = simplex_length;

  // Now iterate over the alignment data; we expect that most of it is perfectly aligned (kEdopMatch).
  // So our main loop finds stretches of perfectly aligned data and writes them out efficiently.
  while (consensus_pos < alignment_length) {
    // Find the next stretch of perfectly aligned data. See simd-functions.cpp for more details.
    // This is a SIMD function that will only copy the matching bases from the alignment data.
    const s32 stretch_length =
        simd::CopyMatchingSeq(alignment + consensus_pos, reference_in + reference_pos, consensus_seq + output_size);
    if (stretch_length > 0) {
      // Copy out the perfectly aligned data.
      if (is_swapped) {  // slow path as we need to go via the alphabet
        for (s32 i = 0; i < stretch_length; ++i) {
          consensus_seq[output_size + i] = kAlphabet[p_query[query_pos + i]];
        }
      }
      query_pos += stretch_length;
      reference_pos += stretch_length;
      consensus_pos += stretch_length;
      output_size += stretch_length;
      // This wraps up the perfectly aligned stretch, which should be the majority of the alignment data.
    } else {  // stretch length was 0, which happens when we ran out of data or encountered a mismatch/indel.
      // TODO: consider masking only after first pass to preserve the original sequence for trimming or perform
      //               Start adaptor trimming before alignment
      switch (alignment[consensus_pos]) {
        case kEdopInsert:
          // Insertion into reference = deletion from query.
          consensus_seq[output_size] = mask_discordant_bases ? 'N' : kAlphabet[p_query[query_pos]];
          ++query_pos;
          break;
        case kEdopDelete:
          // Deletion from reference = insertion to query.
          consensus_seq[output_size] = mask_discordant_bases ? 'N' : kAlphabet[p_reference[reference_pos]];
          ++reference_pos;
          break;
        case kEdopMismatch: {
          // Mismatch = substitution. Note that we might have to swap the query and reference bases.
          const s32 true_ref_base = is_swapped ? p_query[query_pos] : p_reference[reference_pos];
          consensus_seq[output_size] = mask_discordant_bases ? 'N' : kAlphabet[true_ref_base];
          ++query_pos;
          ++reference_pos;
          break;
        }
        default:
          throw error::Error(
              "Alignment::GenerateConsensusSeq: unexpected alignment code {} at consensus position {} for read {}",
              static_cast<s32>(alignment[consensus_pos]), consensus_pos, std::string_view(read.Name(), read.NameLen()));
      }
      ++output_size;
      ++consensus_pos;
    }
  }
  // finish off the consensus sequence if it has remaining bases
  if (reference_pos < reference_length) {
    if (is_swapped) {
      // copy via alphabet representation
      for (; reference_pos < reference_length; ++reference_pos) {
        consensus_seq[output_size++] = kAlphabet[p_reference[reference_pos]];
      }
    } else {
      // fast path as we can directly copy from the reference sequence
      const auto remaining_bases = reference_length - reference_pos;
      std::copy_n(reference_in + reference_pos, remaining_bases, consensus_seq + output_size);
      output_size += remaining_bases;
    }
  }

  read.consensus_seq_len = output_size;

  // End of the alignment, now we trim the consensus sequence if needed.
  switch (_adapter_type) {
    using enum AdapterType;
    case kDuplexStem:
    case kDuplexUMI:
      if (!ProcessExtraTrim(read)) {
        // read failed during extra trimming, due to missing extra midadapter trim location
        return false;
      }
      break;
    case kDuplex:
      // regular trimming
      _endadapter_trim_pos = _trimmer.FindStartAdapterInConsensus(read);
      break;
    default:
      throw error::Error("Alignment::GenerateConsensusSeq: unexpected duplex type {}", static_cast<s32>(_adapter_type));
  }
  // correct for off by one to get the position after the adapter
  _endadapter_trim_pos += 1;
  // This should never happen but if the start position is beyond the consensus sequence length error out
  if (_endadapter_trim_pos == read.consensus_seq_len) {
    // read will be filtered out as it is empty
    read.consensus_seq_len = 0;
  } else if (_endadapter_trim_pos > read.consensus_seq_len) {
    throw error::Error(
        "Alignment::GenerateConsensusSeq: start position {} is beyond consensus sequence length {} for read {}",
        _endadapter_trim_pos, read.consensus_seq_len, std::string_view(read.Name(), read.NameLen()));
  } else {
    // Trim the consensus sequence to the start position
    read.consensus_seq_offset = static_cast<u32>(_endadapter_trim_pos);
    read.consensus_seq_len -= _endadapter_trim_pos;
  }
  return true;
}

/**
 * @brief This function processes the strand of the read, shifting the alignment data appropriately
 * @param read The read to process
 * swapped
 **/
void Alignment::ProcessStrand(const FixedReadRecord& read) {
  switch (_strand) {
    using enum strand::StrandType;
    case kForward:
    case kForwardSig: {  // forward relative to the reference
      // isolate only the non-simplex region (untrimmed consensus region, indexes shared with alignment)
      const auto consensus_only = std::string_view(read.consensus_buffer + simplex_length, alignment_length);
      // move the read indels to the left
      LeftAlignConvert(consensus_only, alignment);
      break;
    }
    case kReverse:
    case kReverseSig: {  // is already left aligned
      const auto consensus_only = std::string_view(read.consensus_buffer + simplex_length, alignment_length);
      GroupRightAlignConvert(consensus_only, alignment);
      break;
    }
    default:
      // do nothing all things being equal
      break;
  }
}

/**
 * @brief Generates the consensus quality string and the IUPAC-like (YC) tag.
 *
 * The consensus quality string is stored in the `Qual()` field of the `FixedReadRecord`,
 * and its size is set to `consensus_seq_len`. The IUPAC-like tag is stored in the `Seq()` field
 * of the `FixedReadRecord`, and its length is stored in the `iupac_length` field.
 */
void Alignment::GenerateConsensusQualAndIUPAC(FixedReadRecord& read) {
  _concordant_duplex_bases_total = 0;
  _discordant_bases_total = 0;

  s32 consensus_pos = 0;
  s32 query_pos = 0;
  // start of duplex data
  s32 reference_pos = simplex_length;
  s32 output_size = 0;

  auto* consensus_qual = read.Qual();

  // these are in alphabet (0..3) representations
  const auto* p_query{is_prefix ? reverse_aln_query : query};  // input data is now in correct output order,
  const auto* p_reference{is_prefix ? reverse_aln_reference : reference};

  // determine where the reference and query start if the read is full duplex
  // If the simplex region length is less than the endadapter trim position, we need to adjust
  // the reference and query start positions. This occurs when the read has a simplex region
  // (non-duplex region) that overlaps with the trimmed duplex region. The alignment must be
  // traversed to account for insertions and deletions, ensuring the correct start position.
  if (simplex_length < _endadapter_trim_pos) {
    const s32 loop_end_condition = _endadapter_trim_pos - simplex_length;
    for (; consensus_pos < loop_end_condition; ++consensus_pos) {
      // Process each alignment operation to adjust query and reference positions.
      switch (alignment[consensus_pos]) {
        case kEdopInsert:
          // Insertion in the query sequence; move the query position forward.
          ++query_pos;
          break;
        case kEdopDelete:
          // Deletion in the reference sequence; move the reference position forward.
          ++reference_pos;
          break;
        case kEdopMatch:
        case kEdopMismatch:
          // Match or mismatch; move both query and reference positions forward.
          ++reference_pos;
          ++query_pos;
          break;
        default:
          // Unexpected alignment code; throw an error with details for debugging.
          throw error::Error(
              "Alignment::GenerateConsensusQualAndIUPAC: unexpected alignment code {} at position {} for read {}",
              static_cast<s32>(alignment[consensus_pos]), consensus_pos, std::string_view(read.Name(), read.NameLen()));
      }
    }
    simplex_length = 0;
  } else {
    // read is not full length, so our start positions are still equal to start_locations[num_locations - 1]
    simplex_length -= _endadapter_trim_pos;
    // fill simplex region with base quality
    std::fill_n(consensus_qual, simplex_length, kSimplexBaseQual);
    output_size = simplex_length;
  }

  u32 iupac_length = 0;

  // Methylation mode
  if (_params.methylation) {
    iupac_cigar[iupac_length++] = ' ';
    // initialize IUPAC cigar string with tag prefix XM:Z:
    iupac_cigar[iupac_length++] = 'X';
    iupac_cigar[iupac_length++] = 'M';
    iupac_cigar[iupac_length++] = ':';
    iupac_cigar[iupac_length++] = 'Z';
    iupac_cigar[iupac_length++] = ':';

    // convert all bases of xm_tag to '.' stored in iupac_cigar
    std::fill_n(iupac_cigar + iupac_length, read.consensus_seq_len, '.');
    // get R2
    auto r2 = is_swapped ? p_reference + reference_pos : p_query + query_pos;
    MethylConvert(read.ConsensusSeq() + simplex_length, alignment + consensus_pos, alignment_length - consensus_pos,
                  iupac_cigar + iupac_length + simplex_length, r2, is_swapped);
    iupac_length += read.consensus_seq_len;
  }
  iupac_cigar[iupac_length++] = ' ';
  // create first part of YC tag
  // initialize IUPAC cigar string with tag prefix YC:Z:
  iupac_cigar[iupac_length++] = 'Y';
  iupac_cigar[iupac_length++] = 'C';
  iupac_cigar[iupac_length++] = ':';
  iupac_cigar[iupac_length++] = 'Z';
  iupac_cigar[iupac_length++] = ':';
  // now add simplex region in decimal
  if (simplex_length > 0) {
    s32 value = simplex_length;
    // buffer to hold the digits in reverse order; max 4 digits for simplex_length <= kMaxAlignLength (1504)
    char buffer[4];
    s32 len = 0;
    do {
      buffer[len++] = static_cast<char>('0' + (value % 10));
      value /= 10;
    } while (value > 0);
    // Append digits in correct order
    for (s32 i = len - 1; i >= 0; --i) {
      iupac_cigar[iupac_length++] = buffer[i];
    }
    iupac_cigar[iupac_length++] = is_swapped ? kSwappedJunction : kDefaultJunction;
  } else {
    iupac_cigar[iupac_length++] = kDefaultJunction;
  }

  // Now iterate over the alignment data; we expect that most of it is perfectly aligned (kEdopMatch).
  // So our main loop finds stretches of perfectly aligned data and writes them out efficiently.
  // The consensus position, reference position, and query position will likely be reassigned to a different value
  // on this second pass due to the start adapter being found (effectively cut off).
  while (consensus_pos < alignment_length) {
    // Find the next stretch of perfectly aligned data. See simd-functions.cpp for more details.
    s32 stretch_length =
        simd::CopyMatchingQual(alignment + consensus_pos, consensus_qual + output_size, kConcordBaseQual);
    if (stretch_length > 0) {
      query_pos += stretch_length;
      reference_pos += stretch_length;
      consensus_pos += stretch_length;
      _concordant_duplex_bases_total += static_cast<u32>(stretch_length);
      output_size += stretch_length;

      // The IUPAC code for a stretch of perfectly aligned data is the number of bases, i.e. the
      // stretch length. Instead of using itos(), itoa() or sprintf(), we will write out the digits
      // ourselves; this code works for values [1..999] and is very fast as I don't have to deal with
      // unnecessary error checking or memory allocation for strings.
      // TODO: Potentially unsafe so consider using a safer version of this code and measuring performance
      constexpr s32 kDecimalBase{10};
      s32 unit = 100;
      while (stretch_length < unit) {
        unit /= kDecimalBase;
      }
      while (unit) {
        const s32 digit = stretch_length / unit;
        iupac_cigar[iupac_length++] = static_cast<char>('0' + digit);
        stretch_length -= digit * unit;
        unit /= kDecimalBase;
      }
      // This wraps up the perfectly aligned stretch, which should be the majority of the alignment data.
    } else {  // stretch length was 0, which happens when we ran out of data or encountered a mismatch/indel.
      switch (alignment[consensus_pos]) {
        case kEdopInsert:
          // Insertion into reference = deletion from query.
          iupac_cigar[iupac_length++] =
              is_swapped ? kIUPACTable[p_query[query_pos]][kIndelIndex] : kIUPACTable[kIndelIndex][p_query[query_pos]];
          ++query_pos;
          consensus_qual[output_size++] = kDiscordBaseQual;
          ++_discordant_bases_total;
          break;
        case kEdopDelete:
          // Deletion from reference = insertion to query.
          iupac_cigar[iupac_length++] = is_swapped ? kIUPACTable[kIndelIndex][p_reference[reference_pos]]
                                                   : kIUPACTable[p_reference[reference_pos]][kIndelIndex];
          ++reference_pos;
          consensus_qual[output_size++] = kDiscordBaseQual;
          ++_discordant_bases_total;
          break;
        case kEdopMismatch: {
          // Mismatch = substitution. Note that we might have to swap the query and reference bases.
          const s32 ref_base = is_swapped ? p_query[query_pos] : p_reference[reference_pos];
          const s32 query_base = is_swapped ? p_reference[reference_pos] : p_query[query_pos];
          iupac_cigar[iupac_length++] = kIUPACTable[ref_base][query_base];
          ++query_pos;
          ++reference_pos;
          consensus_qual[output_size++] = kDiscordBaseQual;
          ++_discordant_bases_total;
          break;
        }
        case kEdopMatchMethyl: {
          // Only used in methylation mode
          // Treat as mismatch for IUPAC but concordant for quality
          const s32 ref_base = is_swapped ? p_query[query_pos] : p_reference[reference_pos];
          const s32 query_base = is_swapped ? p_reference[reference_pos] : p_query[query_pos];
          iupac_cigar[iupac_length++] = kIUPACTable[ref_base][query_base];
          ++query_pos;
          ++reference_pos;
          consensus_qual[output_size++] = kConcordBaseQual;
          ++_concordant_duplex_bases_total;
          break;
        }
        default:
          throw error::Error("Alignment::GenerateConsensusQualAndIUPAC: unexpected alignment code {} for read {}",
                             static_cast<s32>(alignment[consensus_pos]), std::string_view(read.Name(), read.NameLen()));
      }
      ++consensus_pos;
    }
  }
  // finish off the consensus quality string if past end of alignment but only based on consensus read length
  if (output_size < read.consensus_seq_len) {
    const auto remaining_bases = read.consensus_seq_len - output_size;
    std::fill_n(consensus_qual + output_size, remaining_bases, kDiscordBaseQual);
    _discordant_bases_total += remaining_bases;

    for (s32 i = 0; i < remaining_bases; ++reference_pos, ++i) {
      iupac_cigar[iupac_length++] = is_swapped ? kIUPACTable[kIndelIndex][p_reference[reference_pos]]
                                               : kIUPACTable[p_reference[reference_pos]][kIndelIndex];
    }
  }

  // Check iupac_length doesn't exceed the maximum allowed length
  if (iupac_length >= kMaxAlignLength + 2) {
    // +2 for the -0 or +0
    Logging::Warn(
        "Alignment::GenerateConsensusQualAndIUPAC: Comment length {} exceeds maximum allowed length {} for "
        "read {}",
        iupac_length, kMaxAlignLength, std::string_view(read.Name(), read.NameLen()));
    // Set too long failure status
    read.SetStatus(FixedReadRecord::Status::kDuplexTooLongFail);
  }

  // Stored in the read.Seq() field, as we store the consensus sequence elsewhere so we can us to reuse the space.
  const auto* iter = fmt::format_to(read.Seq(), "{}+", std::string_view(iupac_cigar, iupac_cigar + iupac_length));
  read.iupac_length = static_cast<u32>(iter - read.Seq());
}

void PairwiseAlignment::operator()() const {
  try {
    const StopWatch sw;

    const auto& batch{context.GetBatchData(batch_nr)};
    const auto& param{context.DemuxParam()};
    const auto error_rate{param.max_error_rate_percent};

    // Loop over all records in the batch and perform alignment.
    if (batch.records) {
      // Create an aligner object; this object uses a lot of stack space for intermediate alignment results, and
      // allocating space on the stack actually takes quite some time (I benchmarked it to be about 5% of the total
      // time). I obtained a significant speedup by creating the aligner object for every batch and re-using it rather
      // than assuming that stack allocation is "free". Use a max error rate of 10%, do not mask discordant bases.

      Alignment aligner(error_rate, false, AlignmentMode::kInfix, context.GetManager().DemuxObjectDuplex(),
                        context.DemuxParam(), context.GetManager().adapter_type);

      auto& duplex_metrics{DuplexMetrics::Instance()};
      for (size_t i = 0; i < batch.Size(); ++i) {
        if (auto& read = (*(batch.records))[i]; read.GetStatus() == FixedReadRecord::Status::kDemultiplexed) {
          aligner.ComputeConsensus(read, duplex_metrics);
        }
      }
    }

    context.AddToAlignmentTime(sw.ElapsedTime());
  } catch (const std::exception& e) {
    Logging::Error("PairwiseAlignment::operator() failed: {}", e.what());
    SetTaskException(std::current_exception());
  }
}

// NOTE: If you change the parameters here check SIMDNA
Alignment::Alignment(const float max_error, const bool mask_bases, const AlignmentMode mode,  // NOLINT
                     const DemuxAndTrimDuplex& trimmer, const DemuxAndTrimParam& params,      // NOLINT
                     const AdapterType adapter_type)                                          // NOLINT
    : max_error_rate(max_error),
      mask_discordant_bases(mask_bases),
      alignment_mode{mode},
      _trimmer(trimmer),
      _params(params),
      _adapter_type(adapter_type) {}
}  // namespace xoos::demux
