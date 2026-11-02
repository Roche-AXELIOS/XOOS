#include "task/alignment.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <span>
#include <utility>

#include <fmt/format.h>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>

#include "metrics/duplex-metrics.h"
#include "task/flow-context.h"
#include "task/flow-manager.h"
#include "utility/alignment-util.h"
#include "utility/math-util.h"
#include "utility/stop-watch.h"
#include "xoos/util/sequence-functions.h"

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
 * We support both prefix and infix alignment modes. Infix (HW) penalizes gaps at the beginning of the query
 * sequence, but not the end. Prefix (SHW) does not penalize either end. Prefix alignment is approximately
 * 60% faster than infix alignment. Prefix alignment
 * mode anchors the alignment at the loop side of the reversed reference, so the gap at the query end (the simplex
 * region) is not penalized. This is appropriate for all duplex reads — both partial and full — because the simplex
 * region should not contribute to the alignment score. We default to prefix alignment.
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
  if (static_cast<u32>(reference_length) > kMaxAlignLength) {
    return false;
  }

  // reset variables
  edit_distance = -1;
  num_locations = 0;
  alignment_length = 0;
  _had_multi_candidates = false;
  _alignment = _aln_buffer.data();

  /*--------------------- INITIALIZATION ------------------*/
  // bmax in Myers
  const s32 max_num_blocks = CeilDiv(query_length, kWordSize);
  // number of redundant cells in last level blocks
  const s32 num_redundant_cells = (max_num_blocks * kWordSize) - query_length;

  simd::BuildPeq(query, query_length, peq_table);

  /*------------------ MAIN CALCULATION -------------------*/
  // only allow a specific error rate for the alignment. This means that we can calculate the maximum
  // edit distance that we'd allow - and smaller edit distances will be calculated faster.
  const s32 k =
      1 + std::max(0, static_cast<s32>(max_error_rate * static_cast<f32>(std::min(query_length, reference_length))));

  MyersCalcEditDistanceSemiGlobal(peq_table,
                                  num_redundant_cells,
                                  max_num_blocks,
                                  k,
                                  is_prefix,
                                  reference,
                                  reference_length,
                                  edit_distance,
                                  end_locations,
                                  num_locations);

  if (edit_distance >= 0) {
    // If there is solution.
    if (is_prefix) {
      return SelectBestPrefixAlign();
    }

    // Infix: penalize edit distance by unaligned reference tail
    const auto reference_tail_length = reference_length - 1 - end_locations[num_locations - 1];
    if ((edit_distance + reference_tail_length) > k) {
      edit_distance = -1;
      return true;
    }

    // Compute start locations via reverse alignment
    simd::ReverseCopy(reference, reference + reference_length, _infix->reverse_reference.data());
    simd::ReverseCopy(query, query + query_length, _infix->reverse_query.data());
    simd::BuildPeq(_infix->reverse_query.data(), query_length, _infix->reverse_peq_table.data());
    for (u32 i = 0; i < num_locations; ++i) {
      const s32 end_location = end_locations[i];
      if (end_location == -1) {
        _infix->start_locations[i] = 0;
      } else {
        // score is not used here
        s32 score = 0;
        u32 num_reverse_positions = 0;
        MyersCalcEditDistanceSemiGlobal(_infix->reverse_peq_table.data(),
                                        num_redundant_cells,
                                        max_num_blocks,
                                        k,
                                        true,
                                        _infix->reverse_reference.data() + reference_length - end_location - 1,
                                        end_location + 1,
                                        score,
                                        _infix->reverse_positions.data(),
                                        num_reverse_positions);
        _infix->start_locations[i] = end_location - _infix->reverse_positions[num_reverse_positions - 1];
      }
    }

    return SelectBestInfixAlign();
  }
  return true;
}

bool Alignment::SelectBestPrefixAlign() {
  simd::ReverseCopy(reference, reference + reference_length, reverse_aln_reference);
  simd::ReverseCopy(query, query + query_length, reverse_aln_query);

  // assume at first narrowest alignment is best
  u32 best_index = 0;

  if (num_locations > 1) {
    _had_multi_candidates = true;

    // Order candidates by |aln_ref_length - query_length| ascending; tie-break preferring the
    // narrower (smaller end_location) candidate. The intuition is that the candidate whose
    // reference span is closest to the query length tends to need fewer indels, so evaluating
    // those first improves the odds that the zero-gaps early exit fires before we walk every
    // tied candidate.
    auto& order = _prefix->order;
    for (u32 i = 0; i < num_locations; ++i) {
      order[i] = i;
    }
    std::sort(order.begin(), order.begin() + num_locations, [this](const u32 a, const u32 b) {
      const s32 diff_a = std::abs((end_locations[a] + 1) - query_length);
      const s32 diff_b = std::abs((end_locations[b] + 1) - query_length);
      if (diff_a != diff_b) {
        return diff_a < diff_b;
      }
      return end_locations[a] < end_locations[b];
    });

    // Shared NW fill at the widest end_location with band tightening disabled, so the saved DP covers
    // every cell any candidate's traceback might need. Per-candidate work is then just a cheap traceback.
    // Note: end_locations is sorted ascending by the SemiGlobal search, so the widest is at index num_locations-1.
    const s32 max_num_blocks = CeilDiv(query_length, kWordSize);
    const s32 w = (max_num_blocks * kWordSize) - query_length;
    const s32 widest_aln_ref_length = end_locations[num_locations - 1] + 1;
    MyersCalcEditDistanceNW(w,
                            max_num_blocks,
                            std::span<const u8>(reference, static_cast<size_t>(widest_aln_ref_length)),
                            /*disable_band_tightening=*/true);

    // Evaluate the first candidate in sort order — traceback into the primary buffer.
    best_index = order[0];
    ObtainAlignmentTraceback(end_locations[best_index] + 1, _alignment, alignment_length);

    u8* alt_buf = _alt_aln_buffer.data();

    // Count number of indel ops and pick the one with the fewest
    const auto count_indels = [](const u8* const aln, const s32 len) {
      return std::count_if(aln, aln + len, [](const u8 op) { return op == kEdopInsert || op == kEdopDelete; });
    };

    auto best_gaps = count_indels(_alignment, alignment_length);

    for (u32 rank = 1; rank < num_locations; ++rank) {
      const u32 index = order[rank];

      // Traceback only — reuses the shared DP above.
      const s32 aln_ref_length = end_locations[index] + 1;
      auto cur_aln_len = 0;
      ObtainAlignmentTraceback(aln_ref_length, alt_buf, cur_aln_len);

      const auto gaps = count_indels(alt_buf, cur_aln_len);

      if (gaps < best_gaps) {
        best_gaps = gaps;
        best_index = index;
        // Keep the selected alignment in the current buffer and reuse the other buffer
        // for the next candidate. SmearAlignment selects that other buffer dynamically
        // when it needs scratch space, so pointer swapping cannot alias the traceback.
        std::swap(_alignment, alt_buf);
        alignment_length = cur_aln_len;
        // Zero indels is optimal — no need to evaluate further candidates.
        if (best_gaps == 0) {
          break;
        }
      }
    }
  } else {
    // Single candidate: original code path (banded NW fill + traceback).
    ObtainAlignment(std::span<const u8>(reference, static_cast<size_t>(end_locations[best_index] + 1)),
                    _alignment,
                    alignment_length);
  }

  // Compute simplex_length from the winning candidate
  simplex_length = reference_length - end_locations[best_index] - 1;

  const auto consensus_length = static_cast<u32>(simplex_length + alignment_length);
  return (consensus_length <= kMaxAlignLength);
}

bool Alignment::SelectBestInfixAlign() {
  // Pick the alignment closest to the end of the reference
  const s32 aln_start = _infix->start_locations[num_locations - 1];
  const s32 aln_end = end_locations[num_locations - 1];
  const s32 aln_reference_length = aln_end - aln_start + 1;

  simd::ReverseCopy(reference, reference + reference_length, reverse_aln_reference);
  simd::ReverseCopy(query, query + query_length, reverse_aln_query);
  ObtainAlignment(std::span<const u8>(reference + aln_start, static_cast<size_t>(aln_reference_length)),
                  _alignment,
                  alignment_length);

  // Compute simplex_length
  simplex_length = _infix->start_locations[num_locations - 1];

  // Bound the eventual consensus length against kMaxAlignLength
  const auto reference_tail_length = reference_length - 1 - end_locations[num_locations - 1];
  const auto consensus_length = static_cast<u32>(simplex_length + alignment_length + reference_tail_length);
  return (consensus_length <= kMaxAlignLength);
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
inline Delta Alignment::CalculateBlock(const Word pos_v,
                                       const Word minus_v,
                                       Word eq,
                                       const Word pos_hin,
                                       const Word minus_hin,
                                       Word& pos_hout,
                                       Word& minus_hout) {
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

  // equivalent to ph <<= 1 / mh <<= 1, but usually faster
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

static inline s32 Min(const s32 x, const s32 y) {
  return x < y ? x : y;
}

static inline s32 Max(const s32 x, const s32 y) {
  return x > y ? x : y;
}

/**
 * @param [in] block
 * @param [in/out] scores array with scores
 * @return Values of cells in block, starting with bottom cell in block.
 */
void Alignment::GetBlockCellValues(const Block& block, const std::span<s32, static_cast<size_t>(kWordSize)> scores) {
  s32 score = block.score;
  Word mask = kHighBitMask;
  for (size_t i = 0; i < static_cast<size_t>(kWordSize) - 1; ++i) {
    scores[i] = score;
    if ((block.p & mask) != 0U) {
      --score;
    }
    if ((block.m & mask) != 0U) {
      ++score;
    }
    mask >>= 1;
  }
  scores[static_cast<size_t>(kWordSize) - 1] = score;
}

/**
 * Uses Myers' bit-vector algorithm to find edit distance for one of semi-global alignment methods.
 * @param [in] w  Size of padding in last block.
 * @param [in] max_num_blocks  Number of blocks needed to cover the whole query.
 * @param [in] k
 */
void Alignment::MyersCalcEditDistanceSemiGlobal(const Word* const peq_input,
                                                const s32 w,
                                                const s32 max_num_blocks,
                                                s32 k,
                                                const bool do_prefix,
                                                const u8* const target,
                                                const s32 target_length,
                                                s32& score,
                                                s32* const positions_out,
                                                u32& num_positions_out) {
  num_positions_out = 0;

  // first_block is 0-based index of first block in Ukkonen band.
  // last_block is 0-based index of last block in Ukkonen band.
  s32 first_block = 0;
  // y in Myers
  s32 last_block = Min(CeilDiv(k + 1, kWordSize), max_num_blocks) - 1;
  // Empty band: no block to compute means no alignment can be found.
  if (last_block < 0) {
    score = -1;
    return;
  }
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
  const u8* reference_char = target;
  for (s32 c = 0; c < target_length; ++c) {  // for each column
    const Word* peq_c = peq_input + ((*reference_char) * max_num_blocks);

    //----------------------- Calculate column -------------------------//
    auto pos_hout = static_cast<s32>(pos_start_hout);
    s32 minus_hout = 0;
    s32 hout = 0;
    for (s32 b = first_block; b <= last_block; ++b) {
      auto& block{blocks[b]};
      const auto delta_hout = CalculateBlock(
          block.p, block.m, peq_c[b], static_cast<Word>(pos_hout), static_cast<Word>(minus_hout), block.p, block.m);
      pos_hout = delta_hout.pos_hout;
      minus_hout = delta_hout.minus_hout;
      hout = pos_hout - minus_hout;
      block.score += hout;
    }
    Block* bl = blocks + last_block;
    peq_c += last_block;

    //---------- Adjust number of blocks according to Ukkonen ----------//
    // bl is pointing to last block
    if ((last_block < max_num_blocks - 1) && (bl->score - hout <= k) &&
        (((*(peq_c + 1) & kWord1) != 0U) || (hout < 0))) {
      // peq_c is pointing to last block
      // If score of left block is not too big, calculate one more block
      ++last_block;
      ++bl;
      ++peq_c;
      // All 1s
      bl->p = static_cast<Word>(-1);
      bl->m = static_cast<Word>(0);
      const auto delta_hout = CalculateBlock(
          bl->p, bl->m, *peq_c, static_cast<Word>(pos_hout), static_cast<Word>(minus_hout), bl->p, bl->m);
      bl->score = (bl - 1)->score - hout + kWordSize + delta_hout.pos_hout - delta_hout.minus_hout;
    } else {
      // Drive the shrink by index so we never form a pointer before 'blocks' (M23_080); 'bl' is
      // resynced from last_block below (peq_c is not read again this iteration).
      while ((last_block >= first_block) && (blocks[last_block].score >= (k + kWordSize))) {
        --last_block;
      }
    }

    // For infix mode, even if all cells are > k, there still may be a solution in the next column
    // because starting conditions at upper boundary are 0. This means that first block is always a
    // candidate for solution, and we can never end calculation before the last column.
    if (!do_prefix && (last_block == -1)) {
      ++last_block;
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
        num_positions_out = static_cast<u32>(nr_positions);
        std::copy_n(positions, nr_positions, positions_out);
      }
      return;
    }

    // Resync 'bl' only once last_block >= first_block >= 0 is guaranteed, so we never form a
    // pointer before 'blocks' (M23_080).
    bl = blocks + last_block;

    //------------------------- Update best score ----------------------//
    if (last_block == max_num_blocks - 1) {
      const s32 col_score = bl->score;
      // Scores > k don't have correct values (so we cannot use them), but are certainly > k.
      // NOTE: Score that I find in column c is actually score from column c-w.
      if ((col_score <= k) && ((best_score == -1) || (col_score <= best_score))) {
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
    ++reference_char;
  }

  // Obtain results for last w columns from last column.
  if (last_block == max_num_blocks - 1) {
    std::array<s32, static_cast<size_t>(kWordSize)> block_scores{};
    GetBlockCellValues(blocks[last_block], block_scores);
    for (s32 i = 0; i < w; ++i) {
      if (const s32 col_score = block_scores[i + 1];
          (col_score <= k) && (best_score == -1 || (col_score <= best_score))) {
        if (col_score != best_score) {
          nr_positions = 0;
          best_score = col_score;
          k = best_score;
        }
        positions[nr_positions++] = target_length - w + i;
      }
    }
  }

  score = best_score;
  if (best_score != -1) {
    num_positions_out = static_cast<u32>(nr_positions);
    std::copy_n(positions, nr_positions, positions_out);
  }
}

void Alignment::MyersCalcEditDistanceNW(const s32 w,
                                        const s32 max_num_blocks,
                                        const std::span<const u8> align_reference,
                                        const bool disable_band_tightening) {
  const auto align_reference_length = static_cast<s32>(align_reference.size());
  s32 k = edit_distance;
  if (!disable_band_tightening && k < ::abs(align_reference_length - query_length)) {
    return;
  }

  // Upper bound for k
  k = Min(k, Max(query_length, align_reference_length));

  // first_block is 0-based index of first block in Ukkonen band.
  // last_block is 0-based index of last block in Ukkonen band.
  s32 first_block = 0;
  // This is optimal now, by my formula.
  // When band tightening is disabled, force full band height so the saved DP covers every cell any
  // traceback (potentially at a narrower end column) might walk back through.
  s32 last_block =
      disable_band_tightening
          ? max_num_blocks - 1
          : Min(max_num_blocks, CeilDiv(Min(k, (k + query_length - align_reference_length) / 2) + 1, kWordSize)) - 1;
  // Empty band: no block to compute means no alignment can be found.
  if (last_block < 0) {
    return;
  }

  Block blocks[kMaxNrBlocks];

  // Initialize P, M and score
  for (s32 b = 0; b <= last_block; ++b) {
    blocks[b].score = (b + 1) * kWordSize;
    // All 1s
    blocks[b].p = static_cast<Word>(-1);
    blocks[b].m = static_cast<Word>(0);
  }

  // Max diagonal span reachable across two adjacent Myers blocks, used throughout the Ukkonen band bounds below.
  const s32 two_block_diagonal_span = 2 * (kWordSize - 1);

  const u8* reference_char = align_reference.data();
  for (s32 c = 0; c < align_reference_length; ++c) {  // for each column
    const Word* const peq_c = peq_table + (*reference_char * max_num_blocks);

    //----------------------- Calculate column -------------------------//
    s32 hout = 1;
    s32 pos_hout = 1;
    s32 minus_hout = 0;
    // Unroll loop to process two blocks at a time
    for (s32 b = first_block; b <= last_block; b += 2) {
      auto& block1 = blocks[b];
      const auto delta_hout1 = CalculateBlock(
          block1.p, block1.m, peq_c[b], static_cast<Word>(pos_hout), static_cast<Word>(minus_hout), block1.p, block1.m);
      pos_hout = delta_hout1.pos_hout;
      minus_hout = delta_hout1.minus_hout;
      hout = pos_hout - minus_hout;
      block1.score += hout;

      if (b + 1 <= last_block) {
        auto& block2 = blocks[b + 1];
        const auto delta_hout2 = CalculateBlock(block2.p,
                                                block2.m,
                                                peq_c[b + 1],
                                                static_cast<Word>(pos_hout),
                                                static_cast<Word>(minus_hout),
                                                block2.p,
                                                block2.m);
        pos_hout = delta_hout2.pos_hout;
        minus_hout = delta_hout2.minus_hout;
        // NOSONAR(S1854) -- carries the last computed block's horizontal delta into the k-update and
        // block-shrink logic below; the edlib kernel relies on this value across the column iteration.
        hout = pos_hout - minus_hout;  // NOSONAR
        block2.score += hout;
      }
    }

    auto* bl = &blocks[last_block];
    //------------------------------------------------------------------//
    // bl now points to last block

    // Update k. I do it only on end of column because it would slow calculation too much otherwise.
    // NOTICE: I add w when in last block because it is actually result from w cells to the left and w cells up.
    k = Min(k,
            bl->score + Max(align_reference_length - c - 1, query_length - ((1 + last_block) * kWordSize - 1) - 1) +
                (last_block == max_num_blocks - 1 ? w : 0));

    //---------- Adjust number of blocks according to Ukkonen ----------//
    //--- Adjust last block ---//
    // Skipped entirely when band tightening is disabled: we keep first_block=0 and last_block=max_num_blocks-1
    // for every column so any later traceback can read any cell.
    if (!disable_band_tightening) {
      // If block is not beneath band, calculate next block. Only next because others are certainly beneath band.
      if (((last_block + 1) < max_num_blocks) &&
          !(((last_block + 1) * kWordSize - 1) >
            (k - bl->score + two_block_diagonal_span - align_reference_length + c + query_length))) {
        ++last_block;
        ++bl;
        // All 1s
        bl->p = static_cast<Word>(-1);
        bl->m = static_cast<Word>(0);
        const auto delta_hout = CalculateBlock(
            bl->p, bl->m, peq_c[last_block], static_cast<Word>(pos_hout), static_cast<Word>(minus_hout), bl->p, bl->m);
        pos_hout = delta_hout.pos_hout;
        minus_hout = delta_hout.minus_hout;
        const s32 new_hout = pos_hout - minus_hout;
        bl->score = (bl - 1)->score - hout + kWordSize + new_hout;
        hout = new_hout;
      }

      // While block is out of band, move one block up.
      // NOTE: Condition used here is more loose than the one from the article, since I simplified the Max() part of it.
      // I could consider adding that Max part, for optimal performance.
      // Drive the shrink by index so we never form a pointer before 'blocks' (M23_080). 'bl' is not
      // read after this loop; the block is addressed via blocks[last_block] below.
      while ((last_block >= first_block) &&
             ((blocks[last_block].score >= (k + kWordSize)) ||
              ((((last_block + 1) * kWordSize) - 1) >
               ((k - blocks[last_block].score + two_block_diagonal_span - align_reference_length + c + query_length) +
                1)))) {
        --last_block;
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
    }
    //------------------------------------------------------------------//

    //---- Save column so it can be used for reconstruction ----//
    if (c < align_reference_length) {
      for (s32 b = first_block; b <= last_block; ++b) {
        Ps[(max_num_blocks * c) + b] = blocks[b].p;
        Ms[(max_num_blocks * c) + b] = blocks[b].m;
        scores[(max_num_blocks * c) + b] = blocks[b].score;
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
    std::array<s32, static_cast<size_t>(kWordSize)> block_scores{};
    GetBlockCellValues(blocks[last_block], block_scores);
    const s32 best_score_tmp = block_scores[static_cast<size_t>(w)];
    if (best_score_tmp <= k) {
      return;
    }
  }
}

void Alignment::ObtainAlignmentTraceback(const s32 aligned_reference_length, u8* const aln, s32& aln_len) const {
  const auto best_score = edit_distance;
  const s32 max_num_blocks = CeilDiv(query_length, kWordSize);
  const s32 w = (max_num_blocks * kWordSize) - query_length;
  aln_len = 0;
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
  Word curr_p = Ps[(c * max_num_blocks) + b];
  // M of current block
  Word curr_m = Ms[(c * max_num_blocks) + b];
  // True if block to left exists and is in band
  bool there_is_left_block = (c > 0) && (b >= first_blocks[c - 1]) && (b <= last_blocks[c - 1]);
  // We set initial values of l_p and l_m to 0 only to avoid compiler warnings, they should not affect the
  // calculation as both l_p and l_m should be initialized at some moment later (but compiler can not
  // detect it since this initialization is guaranteed by "business" logic).
  Word l_p = 0;
  Word l_m = 0;
  if (there_is_left_block) {
    // P of block to the left
    l_p = Ps[((c - 1) * max_num_blocks) + b];
    // M of block to the left
    l_m = Ms[((c - 1) * max_num_blocks) + b];
  }
  curr_p <<= ToUnsigned(w);
  curr_m <<= ToUnsigned(w);
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
      left_score = scores[((c - 1) * max_num_blocks) + b];

      // Modified code (Roche)
      const s32 nr_bits{kWordSize - block_pos - 1};
      if (nr_bits != 0) {
        // The above loop is equivalent to a popcount of the highest bits of l_p and l_m.
        const auto mask{~((1UL << ToUnsigned(block_pos + 1)) - 1UL)};
        left_score -= std::popcount(l_p & mask);
        left_score += std::popcount(l_m & mask);
        l_p <<= ToUnsigned(nr_bits);
        l_m <<= ToUnsigned(nr_bits);
      }
    }
    if (upper_left_score == -1) {
      if (left_score != -1) {
        upper_left_score = left_score;
        if ((l_p & kHighBitMask) != 0U) {
          --upper_left_score;
        }
        if ((l_m & kHighBitMask) != 0U) {
          ++upper_left_score;
        }
      } else if ((c > 0) && ((b - 1) >= first_blocks[c - 1]) && ((b - 1) <= last_blocks[c - 1])) {
        // This is the case when upper left cell is last cell in block,
        // and block to left is not in band so left_score is -1.
        upper_left_score = scores[((c - 1) * max_num_blocks) + b - 1];
      }
    }
    if (upper_score == -1) {
      upper_score = current_score;
      if ((curr_p & kHighBitMask) != 0U) {
        --upper_score;
      }
      if ((curr_m & kHighBitMask) != 0U) {
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
      upper_score = -1;
      upper_left_score = -1;
      if (block_pos == 0) {  // If entering new (upper) block
        if (b == 0) {        // If there are no cells above (only boundary cells)
          // Move up
          aln[aln_len++] = kEdopInsert;
          for (s32 i = 0; i < c + 1; ++i) {  // Move left until end
            aln[aln_len++] = kEdopDelete;
          }
          break;
        }
        block_pos = kWordSize - 1;
        --b;
        curr_p = Ps[(c * max_num_blocks) + b];
        curr_m = Ms[(c * max_num_blocks) + b];
        if ((c > 0) && (b >= first_blocks[c - 1]) && (b <= last_blocks[c - 1])) {
          there_is_left_block = true;
          l_m = Ms[((c - 1) * max_num_blocks) + b];
          l_p = Ps[((c - 1) * max_num_blocks) + b];
        } else {
          there_is_left_block = false;
        }

      } else {
        --block_pos;
        l_p <<= 1;
        l_m <<= 1;
      }
      // Mark move
      aln[aln_len++] = kEdopInsert;
    } else {
      // Move left - deletion from reference - insertion to query
      if ((left_score != -1) && ((left_score + 1) == current_score)) {
        current_score = left_score;
        upper_score = upper_left_score;
        left_score = -1;
        upper_left_score = -1;
        --c;
        if (c == -1) {  // If there are no cells to the left (only boundary cells)
          // Move left
          aln[aln_len++] = kEdopDelete;
          const s32 num_up = (b * kWordSize) + block_pos + 1;
          for (s32 i = 0; i < num_up; ++i) {  // Move up until end
            aln[aln_len++] = kEdopInsert;
          }
          break;
        }
        curr_p = l_p;
        curr_m = l_m;
        if ((c > 0) && (b >= first_blocks[c - 1]) && (b <= last_blocks[c - 1])) {
          there_is_left_block = true;
          l_p = Ps[((c - 1) * max_num_blocks) + b];
          l_m = Ms[((c - 1) * max_num_blocks) + b];
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
        aln[aln_len++] = kEdopDelete;
      } else {
        // Move up left - (mis)match
        if (upper_left_score != -1) {
          const u8 move_code = upper_left_score == current_score ? kEdopMatch : kEdopMismatch;
          current_score = upper_left_score;
          upper_score = -1;
          left_score = -1;
          upper_left_score = -1;
          --c;
          if (c == -1) {  // If there are no cells to the left (only boundary cells)
            // Move left
            aln[aln_len++] = move_code;
            const s32 num_up = (b * kWordSize) + block_pos;
            for (s32 i = 0; i < num_up; ++i) {  // Move up until end
              aln[aln_len++] = kEdopInsert;
            }
            break;
          }
          if (block_pos == 0) {  // If entering upper left block
            if (b == 0) {        // If there are no more cells above (only boundary cells)
              // Move up left
              aln[aln_len++] = move_code;
              for (s32 i = 0; i < c + 1; ++i) {  // Move left until end
                aln[aln_len++] = kEdopDelete;
              }
              break;
            }
            block_pos = kWordSize - 1;
            --b;
            curr_p = Ps[(c * max_num_blocks) + b];
            curr_m = Ms[(c * max_num_blocks) + b];
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
            l_p = Ps[((c - 1) * max_num_blocks) + b];
            l_m = Ms[((c - 1) * max_num_blocks) + b];
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
          aln[aln_len++] = move_code;
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
    std::reverse(aln, aln + aln_len);
  }
}

void Alignment::ObtainAlignment(const std::span<const u8> align_reference, u8* const aln, s32& aln_len) {
  // Precondition: peq_table must already be populated for the current query.
  const s32 max_num_blocks = CeilDiv(query_length, kWordSize);
  const s32 w = (max_num_blocks * kWordSize) - query_length;

  MyersCalcEditDistanceNW(w, max_num_blocks, align_reference, /*disable_band_tightening=*/false);
  ObtainAlignmentTraceback(static_cast<s32>(align_reference.size()), aln, aln_len);
}

PairwiseAlignment::PairwiseAlignment(FlowContext& exec, size_t batch_number)
    : Task(exec, batch_number, fmt::format("Align{}", batch_number)) {
}

void Alignment::ComputeConsensus(FixedReadRecord& read, DuplexMetrics& metrics) {
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

    // Note: For mismatches in consensus we replace the base with the base from the reference sequence (i.e. R1)
    //
    // The output of the alignment is written to the FixedReadRecord object:
    // - The consensus sequence is built inside the dedicated consensus_seq field
    // - The new quality string is built the quality string Qual() field
    // - The when the consensus sequence is fully constructed the YC tag will be copied into the Seq() field
    //
    // Storing the YC tag string inside read.Seq() is a bit of a hack, but it allows us to reuse the space
    // as we don't need the original sequence anymore.
    // Notes on the YC tag:
    // - It is based on the alignment (showing directionality of changes of R1 relative to R2) in the consensus
    //   allowing for reconstruction of the original sequences (minus adapters)

    // Speed optimization: mark the base at the alignment length as an insert; this will force the loop track the
    // matching sections to stop, but that insert will be ignored as it is past the alignment length. And this
    // then will allow us to count the number of leading zero bytes in the alignment data, which can be implemented
    // efficiently using SIMD.
    // This mark is detected by the CopyMatchingBases() SIMD function, which will stop copying when it detects it
    _alignment[alignment_length] = kEdopInsert;

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

    // Strand detection and indel smearing.
    // Indel smearing is always applied to group scattered indels into fewer, longer gaps — even when
    // strand detection is disabled or the strand is unclassified.
    if (_params.strand_classifier.has_value()) {
      _strand = _params.strand_classifier->Classify(
          std::string_view(read.ConsensusSeq(), static_cast<size_t>(read.consensus_seq_len)));
    }
    // Here we both group indels or strand align them via indel smear (we do it here in either case)
    SmearAlignment(read);

    // Snapshot the raw sequence before GenerateConsensusQualAndIUPAC overwrites Seq() with the IUPAC tag.
    if (HasTrimmingDebugMode(_params.trimming_debug, TrimmingDebugMode::kRawSeq)) [[unlikely]] {
      read.duplex_raw_seq.assign(read.Seq(), read.SeqLen());
    }

    // Here we set quality string and we create the IUPAC-like tag
    GenerateConsensusQualAndIUPAC(read);

    // Store trimmed R1/R2 insert lengths for FormatDuplexRawSeqTag.
    if (HasTrimmingDebugMode(_params.trimming_debug, TrimmingDebugMode::kRawSeq)) [[unlikely]] {
      read.duplex_r1_len = static_cast<u32>(_trimmed_r1_len);
      read.duplex_r2_len = static_cast<u32>(_trimmed_r2_len);
    }

    // Everything is finished and not filtered. What should be true here:
    // - Consensus sequence string is stored in consensus_seq (and offset is set to trim off bases)
    // - Consensus quality string is stored in read.Qual()
    // - YC tag string is stored in read.Seq()
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
      metrics.failed_assigned_counts.raw_bases[sid] += read.SeqLen();
      metrics.RecordRawFailedReadLength(read.SeqLen());
      return;
    case kDuplexTooLongFail:
      ++metrics.failed_assigned_counts.consensus_too_long[sid];
      metrics.failed_assigned_counts.raw_bases[sid] += read.SeqLen();
      metrics.RecordRawFailedReadLength(read.SeqLen());
      return;
    case kTrimmedTooShortFail:
      ++metrics.failed_assigned_counts.trimmed_read_too_short[sid];
      metrics.failed_assigned_counts.raw_bases[sid] += read.SeqLen();
      metrics.RecordRawFailedReadLength(read.SeqLen());
      return;
    case kFailedMidadapterTrimFail:
      ++metrics.failed_assigned_counts.failed_hairpin_stem_trim_reads[sid];
      metrics.failed_assigned_counts.raw_bases[sid] += read.SeqLen();
      metrics.RecordRawFailedReadLength(read.SeqLen());
      return;
    default:
      // This should never happen so log an error
      throw error::Error("Alignment::UpdateMetrics: unexpected read status {} for read {}",
                         std::to_underlying(read.GetStatus()),
                         read.Name());
  }

  // update non-failure but strange cases metrics
  metrics.passing_counts.longer_r2[sid] += is_swapped ? 1 : 0;
  metrics.passing_counts.multi_candidate_alignments[sid] += _had_multi_candidates ? 1 : 0;

  const auto trim_size = static_cast<u32>(_endadapter_trim_pos);
  if (trim_size == 0) {
    ++metrics.passing_counts.no_endadapter[sid];
  }

  metrics.sample_endadapter_position_distr[sid].AddCountToHistogram(trim_size, 1);
  metrics.total_endadapter_position_distr.AddCountToHistogram(trim_size, 1);

  // update base count metrics
  metrics.passing_counts.concordant_bases[sid] += _concordant_and_ambiguous_bases_total - _ambiguous_bases_total;
  metrics.passing_counts.discordant_bases[sid] += static_cast<u64>(_discordant_bases_total);
  metrics.passing_counts.simplex_bases[sid] += static_cast<u64>(simplex_length);
  metrics.passing_counts.raw_bases[sid] += static_cast<u64>(read.SeqLen());
  metrics.passing_counts.ambiguous_bases[sid] += _ambiguous_bases_total;

  const auto report_length = static_cast<u32>(read.consensus_seq_len);
  const auto raw_length = static_cast<u32>(read.SeqLen());
  metrics.untrimmed_passing_length_distr.AddCountToHistogram(raw_length, 1);
  metrics.trimmed_passing_length_distr.AddCountToHistogram(report_length, 1);
  metrics.sample_untrimmed_duplex_passing_read_length_distr[sid].AddCountToHistogram(raw_length, 1);

  // update the metrics on the read length distributions
  if (simplex_length == 0) {
    // full duplex read
    metrics.passing_counts.longer_r2_full_duplex[sid] += is_swapped ? 1 : 0;
    metrics.passing_counts.full_duplex[sid] += 1;
    metrics.sample_trimmed_full_duplex_passing_read_length_distr[sid].AddCountToHistogram(report_length, 1);
    metrics.untrimmed_full_passing_length_distr.AddCountToHistogram(raw_length, 1);
    metrics.trimmed_full_passing_length_distr.AddCountToHistogram(report_length, 1);
    metrics.sample_untrimmed_full_duplex_passing_read_length_distr[sid].AddCountToHistogram(raw_length, 1);
  } else {
    // partial duplex read
    metrics.passing_counts.partial_duplex[sid] += 1;
    metrics.sample_trimmed_partial_duplex_passing_read_length_distr[sid].AddCountToHistogram(report_length, 1);
    metrics.untrimmed_partial_passing_length_distr.AddCountToHistogram(raw_length, 1);
    metrics.trimmed_partial_passing_length_distr.AddCountToHistogram(report_length, 1);
    metrics.sample_untrimmed_partial_duplex_passing_read_length_distr[sid].AddCountToHistogram(raw_length, 1);
  }
  metrics.sample_trimmed_duplex_passing_read_length_distr[sid].AddCountToHistogram(report_length, 1);

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
  const auto [trim_start, trim_end] = _trimmer.FindUMIPos(read);

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
        // This branch is only reached when consensus_seq_len >= simplex_length (see the condition above), so
        // new_alignment_length is always >= 0 and the index below is in bounds. The assert documents that
        // invariant for the static analyzer (SonarQube cpp:S3519) and guards against future regressions.
        assert(alignment_length >= 0 && "alignment_length must be non-negative before indexing 'alignment'");
        // mark the end of the alignment for the SIMD function used afterward
        _alignment[alignment_length] = kEdopInsert;
      }
    }
    // trim off the 5' side
    if (trim_start != -1) {
      // upstream should handle returning valid trim positions
      _endadapter_trim_pos = trim_start;
    } else if (read.consensus_seq_len == 0) {
      // The trim was possible but we shouldn't proceed with partial trimming as
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
  auto* const consensus_seq = read.consensus_buffer;

  // still in alphabet representation
  // original reference sequence IF NOT SWAPPED
  auto* const reference_in{read.Seq() + insert_start1};

  // these are in alphabet (0..3) representations
  const auto* const p_query{is_prefix ? reverse_aln_query : query};  // input data is now in correct output order,
  const auto* const p_reference{is_prefix ? reverse_aln_reference : reference};

  // We first copy non-duplex bases (i.e. unaligned data) into consensus_seq, though likely later assign a shorter
  // consensus_seq_len if we find the endadapter (to effectively trim it off).
  if (is_swapped) {  // The S2 is longer than S1, should only occur with full reads
    for (s32 i = 0; i < simplex_length; ++i) {
      // copy via alphabet representation
      consensus_seq[i] = sequence::kDnaAlphabet[p_reference[i]];
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
        simd::CopyMatchingSeq(_alignment + consensus_pos, reference_in + reference_pos, consensus_seq + output_size);
    if (stretch_length > 0) {
      // Copy out the perfectly aligned data.
      if (is_swapped) {  // slow path as we need to go via the alphabet
        for (s32 i = 0; i < stretch_length; ++i) {
          consensus_seq[output_size + i] = sequence::kDnaAlphabet[p_query[query_pos + i]];
        }
      }
      query_pos += stretch_length;
      reference_pos += stretch_length;
      consensus_pos += stretch_length;
      output_size += stretch_length;
      // This wraps up the perfectly aligned stretch, which should be the majority of the alignment data.
    } else {  // stretch length was 0, which happens when we ran out of data or encountered a mismatch/indel.
      // NOTE: a future optimization could mask only after a first pass to preserve the original
      // sequence for trimming, or perform start-adapter trimming before alignment.
      switch (_alignment[consensus_pos]) {
        case kEdopInsert:
          // Insertion into reference = deletion from query.
          consensus_seq[output_size] = mask_discordant_bases ? 'N' : sequence::kDnaAlphabet[p_query[query_pos]];
          ++query_pos;
          break;
        case kEdopDelete:
          // Deletion from reference = insertion to query.
          consensus_seq[output_size] = mask_discordant_bases ? 'N' : sequence::kDnaAlphabet[p_reference[reference_pos]];
          ++reference_pos;
          break;
        case kEdopMismatch: {
          // Mismatch = substitution. Note that we might have to swap the query and reference bases.
          const auto true_ref_base = is_swapped ? p_query[query_pos] : p_reference[reference_pos];
          consensus_seq[output_size] = mask_discordant_bases ? 'N' : sequence::kDnaAlphabet[true_ref_base];
          ++query_pos;
          ++reference_pos;
          break;
        }
        default:
          throw error::Error(
              "Alignment::GenerateConsensusSeq: unexpected alignment code {} at consensus position {} for read {}",
              static_cast<s32>(_alignment[consensus_pos]),
              consensus_pos,
              std::string_view(read.Name(), read.NameLen()));
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
        consensus_seq[output_size++] = sequence::kDnaAlphabet[p_reference[reference_pos]];
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
      throw error::Error("Alignment::GenerateConsensusSeq: unexpected duplex type {}",
                         std::to_underlying(_adapter_type));
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
        _endadapter_trim_pos,
        read.consensus_seq_len,
        std::string_view(read.Name(), read.NameLen()));
  } else {
    // Trim the consensus sequence to the start position
    read.consensus_seq_offset = static_cast<u32>(_endadapter_trim_pos);
    read.consensus_seq_len -= _endadapter_trim_pos;
  }
  return true;
}

void Alignment::SmearAlignment(const FixedReadRecord& read) {
  // isolate only the non-simplex region (untrimmed consensus region, indexes shared with alignment)
  const auto consensus_only =
      std::string_view(read.consensus_buffer + simplex_length, static_cast<size_t>(alignment_length));
  // Pre-pass: collapse indels toward one extreme, always needs to be done
  const std::span alignment_span{_alignment, consensus_only.size()};
  if (is_prefix) {
    SmearIndelsRight(consensus_only, alignment_span);
  } else {
    SmearIndelsLeft(consensus_only, alignment_span);
  }
  // If set we need a second pass for topology marking, but ultimate directionality doesn't matter
  const auto no_ambiguity = _params.ambiguous_topology_quality_mode == AmbiguousTopologyQualityMode::kNone;
  // If set we need a second pass to set strand direction
  bool strand_pass = false;
  // If we need to preserve direction of smear not
  bool preserve_direction = false;
  switch (_strand) {
    using enum strand::StrandType;
    case kForward:
    case kForwardSig:
      // strand direction need to be moved to the left if it isn't already
      // if prefix mode is turned on this mean it needs to be moved
      strand_pass = is_prefix;
      preserve_direction = true;
      break;
    case kReverse:
    case kReverseSig:
      // strand direction need to be moved to the right if it isn't already
      // if prefix mode is turned on this means it doesn't need to be moved
      strand_pass = !is_prefix;
      preserve_direction = true;
      break;
    default:
      // We don't care about final direction so any smearing doesn't matter
      break;
  }
  // Perform second pass
  if (no_ambiguity) {
    // No marking needed. Only redo the smear when the strand requires flipping the pre-pass
    // direction: the pre-pass smeared prefix->right / infix->left, so the flip is the opposite.
    if (strand_pass) {
      if (is_prefix) {
        SmearIndelsLeft(consensus_only, alignment_span);
      } else {
        SmearIndelsRight(consensus_only, alignment_span);
      }
    }
  } else {
    const std::span<u8> mask_span{_ambiguous_mask.data(), consensus_only.size()};

    // The mark pass always smears in the same direction regardless of strand (prefix marks left,
    // infix marks right); only the target buffer varies. When the strand is known and already in
    // the correct direction (preserve_direction but no strand flip needed), we must not disturb the
    // pre-pass result in _alignment, so mark into a scratch copy instead. Otherwise we can mark in
    // place: either strand smearing is happening anyway (strand_pass) or the strand is undetermined
    // and the final direction doesn't matter.
    std::span<u8> target_span = alignment_span;
    if (preserve_direction && !strand_pass) {
      auto* const smear_scratch = _alignment == _aln_buffer.data() ? _alt_aln_buffer.data() : _aln_buffer.data();
      target_span = std::span<u8>{smear_scratch, consensus_only.size()};
      std::copy_n(_alignment, consensus_only.size(), target_span.data());
    }

    if (is_prefix) {
      SmearMarkAllLeft(consensus_only, target_span, mask_span);
    } else {
      SmearMarkAllRight(consensus_only, target_span, mask_span);
    }
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
  _concordant_and_ambiguous_bases_total = 0;
  _discordant_bases_total = 0;
  _ambiguous_bases_total = 0;

  // Runtime-configurable output quality scores (defaults reproduce the historical constants).
  const auto& quals = _params.base_quality_scores;

  s32 consensus_pos = 0;
  s32 query_pos = 0;
  // start of duplex data
  s32 reference_pos = simplex_length;
  s32 output_size = 0;

  auto* const consensus_qual = read.Qual();

  // these are in alphabet (0..3) representations
  const auto* const p_query{is_prefix ? reverse_aln_query : query};  // input data is now in correct output order,
  const auto* const p_reference{is_prefix ? reverse_aln_reference : reference};

  // determine where the reference and query start if the read is full duplex
  // If the simplex region length is less than the endadapter trim position, we need to adjust
  // the reference and query start positions. This occurs when the read has a simplex region
  // (non-duplex region) that overlaps with the trimmed duplex region. The alignment must be
  // traversed to account for insertions and deletions, ensuring the correct start position.
  if (simplex_length < _endadapter_trim_pos) {
    const s32 loop_end_condition = _endadapter_trim_pos - simplex_length;
    for (; consensus_pos < loop_end_condition; ++consensus_pos) {
      // Process each alignment operation to adjust query and reference positions.
      switch (_alignment[consensus_pos]) {
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
              static_cast<s32>(_alignment[consensus_pos]),
              consensus_pos,
              std::string_view(read.Name(), read.NameLen()));
      }
    }
    simplex_length = 0;
  } else {
    // read is not full length, so our start positions are still equal to start_locations[num_locations - 1]
    simplex_length -= _endadapter_trim_pos;
    // fill simplex region with base quality
    std::fill_n(consensus_qual, simplex_length, quals.simplex);
    output_size = simplex_length;
  }

  // Debug-only snapshot used to format the raw-seq sq:Z tag.
  // Capture the consensus-traversal starting state before GenerateConsensusQualAndIUPAC() overwrites Seq().
  const s32 simplex_len_at_consensus_start = simplex_length;
  const s32 ref_pos_at_consensus_start = reference_pos;
  const s32 query_pos_at_consensus_start = query_pos;

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
    const auto* const r2 = is_swapped ? p_reference + reference_pos : p_query + query_pos;
    MethylConvert(read.ConsensusSeq() + simplex_length,
                  _alignment + consensus_pos,
                  alignment_length - consensus_pos,
                  iupac_cigar + iupac_length + simplex_length,
                  r2,
                  is_swapped);
    iupac_length += static_cast<u32>(read.consensus_seq_len);
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
    constexpr s32 kDecimalBase{10};
    s32 value = simplex_length;
    // buffer to hold the digits in reverse order; max 4 digits for simplex_length <= kMaxAlignLength (1504)
    std::array<char, 4> buffer{};
    size_t len = 0;
    do {
      // NOSONAR(S810) -- writing an ASCII decimal digit into a text buffer, not numeric storage.
      buffer[len++] = static_cast<char>('0' + (value % kDecimalBase));  // NOSONAR
      value /= kDecimalBase;
    } while (value > 0);
    // Append digits in correct order (buffer holds them least-significant-first).
    for (size_t i = 0; i < len; ++i) {
      iupac_cigar[iupac_length++] = buffer[len - 1 - i];
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
        simd::CopyMatchingQual(_alignment + consensus_pos, consensus_qual + output_size, quals.concordant);
    if (stretch_length > 0) {
      // Downgrade concordant bases to ambiguous where the mask is set
      if (_params.ambiguous_topology_quality_mode != AmbiguousTopologyQualityMode::kNone) {
        for (s32 i = 0; i < stretch_length; ++i) {
          if (_ambiguous_mask[static_cast<u32>(consensus_pos + i)] != 0) {
            consensus_qual[output_size + i] = quals.ambiguous;
            ++_ambiguous_bases_total;
          }
        }
      }
      query_pos += stretch_length;
      reference_pos += stretch_length;
      consensus_pos += stretch_length;
      _concordant_and_ambiguous_bases_total += static_cast<u32>(stretch_length);
      output_size += stretch_length;

      // The IUPAC code for a stretch of perfectly aligned data is the number of bases, i.e. the
      // stretch length. Instead of using itos(), itoa() or sprintf(), we will write out the digits
      // ourselves; this code works for values [1..999] and is very fast as I don't have to deal with
      // unnecessary error checking or memory allocation for strings.
      // NOTE: this hand-rolled integer-to-decimal writer trades safety for speed; a safer variant
      // could be substituted if profiling shows it is not on the critical path.
      constexpr s32 kDecimalBase{10};
      // Largest power of ten this fast path emits (values [1..999]); see the comment above.
      constexpr s32 kMaxDecimalUnit{100};
      s32 unit = kMaxDecimalUnit;
      while (stretch_length < unit) {
        unit /= kDecimalBase;
      }
      while (unit != 0) {
        const s32 digit = stretch_length / unit;
        // NOSONAR(S810) -- writing an ASCII decimal digit into the text CIGAR buffer, not numeric storage.
        iupac_cigar[iupac_length++] = static_cast<char>('0' + digit);  // NOSONAR
        stretch_length -= digit * unit;
        unit /= kDecimalBase;
      }
      // This wraps up the perfectly aligned stretch, which should be the majority of the alignment data.
    } else {
      // stretch length was 0, which happens when we ran out of data or encountered a mismatch/indel.
      bool mark_discordant = true;
      switch (_alignment[consensus_pos]) {
        case kEdopInsert:
          // Insertion into reference = deletion from query.
          iupac_cigar[iupac_length++] =
              is_swapped ? kIUPACTable[p_query[query_pos]][kIndelIndex] : kIUPACTable[kIndelIndex][p_query[query_pos]];
          ++query_pos;
          break;
        case kEdopDelete:
          // Deletion from reference = insertion to query.
          iupac_cigar[iupac_length++] = is_swapped ? kIUPACTable[kIndelIndex][p_reference[reference_pos]]
                                                   : kIUPACTable[p_reference[reference_pos]][kIndelIndex];
          ++reference_pos;
          break;
        case kEdopMismatch: {
          // Mismatch = substitution. Note that we might have to swap the query and reference bases.
          const s32 ref_base = is_swapped ? p_query[query_pos] : p_reference[reference_pos];
          const s32 query_base = is_swapped ? p_reference[reference_pos] : p_query[query_pos];
          iupac_cigar[iupac_length++] = kIUPACTable[ref_base][query_base];
          ++query_pos;
          ++reference_pos;
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
          mark_discordant = false;
          break;
        }
        default:
          throw error::Error("Alignment::GenerateConsensusQualAndIUPAC: unexpected alignment code {} for read {}",
                             static_cast<s32>(_alignment[consensus_pos]),
                             std::string_view(read.Name(), read.NameLen()));
      }
      // In kSmearOverride mode we may mark ambiguous, but only where the mask matches.
      const bool mark_ambiguous =
          _params.ambiguous_topology_quality_mode == AmbiguousTopologyQualityMode::kSmearOverride &&
          _ambiguous_mask[static_cast<size_t>(consensus_pos)] != 0;
      if (mark_ambiguous) {
        consensus_qual[output_size] = quals.ambiguous;
        ++_ambiguous_bases_total;
        ++_concordant_and_ambiguous_bases_total;
      } else if (mark_discordant) {
        consensus_qual[output_size] = quals.discordant;
        ++_discordant_bases_total;
      } else {
        consensus_qual[output_size] = quals.concordant;
        ++_concordant_and_ambiguous_bases_total;
      }
      ++output_size;
      ++consensus_pos;
    }
  }
  // finish off the consensus quality string if past end of alignment but only based on consensus read length
  if (output_size < read.consensus_seq_len) {
    // should only happen in non-prefix modes
    assert(!is_prefix);

    const auto remaining_bases = read.consensus_seq_len - output_size;
    std::fill_n(consensus_qual + output_size, remaining_bases, quals.discordant);
    _discordant_bases_total += remaining_bases;

    for (s32 i = 0; i < remaining_bases; ++i) {
      iupac_cigar[iupac_length++] = is_swapped ? kIUPACTable[kIndelIndex][p_reference[reference_pos]]
                                               : kIUPACTable[p_reference[reference_pos]][kIndelIndex];
      ++reference_pos;
    }
  }

  // Check iupac_length doesn't exceed the maximum allowed length
  if (iupac_length >= kMaxAlignLength + 2) {
    // +2 for the -0 or +0
    Logging::Warn(
        "Alignment::GenerateConsensusQualAndIUPAC: Comment length {} exceeds maximum allowed length {} for "
        "read {}",
        iupac_length,
        kMaxAlignLength,
        std::string_view(read.Name(), read.NameLen()));
    // Set too long failure status
    read.SetStatus(FixedReadRecord::Status::kDuplexTooLongFail);
  }

  // Debug-only trimmed insert lengths for sq:Z formatting.
  // Compute how far each strand advanced during consensus traversal, then map the reference/query strands back to R1/R2
  // depending on swap state.
  const s32 ref_trimmed_len = simplex_len_at_consensus_start + (reference_pos - ref_pos_at_consensus_start);
  const s32 query_trimmed_len = query_pos - query_pos_at_consensus_start;
  _trimmed_r1_len = is_swapped ? query_trimmed_len : ref_trimmed_len;
  _trimmed_r2_len = is_swapped ? ref_trimmed_len : query_trimmed_len;

  // Stored in the read.Seq() field, as we store the consensus sequence elsewhere so we can us to reuse the space.
  const auto* const iter = fmt::format_to(read.Seq(), "{}+", std::string_view(iupac_cigar, iupac_cigar + iupac_length));
  read.iupac_length = static_cast<u32>(iter - read.Seq());
}

void PairwiseAlignment::operator()() const {
  try {
    const StopWatch sw;

    const auto& batch{context.GetBatchData(batch_nr)};
    const auto& param{context.DemuxParam()};
    const auto error_rate{param.max_pairwise_error_rate};

    // Loop over all records in the batch and perform alignment.
    if (batch.records) {
      // Create an aligner object; this object uses a lot of stack space for intermediate alignment results, and
      // allocating space on the stack actually takes quite some time (I benchmarked it to be about 5% of the total
      // time). I obtained a significant speedup by creating the aligner object for every batch and re-using it rather
      // than assuming that stack allocation is "free". Use a max error rate of 10%, do not mask discordant bases.

      Alignment aligner(error_rate,
                        false,
                        AlignmentMode::kPrefix,
                        context.GetManager().DemuxObjectDuplex(),
                        context.DemuxParam(),
                        context.GetManager().AdapterTypeValue());

      auto& duplex_metrics{DuplexMetrics::Instance()};
      for (size_t i = 0; i < batch.Size(); ++i) {
        if (auto& read = (*(batch.records))[i]; read.GetStatus() == FixedReadRecord::Status::kDemultiplexed) {
          aligner.ComputeConsensus(read, duplex_metrics);
        }
      }
    }

    context.AddToAlignmentTime(static_cast<u64>(sw.ElapsedTime()));
    // NOSONAR(S1181) -- task boundary: catch every std::exception, log it, and hand it to the flow
    // manager via SetTaskException so one failed read cannot bring down the worker thread.
  } catch (const std::exception& e) {  // NOSONAR(S1181)
    Logging::Error("PairwiseAlignment::operator() failed: {}", e.what());
    SetTaskException(std::current_exception());
  }
}

// NOTE: If you change the parameters here check SIMDNA
Alignment::Alignment(const f32 max_error,
                     const bool mask_bases,
                     const AlignmentMode mode,  // NOLINT
                     const DemuxAndTrimDuplex& trimmer,
                     const DemuxAndTrimParam& params,  // NOLINT
                     const AdapterType adapter_type)   // NOLINT
    : max_error_rate(max_error),
      mask_discordant_bases(mask_bases),
      is_prefix{mode == AlignmentMode::kPrefix},
      _prefix(mode == AlignmentMode::kPrefix ? std::make_unique<PrefixState>() : nullptr),
      _infix(mode == AlignmentMode::kInfix ? std::make_unique<InfixState>() : nullptr),
      _trimmer(trimmer),
      _params(params),
      _adapter_type(adapter_type) {
}
}  // namespace xoos::demux
