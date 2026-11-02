#pragma once

#include <array>
#include <memory>
#include <span>

#include "adapter-design/adapter-design.h"
#include "adapters/duplex/demux-and-trim-duplex.h"
#include "core/demux-and-trim-pipeline.h"
#include "sequence/strand/strand.h"
#include "task/task.h"

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

  See alignment.cpp for implementation details.
*/

namespace xoos::demux {
// forward declaration
struct DuplexMetrics;
class FlowContext;

using Word = u64;
// Size of Word in bits
constexpr s32 kWordSize = sizeof(Word) * 8;
constexpr Word kWord1 = static_cast<Word>(1);
constexpr Word kHighBitMask = kWord1 << (kWordSize - 1);  // 100..00

struct Block {
  // Pvin
  Word p;
  // Mvin
  Word m;
  // score of last cell in block;
  s32 score;

  Block() = default;
};

struct Delta {
  s32 pos_hout;
  s32 minus_hout;
};

// Helper class for alignment. This is an adapted version of the edlib library, see alignment.cpp for details.
class Alignment {
 public:
  // Constructor.
  // - max_error: the maximum error rate (from 0 to 1) that is allowed for the alignment. By default, we allow 1
  //   discordant base for every 10 bases in the alignment (10% error rate).
  // - mask_bases: if true, the discordant bases will be masked in the alignment.
  Alignment(f32 max_error,
            bool mask_bases,
            AlignmentMode mode,
            const DemuxAndTrimDuplex& trimmer,
            const DemuxAndTrimParam& params,
            AdapterType adapter_type);

  // Aligns the two sequences. This is the main entry point for the alignment.
  void ComputeConsensus(FixedReadRecord& read, DuplexMetrics& metrics);

  void UpdateMetrics(FixedReadRecord& read, DuplexMetrics& metrics) const;

 protected:
  // Performs actual alignment. By default, we're using the prefix method (no edit distance penalty for query end);
  // if mode_prefix is false, we're using the infix method (no penalty at query start and end).
  bool Align();

  /**
   * @brief Post-alignment candidate selection for prefix (SHW) mode.
   *
   * When multiple tied optimal alignments exist (same edit distance), their gap topology can
   * differ. This selects the candidate with fewest total indel operations.
   *
   * @return true if alignment succeeded, false if consensus would exceed kMaxAlignLength.
   */
  bool SelectBestPrefixAlign();

  /**
   * @brief Post-alignment candidate selection for infix (HW) mode.
   *
   * Picks the alignment closest to the end of the reference (rightmost candidate).
   *
   * @return true if alignment succeeded, false if consensus would exceed kMaxAlignLength.
   */
  bool SelectBestInfixAlign();

  // Digest the input sequences for the alignment; this includes transforming the sequences into indices into an
  // alphabet and reversing/complementing the sequences as needed.
  bool SetSequences(const FixedReadRecord& read);

  // The use of this function depends on the alignment method; to maximize reuse, we're passing in function arguments.
  void MyersCalcEditDistanceSemiGlobal(const Word peq_input[],
                                       s32 w,
                                       s32 max_num_blocks,
                                       s32 k,
                                       bool do_prefix,
                                       const u8* target,
                                       s32 target_length,
                                       s32& score,
                                       s32* positions_out,
                                       u32& num_positions_out);

  /**
   * @brief Computes the NW DP matrix for the given query/reference using Myers' bit-vector algorithm,
   *        populating Ps/Ms/scores for traceback.
   *
   * When @p disable_band_tightening is true, Ukkonen band shrink/grow logic is skipped: every column writes
   * the full set of blocks (first_block = 0, last_block = max_num_blocks - 1). This is used by
   * SelectBestPrefixAlign so that a single fill at the widest end_location can serve tracebacks from any
   * narrower end_location — the per-endpoint band shapes don't nest, so a banded shared fill would leave
   * cells uninitialized that a narrower traceback would later read.
   *
   * @param[in] w                       Size of padding in last block.
   * @param[in] max_num_blocks          Number of blocks needed to cover the whole query.
   * @param[in] align_reference         Reference sequence slice to align against; its size is the
   *                                    reference length used by the DP.
   * @param[in] disable_band_tightening When true, skip Ukkonen band shrink/grow (see above).
   */
  void MyersCalcEditDistanceNW(s32 w,
                               s32 max_num_blocks,
                               std::span<const u8> align_reference,
                               bool disable_band_tightening);

  /**
   * @brief Finds one possible alignment that gives optimal score.
   *
   * @param[in] align_reference Reference sequence slice to align against.
   * @param[out] aln Destination buffer for the alignment operations.
   * @param[out] aln_len Length of the alignment opt array.
   */
  void ObtainAlignment(std::span<const u8> align_reference, u8* aln, s32& aln_len);

  /**
   * @brief Finds one possible optimal alignment via traceback through the dynamic programming matrix.
   *
   * Consumes large amount of memory: O(query_length * aligned_reference_length).
   *
   * @param[in] aligned_reference_length Normal length, without W.
   * @param[out] aln Destination buffer for the alignment operations.
   * @param[out] aln_len Length of the alignment opt array.
   */
  void ObtainAlignmentTraceback(s32 aligned_reference_length, u8* aln, s32& aln_len) const;

  /**
   * Corresponds to Advance_Block function from Myers.
   */
  static Delta CalculateBlock(Word pos_v,
                              Word minus_v,
                              Word eq,
                              Word pos_hin,
                              Word minus_hin,
                              Word& pos_hout,
                              Word& minus_hout);  // NOLINT

  static void GetBlockCellValues(const Block& block, std::span<s32, static_cast<size_t>(kWordSize)> scores);

  // Result of alignment
  s32 edit_distance{-1};  // If alignment was found, this value is the edit distance (>= 0), -1 otherwise.

  /**
   * Array of zero-based positions in reference where optimal alignment paths end.
   * If gap after query is penalized, gap counts as part of query (NW), otherwise not.
   */
  s32 end_locations[kMaxAlignLength];

  /**
   * Number of end (and start) locations.
   */
  u32 num_locations{0};

  s32 alignment_length{0};  // The alignment length is the number of duplex bases found
  s32 simplex_length{0};    // while simplex length is the number of non-duplex bases found

  // maximum error rate allowed for the alignment
  float max_error_rate;
  // if true, we will mask the discordant bases in the alignment.
  bool mask_discordant_bases;

  // lengths and positions of input strings used during intermediate calculations
  s32 query_length{0};
  s32 reference_length{0};
  s32 insert_start1{0};
  s32 insert_start2{0};
  s32 insert_end1{0};
  s32 insert_end2{0};
  // if true, the query and reference have been swapped.
  bool is_swapped{false};
  // if true, we're using the prefix method (no edit distance penalty for query end)
  bool is_prefix;

  // Various buffers - allocated on the stack to avoid heap memory allocation overhead.
  // These fixed-size scratch buffers are addressed via raw pointer arithmetic in the edlib-derived
  // Myers kernels (e.g. peq_input + b, Ps[(c * max_num_blocks) + b], bl = blocks + n). Switching to
  // std::array would force .data() at every hot-path use without changing the layout or bounds, so
  // the C arrays are kept intentionally (S5945 is scoped out for this file in spec.yaml).
  Word Ps[kMaxNrBlocks * kMaxAlignLength];  // NOLINT - let's keep this name
  Word Ms[kMaxNrBlocks * kMaxAlignLength];  // NOLINT - let's keep this name
  s32 scores[kMaxNrBlocks * kMaxAlignLength];
  s32 first_blocks[kMaxAlignLength];
  s32 last_blocks[kMaxAlignLength];

  u8 query[kMaxAlignLength];
  u8 reference[kMaxAlignLength];
  u8 reverse_aln_query[kMaxAlignLength];
  u8 reverse_aln_reference[kMaxAlignLength];
  Word peq_table[kPeqTableSize];
  s32 positions[kMaxAlignLength];
  char iupac_cigar[kMaxAlignLength];

 private:
  /**
   * Alignment is found for first pair of start and end locations.
   * Alignment is sequence of numbers: 0, 1, 2, 3.
   * 0 stands for match.
   * 1 stands for insertion to reference.
   * 2 stands for insertion to query.
   * 3 stands for mismatch.
   * Alignment aligns query to reference from beginning of query till end of query.
   * If gaps are not penalized, they are not in alignment.
   */
  std::array<u8, kMaxAlignLength> _aln_buffer{};
  /// @brief Second alignment buffer used when evaluating alternative candidates and as smear scratch when unselected.
  std::array<u8, kMaxAlignLength> _alt_aln_buffer{};
  /// @brief Ambiguity mask: 1 = position is topologically ambiguous, 0 = not.
  std::array<u8, kMaxAlignLength> _ambiguous_mask{};
  u8* _alignment{_aln_buffer.data()};

  /**
   * @brief Buffers used exclusively in prefix (SHW) candidate selection.
   *
   * Allocated only when prefix mode is active, avoiding dead memory in infix-only builds.
   */
  struct PrefixState {
    /// @brief Reusable scratch for the candidate-evaluation order in SelectBestPrefixAlign.
    std::array<u32, kMaxAlignLength> order{};
  };

  /// @brief Prefix-only state; null when operating in infix mode.
  std::unique_ptr<PrefixState> _prefix;

  /**
   * @brief Buffers used exclusively in infix (HW) alignment mode.
   *
   * Allocated only when infix mode is active, avoiding dead memory in prefix-only builds.
   */
  struct InfixState {
    /// @brief Zero-based positions in reference where optimal alignment paths start.
    std::array<s32, kMaxAlignLength> start_locations{};
    std::array<u8, kMaxAlignLength> reverse_query{};
    std::array<u8, kMaxAlignLength> reverse_reference{};
    std::array<Word, kPeqTableSize> reverse_peq_table{};
    std::array<s32, kMaxAlignLength> reverse_positions{};
  };

  /// @brief Infix-only state; null when operating in prefix mode.
  std::unique_ptr<InfixState> _infix;

  // Variables for metrics
  u32 _concordant_and_ambiguous_bases_total{0};
  s32 _discordant_bases_total{0};
  u32 _ambiguous_bases_total{0};

  // Debug-only trimmed R1/R2 insert lengths, used to format the raw-seq sq:Z tag
  s32 _trimmed_r1_len{0};
  s32 _trimmed_r2_len{0};

  // for tracking strand information
  strand::StrandType _strand{strand::StrandType::kUnknown};

  // whether this read had multiple tied alignment candidates
  bool _had_multi_candidates{false};

  // trimming position
  s32 _endadapter_trim_pos{0};

  // Use to trim start/end adapter
  const DemuxAndTrimDuplex& _trimmer;
  // Use to get any runtime constants or parameters needed for alignment.
  const DemuxAndTrimParam& _params;
  // The adapter type being used for this alignment
  const AdapterType _adapter_type;

  bool ProcessExtraTrim(FixedReadRecord& read);

  bool GenerateConsensusSeq(FixedReadRecord& read);

  void GenerateConsensusQualAndIUPAC(FixedReadRecord& read);

  /**
   * @brief Applies indel smearing to the alignment array with optional strand-aware directionality.
   *
   * Groups scattered indels into fewer, longer gaps by sliding them along the consensus sequence.
   * Always applies smearing regardless of strand classification — unclassified reads still benefit
   * from gap grouping.
   *
   * The two-pass smear (e.g. left then right) first groups scattered indels together, then the
   * second pass repositions the group to a consistent anchor side. The pass order determines which
   * side the grouped indels end up on. In prefix mode the alignment array is reversed relative to
   * infix mode, so the pass order swaps to produce the same grouping behavior relative to the read
   * structure. When strand detection is enabled, forward and reverse strands get opposite final
   * positioning so that both strands' indels land consistently relative to the read.
   *
   * Directionality depends on alignment mode and strand:
   *
   * Prefix mode (reversed array):
   *   Forward strand: SmearIndelsRight (group) → SmearIndelsLeft (anchor)
   *   Reverse strand: SmearIndelsRight only (groups and anchors in one pass)
   *   Unclassified:   SmearIndelsRight only (group, no anchoring)
   *
   * Infix mode (non-reversed array):
   *   Forward strand: SmearIndelsLeft only (groups and anchors in one pass)
   *   Reverse strand: SmearIndelsLeft (group) → SmearIndelsRight (anchor)
   *   Unclassified:   SmearIndelsLeft only (group, no anchoring)
   *
   * @param read The read containing the consensus buffer.
   */
  void SmearAlignment(const FixedReadRecord& read);
};

class PairwiseAlignment : public Task {
 public:
  PairwiseAlignment(FlowContext& exec, size_t batch_number);
  ~PairwiseAlignment() override = default;

  /// @brief This function is called by the task scheduler to execute the alignment task for a batch of reads.
  void operator()() const;
};
}  // namespace xoos::demux
