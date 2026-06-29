#pragma once

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
  // - max_error: the maximum error rate (in percent) that is allowed for the alignment. By default, we allow 1
  //   discordant base for every 10 bases in the alignment (10% error rate).
  // - mask_bases: if true, the discordant bases will be masked in the alignment.
  Alignment(float max_error, bool mask_bases, AlignmentMode mode, const DemuxAndTrimDuplex& trimmer,
            const DemuxAndTrimParam& params, AdapterType adapter_type);

  // Aligns the two sequences. This is the main entry point for the alignment.
  void ComputeConsensus(FixedReadRecord& read, DuplexMetrics& metrics);

  void UpdateMetrics(FixedReadRecord& read, DuplexMetrics& metrics) const;

  // Reporter functions.
  auto ConcordantDuplexBasesTotal() const { return _concordant_duplex_bases_total; }

  auto DiscordantBasesTotal() const { return _discordant_bases_total; }

 protected:
  // Performs actual alignment. By default, we're using the prefix method (no edit distance penalty for query end);
  // if mode_prefix is false, we're using the infix method (no penalty at query start and end).
  bool Align();

  // Digest the input sequences for the alignment; this includes transforming the sequences into indices into an
  // alphabet and reversing/complementing the sequences as needed.
  bool SetSequences(const FixedReadRecord& read);

  // The use of this function depends on the alignment method; to maximize reuse, we're passing in function arguments.
  void MyersCalcEditDistanceSemiGlobal(const Word peq_input[], s32 w, s32 max_num_blocks, s32 k, bool do_prefix,
                                       const u8* target, s32 target_length, s32& score, s32* positions_out,
                                       s32& num_positions_out);

  void MyersCalcEditDistanceNW(s32 w, s32 max_num_blocks, const u8* align_reference, s32 align_reference_length);

  void ObtainAlignment(const u8* align_reference, s32 aln_reference_length);

  void ObtainAlignmentTraceback(s32 aligned_reference_length);

  /**
   * Corresponds to Advance_Block function from Myers.
   */
  static Delta CalculateBlock(Word pos_v, Word minus_v, Word eq, Word pos_hin, Word minus_hin, Word& pos_hout,
                              Word& minus_hout);  // NOLINT

  static void GetBlockCellValues(const Block& block, s32* scores);

  // Result of alignment
  s32 edit_distance{-1};  // If alignment was found, this value is the edit distance (>= 0), -1 otherwise.

  /**
   * Array of zero-based positions in reference where optimal alignment paths end.
   * If gap after query is penalized, gap counts as part of query (NW), otherwise not.
   */
  s32 end_locations[kMaxAlignLength];

  /**
   * Array of zero-based positions in reference where optimal alignment paths start,
   * they correspond to end_locations.
   */
  s32 start_locations[kMaxAlignLength];

  /**
   * Number of end (and start) locations.
   */
  s32 num_locations{0};

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
  u8 alignment[kMaxAlignLength];

  s32 alignment_length{0};  // The alignment length is the number of duplex bases found
  s32 simplex_length{0};    // while simplex length is the number of non-duplex bases found

  // maximum error rate allowed for the alignment
  float max_error_rate;
  // if true, we will mask the discordant bases in the alignment.
  bool mask_discordant_bases;

  // lengths and positions of input strings used during intermediate calculations
  s32 query_length;
  s32 reference_length;
  s32 insert_start1;
  s32 insert_start2;
  s32 insert_end1;
  s32 insert_end2;
  // if true, the query and reference have been swapped.
  bool is_swapped;
  // if true, we're using the prefix method (no edit distance penalty for query end)
  bool is_prefix;

  // Various buffers - allocated on the stack to avoid heap memory allocation overhead
  Word Ps[kMaxNrBlocks * kMaxAlignLength];  // NOLINT - let's keep this name
  Word Ms[kMaxNrBlocks * kMaxAlignLength];  // NOLINT - let's keep this name
  s32 scores[kMaxNrBlocks * kMaxAlignLength];
  s32 first_blocks[kMaxAlignLength];
  s32 last_blocks[kMaxAlignLength];

  u8 query[kMaxAlignLength];
  u8 reference[kMaxAlignLength];
  u8 reverse_query[kMaxAlignLength];
  u8 reverse_reference[kMaxAlignLength];
  u8 reverse_aln_query[kMaxAlignLength];
  u8 reverse_aln_reference[kMaxAlignLength];
  Word peq_table[kPeqTableSize];
  Word reverse_peq_table[kPeqTableSize];
  s32 positions[kMaxAlignLength];
  s32 reverse_positions[kMaxAlignLength];
  char iupac_cigar[kMaxAlignLength];
  AlignmentMode alignment_mode;

 private:
  // Variables for metrics
  u32 _concordant_duplex_bases_total{0};
  s32 _discordant_bases_total{0};

  strand::StrandType _strand{strand::StrandType::kUnknown};  // for tracking strand information

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

  void ProcessStrand(const FixedReadRecord& read);
};

/// @brief Pairwise alignment of HD Mutex data.
class PairwiseAlignment : public Task {
 public:
  PairwiseAlignment(FlowContext& exec, size_t batch_nr);
  ~PairwiseAlignment() override = default;

  /// @brief This function is called by the task scheduler to execute the alignment task for a batch of reads.
  void operator()() const;
};
}  // namespace xoos::demux
