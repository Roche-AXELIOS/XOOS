#pragma once

#include <array>

#include "consensus/base-encoder.h"
#include "consensus/consensus-matrix.h"
#include "core/read-collapser-options.h"
#include "xoos/types/int.h"

namespace xoos::read_collapser {

using BaseCounts = std::array<u32, kBaseIndexCount>;

// Struct to hold base counts for a single position (column) in the consensus matrix
struct ColumnBaseCounts {
  // Base counts of A, C, G, T, Gap for forward and reverse strands and overall
  BaseCounts fwd_base_counts{};
  BaseCounts rev_base_counts{};
  BaseCounts overall_base_counts{};

  // Only used for parent-daughter
  BaseCounts r1_base_counts{};
  BaseCounts r2_base_counts{};

  // Number of reads supporting this position
  u32 read_support{};
  u32 fwd_read_support{};
  u32 rev_read_support{};

  // Helper to reset all counts for the next column iteration
  void Reset();
};

// Result struct for majority voting on a single column
struct ColumnConsensusResult {
  char majority_base{kBaseN};
  u8 qscore{};
  u32 read_support{};

  struct ColumnMetrics {
    // Metrics needed for optional FASTQ tags and quality score calculation
    u32 fwd_majority_count{};
    u32 rev_majority_count{};
    u32 fwd_read_support{};
    u32 rev_read_support{};

    // Additional metrics needed for majority voting
    char fwd_majority_base{kBaseN};
    char rev_majority_base{kBaseN};
  };

  // Status flags needed for majority voting and Q-score fallback logic
  struct ColumnVotingFlags {
    // General voting flags
    bool is_singleton{false};
    bool not_majority{false};
    bool not_gap_majority{false};
    bool has_strand_bias{false};

    // Flags related to deconvolution
    bool deconvolution_enabled{false};
    bool apply_hp_quality_override{false};
    bool is_hp_insertion_tie{false};
    bool parent_daughter_discordant{false};
    // Assumed true until proven duplex concordant
    bool is_simplex_base{true};
    bool all_duplex_bases_are_discordant_or_simplex{true};
  };

  ColumnMetrics metrics{};
  ColumnVotingFlags flags{};
};

// Class to perform majority voting on a single column of the consensus matrix
class ColumnMajorityVotingWorker {
 public:
  // Takes options to use thresholds and quality models
  explicit ColumnMajorityVotingWorker(const ReadCollapserOptions& options);

  // Computes the consensus result for one column in the consensus matrix
  void ComputeColumnConsensusResult(const ConsensusMatrix& matrix, size_t column_index, ColumnConsensusResult& result);

 private:
  const f64 _gap_majority_ratio;
  const f64 _majority_ratio;
  const HDDeconvolutionType _hd_deconvolution_type;
  const bool _enable_legacy_qscore_model;

  // Internal state to avoid reallocating count buffers repeatedly
  ColumnBaseCounts _counts;

  // Step 1: Accumulate counts for the current column
  void AccumulateCounts(const ConsensusMatrix& matrix, size_t column_index);

  // Step 2: Determine the majority base and check initial flags (singleton, majority ratio)
  void DetermineMajorityAndInitialFlags(const ConsensusMatrix& matrix,
                                        size_t column_index,
                                        ColumnConsensusResult& result) const;

  // Step 3: Revert the gap call if low confidence or strand bias is detected
  void ApplyGapRules(ColumnConsensusResult& result);

  // Step 4: Set duplex and homopolymer related flags
  void SetDuplexAndHomopolymerFlags(const ConsensusMatrix& matrix,
                                    size_t column_index,
                                    ColumnConsensusResult& result) const;

  // Step 5: Calculate final quality score based on all conditions
  void CalculateQualityScore(ColumnConsensusResult& result);
};

}  // namespace xoos::read_collapser
