#pragma once

#include "xoos/types/float.h"
#include "xoos/types/int.h"
#include "xoos/types/vec.h"

namespace xoos::read_collapser {

// Q-score for positions that have only simplex support
constexpr u8 kSimplexFallbackQscore = 22;
// Q-score for positions that have low confidence majority bases (e.g. strand bias, low depth, etc.)
constexpr u8 kLowConfidenceFallbackQscore = 5;
// Q-score for singleton bases
constexpr u8 kSingletonFallbackQscore = 22;

// Legacy Q-score for simplex bases
constexpr u8 kLegacySimplexFallbackQscore = 18;
// Legacy Q-score for low confidence bases
constexpr u8 kLegacyLowConfidenceFallbackQscore = 0;
// Legacy Q-score for singleton bases
constexpr u8 kLegacySingletonFallbackQscore = 0;

using BaseQualMap = vec<vec<u8>>;

/// Compute base quality score from a pileup of reads.
class QscoreCalculator {
 public:
  static constexpr s32 kMaxReadsInCluster = 200;
  /**
   * @brief Compute base quality given the total number of reads and the number of majority reads.
   * @param num_total_reads Total number of reads
   * @param num_majority_reads Number of majority reads
   * @return u8 A numeric representation of the base quality
   */
  static u8 ComputeQualLegacy(u32 num_total_reads, u32 num_majority_reads);

 private:
  // Static variables
  static const BaseQualMap kQualMap;          //!< Cache hashmap
  static constexpr f64 kSbxErrorRate = 0.01;  //!< SBX error rate
  static constexpr u8 kMaxPhred = 40;         //!< Maximum phred score that we handle
  /**
   * @brief Compute and store base qualities for all possible configurations.
   * @return BaseQualMap All possible base qualities for different read numbers
   */
  static BaseQualMap InitializeQualMap();
};

// Minimum ratio of majority base to total bases to be considered a majority
constexpr f64 kQscoreCalculatorMinMajorityRatio = 0.5;
// Minimum Phred score for simplex base
constexpr f64 kQscoreCalculatorSimplexMinQscore = 10.0;
// Maximum Phred score for simplex base
constexpr f64 kQscoreCalculatorSimplexMaxQscore = 40.0;
// Minimum Phred score for duplex base
constexpr f64 kQscoreCalculatorDuplexMinQscore = 20.0;
// Maximum Phred score for duplex base
constexpr f64 kQscoreCalculatorDuplexMaxQscore = 50.0;

constexpr f64 kQscoreCalculatorSimplexHorizontalShift = 0.03;
constexpr f64 kQscoreCalculatorDuplexHorizontalShift = 0.06;
constexpr f64 kQscoreCalculatorSimplexVerticalShift = 45.0;
constexpr f64 kQscoreCalculatorDuplexVerticalShift = 45.0;

// Constant factor added to the denominator when calculating the majority ratio
// to penalize small simplex cluster sizes
constexpr f64 kQscoreCalculatorSimplexClusterSizePenalty = 1.4;
// Constant factor added to the denominator when calculating the majority ratio
// to penalize small duplex cluster sizes
constexpr f64 kQscoreCalculatorDuplexClusterSizePenalty = 0.7;

constexpr f64 kQscoreCalculatorSimplexCoefficient = 3;
constexpr f64 kQscoreCalculatorDuplexCoefficient = 3;
constexpr f64 kQscoreCalculatorSimplexConstant = 70.0;
constexpr f64 kQscoreCalculatorDuplexConstant = 70.0;

constexpr f64 kQscoreCalculatorSimplexMinObservedRatioA = 0.4;
constexpr f64 kQscoreCalculatorSimplexMinObservedRatioB = 0.3;
constexpr f64 kQscoreCalculatorSimplexMinObservedRatioC = 0.8;
constexpr f64 kQscoreCalculatorDuplexMinObservedRatioA = 0.4;
constexpr f64 kQscoreCalculatorDuplexMinObservedRatioB = 0.3;
constexpr f64 kQscoreCalculatorDuplexMinObservedRatioC = 0.8;

constexpr f64 kQscoreCalculatorMaxStrandBiasPenalty = 5.0;

u8 ComputeQualSimplex(u32 num_total_reads, u32 num_majority_reads);

u8 ComputeQualDuplex(u32 num_fwd_reads, u32 num_fwd_majority_reads, u32 num_rev_reads, u32 num_rev_majority_reads);

u8 ComputeQual(u32 fwd_depth, u32 fwd_majority, u32 rev_depth, u32 rev_majority);

}  // namespace xoos::read_collapser
