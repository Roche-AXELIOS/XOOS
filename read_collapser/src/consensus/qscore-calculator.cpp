#include "consensus/qscore-calculator.h"

#include <algorithm>
#include <cmath>

#include <boost/math/distributions/binomial.hpp>

namespace xoos::read_collapser {

const BaseQualMap QscoreCalculator::kQualMap = QscoreCalculator::InitializeQualMap();

BaseQualMap QscoreCalculator::InitializeQualMap() {
  BaseQualMap qual_map(kMaxReadsInCluster + 1);
  for (u32 n = 1; n <= kMaxReadsInCluster; n++) {
    qual_map[n].resize(n + 1);
    for (u32 k = 1; k <= n; ++k) {
      auto d = boost::math::binomial_distribution<>{static_cast<f64>(n), kSbxErrorRate};
      auto probability = static_cast<u32>(-10 * std::log10(boost::math::pdf(d, k)));
      probability = std::min(probability, static_cast<u32>(std::numeric_limits<u8>::max()));
      qual_map[n][k] = std::min(kMaxPhred, static_cast<u8>(probability));
    }
  }
  return qual_map;
}

u8 QscoreCalculator::ComputeQualLegacy(u32 num_total_reads, u32 num_majority_reads) {
  // Result is not cached for large values of num_total_reads and num_majority_reads
  // so we compute it on the fly.
  if (num_total_reads > kMaxReadsInCluster || num_majority_reads > kMaxReadsInCluster) {
    auto d = boost::math::binomial_distribution<>{static_cast<f64>(num_total_reads), kSbxErrorRate};
    auto probability = static_cast<u32>(-10 * std::log10(boost::math::pdf(d, num_majority_reads)));
    probability = std::min(probability, static_cast<u32>(std::numeric_limits<u8>::max()));
    return std::min(kMaxPhred, static_cast<u8>(probability));
  }
  return kQualMap.at(num_total_reads).at(num_majority_reads);
}

u8 ComputeQualSimplex(const u32 num_total_reads, const u32 num_majority_reads) {
  const f64 ratio = static_cast<f64>(num_majority_reads) /
                    static_cast<f64>(num_total_reads + kQscoreCalculatorSimplexClusterSizePenalty);
  if (ratio < kQscoreCalculatorMinMajorityRatio) {
    // No reliable majority if the ratio is below 0.5
    return kLowConfidenceFallbackQscore;
  }
  // The majority only becomes reliable when the ratio of majority reads to total reads is above
  // a certain threshold, and this threshold scales with the total number of reads
  // This minimum ratio is the minimum majority ratio observed in actual data (i.e. there are either no
  // or very few consensus base calls with a majority ratio below this threshold in the actual data)
  const f64 min_observed_ratio =
      std::max(kQscoreCalculatorMinMajorityRatio,
               std::min(kQscoreCalculatorSimplexMinObservedRatioA * std::log10(num_total_reads) +
                            kQscoreCalculatorSimplexMinObservedRatioB,
                        kQscoreCalculatorSimplexMinObservedRatioC));

  // This is a constant derived from empirical data
  const f64 m = kQscoreCalculatorSimplexCoefficient * num_total_reads + kQscoreCalculatorSimplexConstant;

  f64 qscore = 0;
  if (ratio < min_observed_ratio) {
    const f64 qscore_at_min_ratio = m * std::log10(min_observed_ratio + kQscoreCalculatorSimplexHorizontalShift) +
                                    kQscoreCalculatorSimplexVerticalShift;
    // If the ratio is below the min observed ratio threshold, we extrapolate the quality score linearly using the
    // following formula q = (q_at_min_ratio / (min_ratio - 0.5)) * (ratio - 0.5)
    qscore = (qscore_at_min_ratio / (min_observed_ratio - kQscoreCalculatorMinMajorityRatio)) *
             (ratio - kQscoreCalculatorMinMajorityRatio);
  } else {
    // If the ratio is above the threshold, we use a different formula
    // q = m * log10(ratio + 0.03) + 45
    // Once a reliable majority is established, each additional read that supports the majority base
    // increases the confidence in the majority base, but the rate of increase decreases with the number of
    // additional reads, which can be modeled as a logarithmic function
    qscore = m * std::log10(ratio + kQscoreCalculatorSimplexHorizontalShift) + kQscoreCalculatorSimplexVerticalShift;
  }
  // Round the quality score to the nearest multiple of 10 and ensure it is within the valid range of [10, 40]
  qscore = std::clamp(qscore, kQscoreCalculatorSimplexMinQscore, kQscoreCalculatorSimplexMaxQscore);
  // Round to the nearest multiple of 10
  return static_cast<u8>(std::round(qscore / 10.0) * 10);
}

u8 ComputeQualDuplex(const u32 num_fwd_reads,
                     const u32 num_fwd_majority_reads,
                     const u32 num_rev_reads,
                     const u32 num_rev_majority_reads) {
  const u32 num_total_reads = num_fwd_reads + num_rev_reads;
  const u32 num_majority_reads = num_fwd_majority_reads + num_rev_majority_reads;
  const f64 ratio = static_cast<f64>(num_majority_reads) /
                    static_cast<f64>(num_total_reads + kQscoreCalculatorDuplexClusterSizePenalty);
  if (ratio < kQscoreCalculatorMinMajorityRatio) {
    // No reliable majority if the ratio is below 0.5
    return kLowConfidenceFallbackQscore;
  }
  // The majority only becomes reliable when the ratio of majority reads to total reads is above
  // a certain threshold, and this threshold scales with the total number of reads
  // This minimum ratio is the minimum majority ratio observed in actual data (i.e. there are either no
  // or very few consensus base calls with a majority ratio below this threshold in the actual data)
  const f64 min_observed_ratio =
      std::max(kQscoreCalculatorMinMajorityRatio,
               std::min(kQscoreCalculatorDuplexMinObservedRatioA * std::log10(num_total_reads) +
                            kQscoreCalculatorDuplexMinObservedRatioB,
                        kQscoreCalculatorDuplexMinObservedRatioC));

  // This is a constant derived from empirical data
  const f64 m = kQscoreCalculatorDuplexCoefficient * num_total_reads + kQscoreCalculatorDuplexConstant;

  f64 qscore = 0;
  if (ratio < min_observed_ratio) {
    // If the ratio is below the min observed ratio threshold, we extrapolate the quality score linearly using the
    // following formula q = (q_at_min_ratio / (min_ratio - 0.5)) * (ratio - 0.5)
    const f64 qscore_at_min_ratio = m * std::log10(min_observed_ratio + kQscoreCalculatorDuplexHorizontalShift) +
                                    kQscoreCalculatorDuplexVerticalShift;
    qscore = (qscore_at_min_ratio / (min_observed_ratio - kQscoreCalculatorMinMajorityRatio)) *
             (ratio - kQscoreCalculatorMinMajorityRatio);
  } else {
    // If the ratio is above the threshold, we use a different formula
    // q = m * log10(ratio + 0.06) + 45
    // Once a reliable majority is established, each additional read that supports the majority base
    // increases the confidence in the majority base, but the rate of increase decreases with the number of
    // additional reads, which can be modeled as a logarithmic function
    qscore = m * std::log10(ratio + kQscoreCalculatorDuplexHorizontalShift) + kQscoreCalculatorDuplexVerticalShift;
    // Penalize strand bias. This accounts for the fact that clusters with higher strand bias often have lower quality
    // scores than those with balanced strand support. This `strand_bias` value is a number between 0 and 1 that
    // measures how unbalanced the support is between the forward and reverse strands.
    const f64 strand_bias = (static_cast<f64>(std::max(num_fwd_reads, num_rev_reads) * 2) /
                             static_cast<f64>(num_fwd_reads + num_rev_reads)) -
                            1;
    qscore -= kQscoreCalculatorMaxStrandBiasPenalty * strand_bias;
  }
  // Round the quality score to the nearest multiple of 10 and ensure it is within the valid range of [20, 50]
  qscore = std::clamp(qscore, kQscoreCalculatorDuplexMinQscore, kQscoreCalculatorDuplexMaxQscore);
  // Round to the nearest multiple of 10
  return static_cast<u8>(std::round(qscore / 10.0) * 10);
}

u8 ComputeQual(const u32 fwd_depth, const u32 fwd_majority, const u32 rev_depth, const u32 rev_majority) {
  const f64 fwd_ratio = static_cast<f64>(fwd_majority) / static_cast<f64>(fwd_depth);
  const f64 rev_ratio = static_cast<f64>(rev_majority) / static_cast<f64>(rev_depth);
  // Use simplex qscore model if
  // 1. Not both strands support the majority base with a frequency of at least 50%
  // 2. Cluster size is less than 5
  if (fwd_ratio >= 0.5 && rev_ratio >= 0.5 && fwd_depth + rev_depth > 5) {
    // Both strands support the majority base with a frequency of at least 50%
    // We can use the duplex quality score calculation
    return ComputeQualDuplex(fwd_depth, fwd_majority, rev_depth, rev_majority);
  } else {
    return ComputeQualSimplex(fwd_depth + rev_depth, fwd_majority + rev_majority);
  }
}

}  // namespace xoos::read_collapser
