#pragma once

#include <optional>
#include <string>

#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::read_collapser {

/**
 * Depth/error summary for a single consensus strand: the maximum raw-read depth (tag D), the minimum raw-read depth
 * over covered positions (tag M), and the fraction of raw-read bases disagreeing with the consensus call (tag E).
 */
struct ConsensusStrandSummary {
  u32 max_depth{};
  u32 min_depth{};
  f32 error_rate{};
};

/**
 * The per-strand and combined summary auxiliary tags derived from the per-base depth/error vectors. `forward` is the
 * `a` strand, `reverse` the `b` strand, and `combined` the `c` (both-strand) consensus. This holds only the summary
 * tags; the per-base array tags (`ad`/`bd`/`ae`/`be`) are the ConsensusResult per-base vectors themselves.
 */
struct ConsensusSummary {
  ConsensusStrandSummary forward;
  ConsensusStrandSummary reverse;
  ConsensusStrandSummary combined;
};

/**
 * The ConsensusResult struct holds the final consensus sequence and quality scores.
 */
struct ConsensusResult {
  std::string sequence;
  vec<u8> quality_scores;
  vec<u32> depths;

  // Per-position raw-read depth contributing to each single-strand consensus (ad/bd tags).
  std::optional<vec<u32>> forward_per_base_depth;
  std::optional<vec<u32>> reverse_per_base_depth;
  // Per-position raw-read bases disagreeing with each single-strand consensus base (ae/be tags).
  std::optional<vec<u32>> forward_per_base_errors;
  std::optional<vec<u32>> reverse_per_base_errors;

  explicit ConsensusResult(size_t max_length);

  /**
   * Compute the per-strand and combined summary auxiliary tags (max/min depth and error rate) from the per-base
   * depth and error vectors. Requires the per-base vectors to be populated; throws xoos::error::Error otherwise.
   */
  ConsensusSummary ComputeSummary() const;

  /**
   * Compute the average family size of the cluster from which the consensus
   * sequence is called. It is the sum of the depths divided by the length of
   * the consensus sequence.
   */
  u32 MeanClusterSize() const;

  // Remove low depth bases from the beginning and end of the consensus sequence.
  void TrimEnds(u32 min_trim_read_support);

  ConsensusResult& operator+=(const ConsensusResult& other);

  friend ConsensusResult operator+(const ConsensusResult& left, const ConsensusResult& right) {
    ConsensusResult result = left;
    result += right;
    return result;
  }
};

}  // namespace xoos::read_collapser
