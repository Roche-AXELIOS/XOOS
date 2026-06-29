#pragma once

#include <optional>
#include <string>

#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::read_collapser {

/**
 * The ConsensusResult struct holds the final consensus sequence and quality scores.
 */
struct ConsensusResult {
  std::string sequence;
  vec<u8> quality_scores;
  vec<u32> depths;

  std::optional<vec<u32>> forward_per_base_depth{};
  std::optional<vec<u32>> reverse_per_base_depth{};
  std::optional<vec<u32>> forward_per_base_majority_count{};
  std::optional<vec<u32>> reverse_per_base_majority_count{};

  explicit ConsensusResult(size_t max_length);

  /**
   * Compute the average family size of the cluster from which the consensus
   * sequence is called. It is the sum of the depths divided by the length of
   * the consensus sequence.
   */
  u32 MeanClusterSize();

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
