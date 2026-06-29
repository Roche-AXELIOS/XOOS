#include "consensus/consensus-result.h"

#include <numeric>

#include <xoos/types/int.h>

namespace xoos::read_collapser {

ConsensusResult::ConsensusResult(size_t max_length) {
  sequence.reserve(max_length);
  quality_scores.reserve(max_length);
  depths.reserve(max_length);
}

u32 ConsensusResult::MeanClusterSize() {
  u32 counted_positions = depths.size();
  if (counted_positions == 0) {
    return 0;
  }
  u32 total_bases = std::accumulate(depths.begin(), depths.end(), 0u);
  return (total_bases + counted_positions / 2) / counted_positions;
}

void ConsensusResult::TrimEnds(const u32 min_trim_read_support) {
  if (depths.empty() || min_trim_read_support <= 1) {
    return;
  }
  size_t trimmed_start = 0;
  size_t trimmed_end = depths.size() - 1;
  while (trimmed_start < depths.size() && depths[trimmed_start] < min_trim_read_support) {
    ++trimmed_start;
  }
  while (trimmed_end > trimmed_start && depths[trimmed_end] < min_trim_read_support) {
    --trimmed_end;
  }
  // If the trimmed range is empty, clear the sequence and return
  if (trimmed_start >= trimmed_end) {
    sequence.clear();
    quality_scores.clear();
    depths.clear();
    return;
  }
  // Create new vectors for the trimmed sequence, quality scores, and depths
  std::string trimmed_sequence = sequence.substr(trimmed_start, trimmed_end - trimmed_start + 1);
  vec<u8> trimmed_quality_scores(quality_scores.begin() + static_cast<s32>(trimmed_start),
                                 quality_scores.begin() + static_cast<s32>(trimmed_end) + 1);
  vec<u32> trimmed_depths(depths.begin() + static_cast<s32>(trimmed_start),
                          depths.begin() + static_cast<s32>(trimmed_end) + 1);
  // Update the member variables with the trimmed values
  sequence = std::move(trimmed_sequence);
  quality_scores = std::move(trimmed_quality_scores);
  depths = std::move(trimmed_depths);
  if (forward_per_base_depth.has_value()) {
    forward_per_base_depth = vec<u32>(forward_per_base_depth->begin() + static_cast<s32>(trimmed_start),
                                      forward_per_base_depth->begin() + static_cast<s32>(trimmed_end) + 1);
  }
  if (reverse_per_base_depth.has_value()) {
    reverse_per_base_depth = vec<u32>(reverse_per_base_depth->begin() + static_cast<s32>(trimmed_start),
                                      reverse_per_base_depth->begin() + static_cast<s32>(trimmed_end) + 1);
  }
  if (forward_per_base_majority_count.has_value()) {
    forward_per_base_majority_count =
        vec<u32>(forward_per_base_majority_count->begin() + static_cast<s32>(trimmed_start),
                 forward_per_base_majority_count->begin() + static_cast<s32>(trimmed_end) + 1);
  }
  if (reverse_per_base_majority_count.has_value()) {
    reverse_per_base_majority_count =
        vec<u32>(reverse_per_base_majority_count->begin() + static_cast<s32>(trimmed_start),
                 reverse_per_base_majority_count->begin() + static_cast<s32>(trimmed_end) + 1);
  }
}

ConsensusResult& ConsensusResult::operator+=(const ConsensusResult& other) {
  sequence += other.sequence;
  quality_scores.insert(quality_scores.end(), other.quality_scores.begin(), other.quality_scores.end());
  depths.insert(depths.end(), other.depths.begin(), other.depths.end());
  return *this;
}

}  // namespace xoos::read_collapser
