#include "consensus/consensus-result.h"

#include <algorithm>
#include <numeric>

#include <xoos/error/error.h>
#include <xoos/types/int.h>

namespace xoos::read_collapser {

ConsensusResult::ConsensusResult(size_t max_length) {
  sequence.reserve(max_length);
  quality_scores.reserve(max_length);
  depths.reserve(max_length);
}

ConsensusSummary ConsensusResult::ComputeSummary() const {
  if (!forward_per_base_depth.has_value() || !reverse_per_base_depth.has_value() ||
      !forward_per_base_errors.has_value() || !reverse_per_base_errors.has_value()) {
    throw error::Error("ComputeSummary requires per-base depth and error vectors to be populated");
  }

  // Summarize one strand's per-base depth and error vectors into a max/min depth and error rate. Positions with zero
  // depth are skipped for the min/max so a strand that does not cover part of the consensus does not force the minimum
  // to zero; the error rate is total errors divided by total depth.
  const auto summarize = [](const vec<u32>& depths, const vec<u32>& errors) {
    ConsensusStrandSummary summary{};
    u64 total_depth = 0;
    u64 total_errors = 0;
    for (size_t i = 0; i < depths.size(); ++i) {
      total_depth += depths[i];
      total_errors += errors[i];
      if (depths[i] == 0) {
        continue;
      }
      summary.max_depth = std::max(summary.max_depth, depths[i]);
      summary.min_depth = summary.min_depth == 0 ? depths[i] : std::min(summary.min_depth, depths[i]);
    }
    if (total_depth > 0) {
      summary.error_rate = static_cast<f32>(static_cast<f64>(total_errors) / static_cast<f64>(total_depth));
    }
    return summary;
  };

  // The combined summary sums forward and reverse depth/errors at each position to derive the cD/cM/cE tags from both
  // single strands.
  const size_t length = forward_per_base_depth->size();
  if (reverse_per_base_depth->size() != length || forward_per_base_errors->size() != length ||
      reverse_per_base_errors->size() != length) {
    throw error::Error("ComputeSummary requires all per-base vectors to have equal lengths");
  }
  vec<u32> combined_depths(length);
  vec<u32> combined_errors(length);
  for (size_t i = 0; i < length; ++i) {
    combined_depths[i] = (*forward_per_base_depth)[i] + (*reverse_per_base_depth)[i];
    combined_errors[i] = (*forward_per_base_errors)[i] + (*reverse_per_base_errors)[i];
  }

  return ConsensusSummary{.forward = summarize(*forward_per_base_depth, *forward_per_base_errors),
                          .reverse = summarize(*reverse_per_base_depth, *reverse_per_base_errors),
                          .combined = summarize(combined_depths, combined_errors)};
}

u32 ConsensusResult::MeanClusterSize() const {
  u32 counted_positions = depths.size();
  if (counted_positions == 0) {
    return 0;
  }
  u32 total_bases = std::accumulate(depths.begin(), depths.end(), 0U);
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
  if (trimmed_start > trimmed_end) {
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
  if (forward_per_base_errors.has_value()) {
    forward_per_base_errors = vec<u32>(forward_per_base_errors->begin() + static_cast<s32>(trimmed_start),
                                       forward_per_base_errors->begin() + static_cast<s32>(trimmed_end) + 1);
  }
  if (reverse_per_base_errors.has_value()) {
    reverse_per_base_errors = vec<u32>(reverse_per_base_errors->begin() + static_cast<s32>(trimmed_start),
                                       reverse_per_base_errors->begin() + static_cast<s32>(trimmed_end) + 1);
  }
}

ConsensusResult& ConsensusResult::operator+=(const ConsensusResult& other) {
  sequence += other.sequence;
  quality_scores.insert(quality_scores.end(), other.quality_scores.begin(), other.quality_scores.end());
  depths.insert(depths.end(), other.depths.begin(), other.depths.end());
  auto merge_optional_vec =
      [](std::optional<vec<u32>>& dst, const std::optional<vec<u32>>& src, const std::string_view name) {
        if (dst.has_value() != src.has_value()) {
          throw error::Error("Cannot concatenate ConsensusResults with inconsistent {} presence", name);
        }
        if (dst.has_value()) {
          dst->insert(dst->end(), src->begin(), src->end());
        }
      };
  merge_optional_vec(forward_per_base_depth, other.forward_per_base_depth, "forward_per_base_depth");
  merge_optional_vec(reverse_per_base_depth, other.reverse_per_base_depth, "reverse_per_base_depth");
  merge_optional_vec(forward_per_base_errors, other.forward_per_base_errors, "forward_per_base_errors");
  merge_optional_vec(reverse_per_base_errors, other.reverse_per_base_errors, "reverse_per_base_errors");
  return *this;
}

}  // namespace xoos::read_collapser
