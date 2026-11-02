#pragma once

#include <algorithm>
#include <array>
#include <concepts>
#include <optional>
#include <queue>
#include <string_view>
#include <utility>

#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "core/interval.h"
#include "core/region-labels-set.h"
#include "io/bed-reader.h"

namespace xoos::alignment_metrics {

namespace detail {

/**
 * Represents the end of an active interval along with its associated labels.
 * Used to keep track of active intervals in a priority queue during merging of labeled regions.
 */
struct ActiveIntervalEnd {
  s64 end;
  u16 labels;
};

/**
 * Update the active labels and mask when we encounter an endpoint (start or end) of an interval
 * for the sweep-line merge algorithm.
 *
 * @param labels The labels of the interval being processed.
 * @param starting True if the endpoint corresponds to the start of an interval, false if it corresponds to the end.
 * @param mask The current active label mask before processing this endpoint.
 * @param active_counts_per_label Array tracking the count of active intervals for each label. The array is updated
 * in place. If the count for a label drops to zero, the corresponding bit is cleared from the mask indicating the
 * from zero, the corresponding bit is set in the mask indicating the running interval now has that label.
 * @return The updated active label mask after processing this endpoint.
 */
inline u16 UpdateActiveLabels(const u16 labels,
                              const bool starting,
                              const u16 mask,
                              std::array<size_t, kMaxLabels>& active_counts_per_label) {
  u16 new_mask = mask;
  for (size_t bit = 0; bit < kMaxLabels; ++bit) {
    auto& active_count = active_counts_per_label.at(bit);
    if ((labels & (1U << bit)) == 0) {
      continue;
    }
    if (starting) {
      ++active_count;
    } else {
      --active_count;
    }
    const auto bit_flag = static_cast<u16>(1U << bit);
    new_mask = active_counts_per_label.at(bit) > 0 ? static_cast<u16>(new_mask | bit_flag)
                                                   : static_cast<u16>(new_mask & static_cast<u16>(~bit_flag));
  }
  return new_mask;
}

}  // namespace detail

template <typename F>
concept MergedRegionCallback = std::invocable<F&, std::string_view, Interval&>;

/**
 * Merge labeled intervals on a given chromosome using a sweep-line algorithm.
 * The algorithm maintains a min-heap of active interval ends and a count of active intervals per label
 * to efficiently determine the merged intervals with their active labels.
 *
 * @param chromosome The name of the chromosome for which intervals are being merged.
 * @param intervals The list of intervals on the chromosome to be merged.
 * @param callback The callback function to be invoked for each merged interval.
 */
void MergeLabeledIntervals(const std::string_view chromosome,
                           vec<Interval>& intervals,
                           MergedRegionCallback auto&& callback) {
  std::ranges::sort(intervals);

  const auto earlier_end = [](const detail::ActiveIntervalEnd& left, const detail::ActiveIntervalEnd& right) {
    return left.end > right.end;
  };
  std::priority_queue<detail::ActiveIntervalEnd, vec<detail::ActiveIntervalEnd>, decltype(earlier_end)> active_ends(
      earlier_end);
  std::array<size_t, kMaxLabels> active_counts_per_label{};
  u16 mask = 0;
  std::optional<Interval> pending_merged_interval;
  size_t interval_index = 0;
  s64 segment_start = intervals.front().start;

  while (interval_index < intervals.size() || !active_ends.empty()) {
    const s64 next_start =
        interval_index < intervals.size() ? intervals.at(interval_index).start : std::numeric_limits<s64>::max();
    const s64 next_end = active_ends.empty() ? std::numeric_limits<s64>::max() : active_ends.top().end;
    const s64 position = std::min(next_start, next_end);

    if (!active_ends.empty() && position > segment_start) {
      const bool extends_pending_merged_interval = pending_merged_interval.has_value() &&
                                                   pending_merged_interval->end == segment_start &&
                                                   pending_merged_interval->labels == mask;
      if (extends_pending_merged_interval) {
        pending_merged_interval->end = position;
      } else {
        const std::optional<Interval> completed =
            std::exchange(pending_merged_interval, Interval{segment_start, position, mask});
        if (completed.has_value()) {
          callback(chromosome, completed.value());
        }
      }
    }

    while (!active_ends.empty() && active_ends.top().end == position) {
      mask = detail::UpdateActiveLabels(active_ends.top().labels, false, mask, active_counts_per_label);
      active_ends.pop();
    }

    while (interval_index < intervals.size() && intervals.at(interval_index).start == position) {
      const auto& interval = intervals.at(interval_index);
      mask = detail::UpdateActiveLabels(interval.labels, true, mask, active_counts_per_label);
      active_ends.push(detail::ActiveIntervalEnd{.end = interval.end, .labels = interval.labels});
      ++interval_index;
    }
    segment_start = position;
  }
  if (pending_merged_interval.has_value()) {
    callback(chromosome, pending_merged_interval.value());
  }
}

/**
 * Merge labeled regions and invoke `callback` for each result in chromosome and position order.
 * This overload does not retain merged intervals in memory.
 *
 * The function needs to sort the regions internally so it needs to take ownership of the input
 * either moved from an rvalue or copied from an lvalue.
 *
 * @param regions Possibly overlapping regions grouped by chromosome; each interval's `labels` field
 * is a bit vector of the labels attached to it. Pass an rvalue to avoid unnecessary copies.
 * @param callback Function to invoke for each merged region, in chromosome and position order.
 */
template <MergedRegionCallback F>
void MergeRegionsByLabel(RegionsByChromosome regions, F&& callback) {
  for (auto& [chromosome, intervals] : std::move(regions)) {
    if (!intervals.empty()) {
      MergeLabeledIntervals(chromosome, intervals, std::forward<F>(callback));
    }
  }
}

/**
 * Merge possibly overlapping labeled regions into a set of disjoint, non-overlapping intervals where
 * each output interval carries the union (bitwise-OR) of the labels of every input region covering
 * it.
 *
 * This is a sweep-line merge: at every position where the set of covering regions changes, a new
 * output interval begins. Overlapping regions with different labels are therefore split at their
 * boundaries so that each label class can be attributed independently. For example, region A over
 * [1,10) and region B over [5,20) produce [1,5)=A, [5,10)=A|B, [10,20)=B. Adjacent output intervals
 * that would carry an identical label set are coalesced. Input regions with no labels (bit vector 0)
 * are still merged by coverage, yielding unlabeled output intervals, so this function is equivalent
 * to MergeOverlappingRegions when no labels are present.
 *
 * The function needs to sort the regions internally so it needs to take ownership of the input
 * either moved from an rvalue or copied from an lvalue.
 *
 * @param regions Possibly overlapping regions grouped by chromosome; each interval's `labels` field
 *        is a bit vector of the labels attached to it. Pass an rvalue to avoid unnecessary copies.
 * @return Disjoint intervals grouped by chromosome, each carrying the union label bit vector.
 */
RegionsByChromosome MergeRegionsByLabel(RegionsByChromosome regions);

}  // namespace xoos::alignment_metrics
