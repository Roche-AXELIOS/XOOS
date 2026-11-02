#pragma once

#include <concepts>
#include <filesystem>
#include <numeric>

#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "core/region.h"

namespace xoos::read_collapser {

template <typename Func>
concept RegionPartitionCostFunction =
    std::regular_invocable<Func, const Region&> && std::same_as<std::invoke_result_t<Func, const Region&>, u64>;

/**
 * Find the optimal work size per partition for the given number of partitions and cost function.
 * This provides a better estimate than the naive approach of total_cost / target_partition_count
 * when the cost distribution is uneven.
 *
 * Callers are expected to pass a `target_partition_count` greater than zero. As a defensive guard
 * (the count is used as a divisor), 0 is returned when `target_partition_count` is zero or when the
 * total work across all regions is zero.
 */
u64 FindOptimalWorkSize(const vec<Region>& regions,
                        const u64 target_partition_count,
                        const RegionPartitionCostFunction auto& cost_function) {
  // Calculate the total work for all intervals
  u64 total_work = 0;
  for (const auto& interval : regions) {
    total_work += cost_function(interval);
  }

  if ((total_work == 0) || (target_partition_count == 0)) {
    return 0;
  }

  // The goal is to find the largest possible work size for a partition
  // that still yields `target_partition_count` non-empty partitions.
  u64 low = 0;
  u64 high = (total_work / target_partition_count) + 1;
  u64 optimal_work_size = high;

  while (low <= high) {
    const u64 mid = std::midpoint(low, high);
    // Avoid division by zero or infinite loops
    if (mid == 0) {
      low = 1;
      continue;
    }

    // Simulate partitioning with the current `mid` as the target work size
    u32 partitions_created = 0;
    u64 current_partition_work = 0;
    for (const auto& interval : regions) {
      if (current_partition_work >= mid) {
        ++partitions_created;
        current_partition_work = 0;
      }
      current_partition_work += cost_function(interval);
    }
    // Account for the last partition
    if (current_partition_work > 0) {
      ++partitions_created;
    }

    if (partitions_created >= target_partition_count) {
      // This work size is feasible, try for a larger one to make partitions bigger
      optimal_work_size = mid;
      low = mid + 1;
    } else {
      // This work size is too large, it creates too few partitions.
      high = mid - 1;
    }
  }

  return optimal_work_size;
}

/**
 * Parse the BAM linear index to extract 16kb intervals with estimated workload.
 *
 * @param bam_path Path to the BAM file.
 * @param include_unmapped Whether to include unmapped reads as a separate interval.
 * @return A vector of intervals with estimated workload.
 */
vec<Region> ParseBamLinearIndex(const fs::path& bam_path, bool include_unmapped);

/**
 * Estimate and assign workload for the given regions based on the BAM index.
 *
 * @param regions Vector of intervals for which to estimate workload. Workload will be updated in place.
 * Any existing `estimated_workload` values will be overwritten. Should not contain intervals for unmapped reads
 * as those will be added separately.
 * @param bam_intervals_with_estimated_workload Vector of intervals with estimated workload from the BAM index.
 */
void EstimateWorkloadForRegions(vec<Region>& regions, const vec<Region>& bam_intervals_with_estimated_workload);

/**
 * Estimate and assign workload for the given regions based on the BAM index.
 *
 * @param regions Vector of intervals for which to estimate workload. Workload will be updated in place.
 * Any existing `estimated_workload` values will be overwritten. Should not contain intervals for unmapped reads
 * as those will be added separately.
 * @param bam_path Path to the BAM file. Used to read the BAM index for workload estimation.
 * @param include_unmapped Whether to include unmapped reads as a separate interval.
 */
void EstimateWorkloadForRegions(vec<Region>& regions, const fs::path& bam_path, bool include_unmapped);

/**
 * Partition the given intervals into at most `target_partition_count` partitions,
 * attempting to balance the total cost of intervals in each partition with respect to
 * the provided `cost_function`.
 */
vec<SuperRegion> PartitionRegionsByCustomCostFunction(const vec<Region>& intervals,
                                                      u64 target_partition_count,
                                                      const RegionPartitionCostFunction auto& cost_function);

/**
 * Partition the given intervals into at most `target_partition_count` partitions,
 * attempting to balance the total width of intervals in each partition.
 */
vec<SuperRegion> PartitionRegionsByWidth(const vec<Region>& intervals, u64 target_partition_count);

/**
 * Partition the given intervals into at most `target_partition_count` partitions,
 * attempting to balance the total estimated workload of intervals in each partition.
 */
vec<SuperRegion> PartitionRegionsByEstimatedWorkload(const vec<Region>& intervals, u64 target_partition_count);

/**
 * Partition the given intervals into partitions such that each partition has total width
 * approximately equal to `target_partition_width`. The number of partitions is not limited.
 */
vec<SuperRegion> PartitionRegionsByWidthWithoutCountLimit(const vec<Region>& intervals, u64 target_partition_width);

}  // namespace xoos::read_collapser
