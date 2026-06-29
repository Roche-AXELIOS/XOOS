#pragma once

#include <concepts>
#include <filesystem>

#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "core/region.h"

namespace xoos::read_collapser {

template <typename Func>
concept RegionPartitionCostFunction =
    std::regular_invocable<Func, const Region&> && std::same_as<std::invoke_result_t<Func, const Region&>, u64>;

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
