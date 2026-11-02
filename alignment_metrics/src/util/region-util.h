#pragma once

#include <optional>

#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "io/bed-reader.h"

namespace xoos::alignment_metrics {

/**
 * Given a list of possibly overlapping regions grouped by chromosome, merge overlapping regions
 * into non-overlapping regions. Region labels are ignored (the result carries no labels); use
 * MergeRegionsByLabel when label stratification is required.
 *
 * The function needs to sort the regions internally so it needs to take ownership of the input
 * either moved from an rvalue or copied from an lvalue.
 *
 * @param regions Possibly overlapping regions grouped by chromosome.
 * Pass an rvalue to avoid unnecessary copies.
 */
RegionsByChromosome MergeOverlappingRegions(RegionsByChromosome regions);

/**
 * Partition the given regions into smaller regions of approximately `region_size` length.
 * If `reference_path` is provided, we will use the reference sequence to find homopolymers
 * and ensure that homopolymers are not split across regions if `hp_aware` is true.
 *
 * If `hp_aware` is true, we will find homopolymers in the reference sequence and adjust the
 * region boundaries to avoid splitting homopolymers. This may result in regions that are
 * larger than `region_size`.
 *
 * If `reference_path` is not provided, `hp_aware` will be ignored and regions will be
 * partitioned solely based on `region_size`.
 *
 * @param regions Regions to partition, grouped by chromosome
 * @param region_size Desired size of each partitioned region
 * @param reference_path Optional path to the reference FASTA file
 * @param hp_aware If true, adjust region boundaries to avoid splitting homopolymers
 * @return Partitioned regions grouped by chromosome
 */
RegionsByChromosome PartitionRegions(const RegionsByChromosome& regions,
                                     u32 region_size,
                                     const std::optional<fs::path>& reference_path = std::nullopt,
                                     bool hp_aware = false);

}  // namespace xoos::alignment_metrics
