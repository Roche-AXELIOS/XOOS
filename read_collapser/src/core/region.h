#pragma once

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "core/read-collapser-options.h"

namespace xoos::read_collapser {

/**
 * A genomic region defined by a reference sequence ID, start position, end position, and estimated workload
 * associated with the region.
 *
 * The start position is inclusive, and the end position is exclusive.
 * Both positions are 0-based.
 */
struct Region {
  // Reference sequence ID (tid) as per BAM format.
  s32 tid{};
  // Start position (0-based inclusive).
  s64 start{};
  // End position (0-based exclusive).
  s64 end{};
  // Estimated workload associated with this interval, e.g., in bytes of BGZF compressed data.
  // Can be set to 0 if not applicable.
  u64 estimated_workload{};
  u64 Width() const;
  bool Overlaps(const Region& other) const;
  auto operator<=>(const Region& other) const = default;
};

/**
 * Read target regions from the BAM file header. Each target found in the header
 * is represented as a `Region` with the start position set to 0 and the end
 * position set to the length of the target.
 *
 * Regions are sorted by reference sequence ID (tid) and position.
 */
vec<Region> ReadBamTargetRegions(const fs::path& bam_filename);

/**
 * Read target regions from a BED file, optionally padding each region by a specified number of bases.
 * The chromosome names in the BED file are mapped to reference sequence IDs (tids) using the BAM file header.
 *
 * Regions are sorted by reference sequence ID (tid) and position.
 */
vec<Region> ReadBedTargetRegions(const fs::path& bed_filename, const fs::path& bam_filename, s32 padding);

/**
 * Split large regions into smaller regions of a specified maximum size.
 * The regions are split into fixed-size segments, except for the last segment,
 * which may be limited by reaching the end of the original region.
 */
vec<Region> SplitLargeRegions(const vec<Region>& regions, u64 region_size);

/**
 * A super region is a collection of regions that are grouped together for processing.
 * Partitioning creates a vector of super regions from a vector of individual regions.
 */
using SuperRegion = vec<Region>;

/**
 * Create partitioned super regions based on the provided read collapser options.
 * The number of super regions (partitions) is not limited, but the width of each partition
 * is approximately equal to the specified region size in the options.
 *
 * Each SuperRegion contains one or more Regions from the same contig.
 * A SuperRegion containing a single region with tid set to HTS_IDX_NOCOOR is added
 * to the end of the returned vector for processing unmapped reads.
 *
 * @param options The ReadCollapserOptions containing region size information. `bam_input`
 * is required for mapping chromosome names in BED files to reference sequence IDs (tids).
 */
vec<SuperRegion> DetermineSuperRegions(const ReadCollapserOptions& options);

/**
 * Create partitioned super regions based on the provided read collapser options,
 * limiting the number of super regions to `max_super_region_count` and attempting
 * to balance the width of each partition.
 *
 * Each SuperRegion may contain Regions from multiple contigs in order to limit
 * the total number of super regions created.
 *
 * @param options The ReadCollapserOptions containing region size information. `bam_input`
 * is required for mapping chromosome names in BED files to reference sequence IDs (tids)
 * as well as estimating workloads from the BAM index.
 * @param max_super_region_count The maximum number of super regions to create.
 */
vec<SuperRegion> DetermineSuperRegions(const ReadCollapserOptions& options, size_t max_super_region_count);

}  // namespace xoos::read_collapser
