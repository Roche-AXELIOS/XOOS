#pragma once

#include <unordered_set>

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/float.h>
#include <xoos/util/hash.h>

#include "core/read-collapser-options.h"
#include "core/region.h"
#include "core/thread-context.h"

namespace xoos::read_collapser {

/// Per-thread set of read names (qnames) of primary alignments that were marked as duplicate
/// and have an SA tag. Used to propagate duplicate status to supplementary alignments.
using DuplicateReadNames = std::unordered_set<std::string, util::hash::TransparentStringHash, std::equal_to<>>;

/**
 * Performed optimized duplicate marking on a BAM file. This mode will:
 *  1. Partition the BAM file into super regions
 *  2. Mark positional duplicates according to start and end position and wiggle room
 *  3. Produce multiple BAM files, one per super region
 *  4. Optionally add cluster ID and cluster size to the BAM files
 *  5. Produce clustering metrics
 *  6. Merge the output BAM files into a single sorted BAM file if requested
 *
 * @param options Configuration for duplicate marking
 * @return The number of super regions processed
 */
u32 MarkDuplicate(const ReadCollapserOptions& options);

/**
 * Mark duplicates on the BAM file specified in @p options, and merge the per-thread output BAM files
 * into a single sorted BAM file if specified in @p options.
 *
 * @param options Configuration for duplicate marking
 */
void MarkDuplicateAndMergeOutput(const ReadCollapserOptions& options);

/**
 * Perform duplicate marking on @p super_region using @p context, reading alignments and writing the marked clusters.
 * Use @p super_regions to check whether any reads overlap other super regions.
 *
 * When mark_supplementary_alignments is enabled, the qnames of primary alignments that are marked as duplicate
 * and have an SA tag are inserted into @p duplicate_read_names for later supplementary duplicate marking.
 *
 * Clustering and discard metrics are accumulated into `context.GetMetrics()`, which must be dedicated to a single
 * super-region task so concurrently executing tasks do not race on it.
 */
void DuplicateMarkSuperRegion(ThreadContext& context,
                              const ReadCollapserOptions& options,
                              const SuperRegion& super_region,
                              u32 super_region_id,
                              const vec<SuperRegion>& super_regions,
                              DuplicateReadNames& duplicate_read_names);

/**
 * Determine if the alignment record should be skipped because it overlaps with a prior region.
 * This includes regions within the same super region or the last region of a prior super region.
 *
 * An alignment is owned by the region in which its start position falls. When an alignment starts
 * before the current region (because it overlaps the current region boundary), this function checks
 * whether it overlaps the prior region, in which case the prior region already processed it and the
 * current region should skip it to avoid double-counting across region boundaries.
 *
 * @param record The alignment record to check.
 * @param super_regions The vector of all super regions.
 * @param super_region_id The ID of the current super region.
 * @param subregion_index The index of the current subregion within the super region.
 * @return True if the alignment overlaps with a prior region and should be skipped, false otherwise.
 */
bool ShouldSkipAlignment(const bam1_t* record,
                         const vec<vec<Region>>& super_regions,  // NOSONAR(M23_225) — std::vector is idiomatic here
                         size_t super_region_id,
                         size_t subregion_index);

}  // namespace xoos::read_collapser
