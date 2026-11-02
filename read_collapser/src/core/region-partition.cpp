#include "core/region-partition.h"

#include <algorithm>
#include <iterator>

#include <htslib/hts.h>

#include <xoos/error/error.h>
#include <xoos/io/bed-region.h>
#include <xoos/io/htslib-util/hts-idx-internal.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "core/region.h"

namespace xoos::read_collapser {

/**
 * 64-bit virtual file offset used in BAM index
 * The higher 48 bits represent the byte offset of the BGZF block in the file
 * The lower 16 bits represent the offset within the uncompressed BGZF block
 */
using VOffset = u64;

// Linear index covers genome in 16kb windows
constexpr auto kLinearIndexGenomeWidth = 16'384;
// Virtual file offsets have 16 bits for the uncompressed offset indicating the offset within the BGZF block
constexpr auto kUncompressedFileOffsetWidth = 16U;
// Size of the BGZF EOF block, used to estimate unmapped read data size,
constexpr auto kBgzfEofBlockSize = 28;

namespace {

/**
 * Convert file offset to virtual file offset by shifting left to make space for uncompressed offset.
 *
 * @param file_offset The file offset in bytes. Should be less than 2^48.
 */
VOffset ToVirtualFileOffset(const u64 file_offset) {
  return (file_offset << kUncompressedFileOffsetWidth);
}

/**
 * Find the start virtual file offset of the next 16kb window following the window at index `i` for
 * reference sequence ID `tid`.
 * It uses the binning (hierarchical) index to find the end offset if the current window is the last
 * one for the reference sequence.
 */
VOffset FindNextLinearIndexOffset(const hts_idx_t* const idx,
                                  const s32 tid,
                                  const vec<u64>& linear_offsets,
                                  const s64 i) {
  // If the next window exists in the linear index, return its start offset.
  if (i + 1 < std::ssize(linear_offsets)) {
    return linear_offsets[i + 1];
  }
  // Otherwise, return the maximum mapped offset for this reference sequence.
  return io::BamIndexMaxMappedOffset(idx, tid);
}

/**
 * Calculate the estimated workload for a given range of virtual file offsets in terms of
 * bytes of BGZF compressed data.
 */
u64 CalculateEstimatedWorkload(const VOffset start_offset, const VOffset end_offset) {
  if (end_offset < start_offset) {
    return 0;
  }
  const u64 end_file_offset = end_offset >> kUncompressedFileOffsetWidth;
  const u64 start_file_offset = start_offset >> kUncompressedFileOffsetWidth;
  return (end_file_offset > start_file_offset) ? (end_file_offset - start_file_offset) : 0;
}

}  // namespace

vec<Region> ParseBamLinearIndex(const fs::path& bam_path, const bool include_unmapped) {
  const auto hdr = io::SamHdrRead(io::SamOpen(bam_path, "r").get());
  const auto idx = io::HtsIdxLoad(bam_path, HTS_FMT_BAI);

  vec<Region> intervals;

  VOffset end_offset = 0;
  const auto n_ref = io::SamHdrNRef(hdr.get());
  const auto idx_n_ref = io::BamIndexNumRef(idx.get());
  for (s32 tid = 0; tid < n_ref; ++tid) {
    const auto linear_offsets = io::BamIndexLinearOffsets(idx.get(), tid);
    // If there is no data for this tid, create an interval which contains
    // the entire genome length with no workload for this tid.
    if (tid >= idx_n_ref || linear_offsets.empty()) {
      intervals.emplace_back(
          Region{.tid = tid, .start = 0, .end = io::SamHdrTid2Length(hdr.get(), tid), .estimated_workload = 0});
      continue;
    }

    const auto min_offset = io::BamIndexMinMappedOffset(idx.get(), tid);
    for (s64 i = 0; i < std::ssize(linear_offsets); ++i) {
      // If the start offset is less than the minimum mapped offset for this tid,
      // the start offset is not valid, and we use the minimum mapped offset instead.
      const auto start_offset = std::max(min_offset, linear_offsets[i]);
      end_offset = FindNextLinearIndexOffset(idx.get(), tid, linear_offsets, i);
      intervals.emplace_back(Region{.tid = tid,
                                    .start = i * kLinearIndexGenomeWidth,
                                    .end = (i + 1) * kLinearIndexGenomeWidth,
                                    .estimated_workload = CalculateEstimatedWorkload(start_offset, end_offset)});
    }
  }

  if (include_unmapped) {
    const auto nocoor_start_offset = end_offset;
    // Get the virtual offset corresponding to the end offset of the last block containing unmapped reads (if any)
    const auto file_size = fs::file_size(bam_path);
    // The last byte of the block containing unmapped reads (if any) is before the BGZF EOF block
    // so to get the end offset of unmapped reads, we need to subtract the size of the BGZF EOF block
    const auto nocoor_end_offset =
        ToVirtualFileOffset(file_size > kBgzfEofBlockSize ? file_size - kBgzfEofBlockSize : file_size);
    if (io::BamIndexNumNoCoor(idx.get()) > 0) {
      intervals.emplace_back(
          Region{.tid = HTS_IDX_NOCOOR,
                 .start = 0,
                 .end = 0,
                 .estimated_workload = CalculateEstimatedWorkload(nocoor_start_offset, nocoor_end_offset)});
    }
  }

  return intervals;
}

void EstimateWorkloadForRegions(vec<Region>& regions, const vec<Region>& bam_intervals_with_estimated_workload) {
  // Use the BAM intervals to extrapolate estimated workload for each region
  if (regions.empty()) {
    return;
  }
  // With no BAM intervals there is no overlap to attribute, so every region's workload is reset to
  // zero below (honoring the "existing values are overwritten" contract) rather than left stale.
  auto interval_it = bam_intervals_with_estimated_workload.begin();
  for (auto& [tid, start, end, estimated_workload] : regions) {
    u64 region_workload = 0;

    // Advance interval_it to the first BAM interval that could possibly overlap this interval.
    // We can skip any BAM intervals that end before this region starts.
    while (interval_it != bam_intervals_with_estimated_workload.end() &&
           (interval_it->tid < tid || (interval_it->tid == tid && interval_it->end <= start))) {
      ++interval_it;
    }

    // From interval_it, check all subsequent BAM intervals that could overlap with the current region.
    auto current_check_it = interval_it;
    while (current_check_it < bam_intervals_with_estimated_workload.end() && current_check_it->tid == tid &&
           current_check_it->start < end) {
      const auto& interval = *current_check_it;

      // Calculate the length of the overlap.
      const s64 overlap_start = std::max(start, interval.start);
      const s64 overlap_end = std::min(end, interval.end);
      const s64 overlap_length = overlap_end - overlap_start;

      if (overlap_length > 0) {
        const s64 interval_length = interval.end - interval.start;
        if (interval_length > 0) {
          region_workload += static_cast<u64>((static_cast<f64>(overlap_length) / static_cast<f64>(interval_length)) *
                                              static_cast<f64>(interval.estimated_workload));
        }
      }
      ++current_check_it;
    }
    estimated_workload = region_workload;
  }

  for (const auto& [tid, start_pos, end_pos, bgzf_bytes_covered] : bam_intervals_with_estimated_workload) {
    if (tid == HTS_IDX_NOCOOR) {
      // If the interval is for unmapped reads, we can add it directly.
      regions.emplace_back(tid, start_pos, end_pos, bgzf_bytes_covered);
    }
  }
}

void EstimateWorkloadForRegions(vec<Region>& regions, const fs::path& bam_path, const bool include_unmapped) {
  EstimateWorkloadForRegions(regions, ParseBamLinearIndex(bam_path, include_unmapped));
}

vec<vec<Region>> PartitionRegionsByCustomCostFunction(const vec<Region>& intervals,
                                                      const u64 target_partition_count,
                                                      const RegionPartitionCostFunction auto& cost_function) {
  if (target_partition_count == 0) {
    throw error::Error("Number of partitions must be greater than zero.");
  }
  if (intervals.empty()) {
    return {};
  }
  const u64 target_cost_per_partition = FindOptimalWorkSize(intervals, target_partition_count, cost_function);

  vec<vec<Region>> partitions;

  partitions.emplace_back();
  u64 current_cost = 0;

  for (const auto& interval : intervals) {
    // If the current partition exceeds the optimal size, and we still have partitions to create,
    // start a new partition.
    if (current_cost >= target_cost_per_partition && partitions.size() < target_partition_count) {
      partitions.emplace_back();
      current_cost = 0;
    }
    partitions.back().emplace_back(interval);
    current_cost += cost_function(interval);
  }

  // If we have fewer partitions than requested, fill with empty partitions.
  while (partitions.size() < target_partition_count) {
    partitions.emplace_back();
  }

  return partitions;
}

vec<vec<Region>> PartitionRegionsByWidth(const vec<Region>& intervals, const u64 target_partition_count) {
  return PartitionRegionsByCustomCostFunction(
      intervals, target_partition_count, [](const Region& interval) { return interval.Width(); });
}

vec<vec<Region>> PartitionRegionsByEstimatedWorkload(const vec<Region>& intervals, const u64 target_partition_count) {
  return PartitionRegionsByCustomCostFunction(
      intervals, target_partition_count, [](const Region& interval) { return interval.estimated_workload; });
}

vec<vec<Region>> PartitionRegionsByWidthWithoutCountLimit(const vec<Region>& intervals,
                                                          const u64 target_partition_width) {
  if (target_partition_width == 0) {
    throw error::Error("Target partition width must be greater than zero.");
  }
  if (intervals.empty()) {
    return {};
  }
  vec<vec<Region>> partitions;
  vec<Region> current_partition;
  for (const auto& interval : intervals) {
    if (current_partition.empty()) {
      current_partition.emplace_back(interval);
    } else {
      const auto& first_region = current_partition.front();
      // Current interval is from a different contig or adding it would exceed the target width
      if (interval.tid != first_region.tid ||
          (interval.end - first_region.start) > static_cast<s64>(target_partition_width)) {
        partitions.emplace_back(current_partition);
        current_partition.clear();
      }
      current_partition.emplace_back(interval);
    }
  }
  if (!current_partition.empty()) {
    partitions.emplace_back(current_partition);
  }
  return partitions;
}

}  // namespace xoos::read_collapser
