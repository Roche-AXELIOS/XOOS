#include "core/region-partition.h"

#include <numeric>

#include <htslib/hts.h>

#include <xoos/io/bed-region.h>
#include <xoos/io/htslib-util/htslib-internal.h>
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
constexpr auto kUncompressedFileOffsetWidth = 16;
// Size of the BGZF EOF block, used to estimate unmapped read data size,
constexpr auto kBgzfEofBlockSize = 28;

/**
 * Find the maximum virtual file offset for the given reference sequence ID (tid) in the BAM index.
 * This is the virtual file offset of the end of reads mapped to this reference sequence.
 *
 * @param idx Pointer to the BAM index structure.
 * @param tid Reference sequence ID.
 * @return Maximum virtual file offset for the given tid, or 0 if tid is invalid or has no data.
 */
static u64 FindMaxMappedOffset(const io::htslib::hts_idx_t* const idx, const s32 tid) {
  u64 max_offset = 0;
  if (tid < 0 || tid >= idx->n) {
    return 0;
  }
  const auto* const bidx = idx->bidx[tid];
  if (bidx == nullptr) {
    return 0;
  }
  // If the binning (hierarchical) index has the metadata bin, read the max offset directly from there
  const auto meta_bin = kh_get(bin, bidx, META_BIN(idx));
  if (meta_bin != kh_end(bidx) && kh_val(bidx, meta_bin).n > 0) {
    return kh_val(bidx, meta_bin).list[0].v;
  }
  // Otherwise, manually search through all bins and the linear index to find the maximum offset
  // First find the maximum offset in the hierarchical index for this reference sequence
  for (khint_t k = 0; k < kh_end(bidx); ++k) {
    if (!kh_exist(bidx, k)) {
      continue;
    }
    const auto* const bin = &kh_val(bidx, k);
    if (bin->list == nullptr || bin->n <= 0) {
      continue;
    }
    for (s32 j = 0; j < bin->n; ++j) {
      max_offset = std::max(max_offset, bin->list[j].v);
    }
  }
  // Now check the linear index for this reference sequence
  if (idx->lidx != nullptr) {
    const auto& lidx = idx->lidx[tid];
    for (hts_pos_t j = 0; j < lidx.n; ++j) {
      max_offset = std::max(max_offset, lidx.offset[j]);
    }
  }
  return max_offset;
}

/**
 * Find the minimum virtual file offset for the given reference sequence ID (tid) in the BAM index.
 * This is the virtual file offset of the start of reads mapped to this reference sequence.
 *
 * @param idx Pointer to the BAM index structure.
 * @param tid Reference sequence ID.
 * @return Minimum virtual file offset for the given tid, or 0 if tid is invalid or has no data.
 */
static u64 FindMinMappedOffset(const io::htslib::hts_idx_t* const idx, const s32 tid) {
  u64 min_offset = std::numeric_limits<u64>::max();
  if (tid < 0 || tid >= idx->n) {
    return 0;
  }
  const auto* const bidx = idx->bidx[tid];
  if (bidx == nullptr) {
    return 0;
  }
  // If the binning (hierarchical) index has the metadata bin, read the min offset directly from there
  const auto meta_bin = kh_get(bin, bidx, META_BIN(idx));
  if (meta_bin != kh_end(bidx) && kh_val(bidx, meta_bin).n > 0) {
    return kh_val(bidx, meta_bin).list[0].u;
  }
  // Otherwise, manually search through all bins and the linear index to find the minimum offset
  // First find the maximum offset in the hierarchical index for this reference sequence
  for (khint_t k = 0; k < kh_end(bidx); ++k) {
    if (!kh_exist(bidx, k)) {
      continue;
    }
    const auto* const bin = &kh_val(bidx, k);
    if (bin->list == nullptr || bin->n <= 0) {
      continue;
    }
    for (s32 j = 0; j < bin->n; ++j) {
      min_offset = std::min(min_offset, bin->list[j].u);
    }
  }
  // Now check the linear index for this reference sequence
  if (idx->lidx != nullptr) {
    const auto& lidx = idx->lidx[tid];
    for (hts_pos_t j = 0; j < lidx.n; ++j) {
      min_offset = std::min(min_offset, lidx.offset[j]);
    }
  }
  return min_offset == std::numeric_limits<u64>::max() ? 0 : min_offset;
}

/**
 * Convert file offset to virtual file offset by shifting left to make space for uncompressed offset.
 *
 * @param file_offset The file offset in bytes. Should be less than 2^48.
 */
static VOffset ToVirtualFileOffset(const u64 file_offset) {
  return (file_offset << kUncompressedFileOffsetWidth);
}

/**
 * Find the start virtual file offset of the next 16kb window following the window at index `i` for
 * reference sequence ID `tid`.
 * It uses the binning (hierarchical) index to find the end offset if the current window is the last
 * one for the reference sequence.
 */
static VOffset FindNextLinearIndexOffset(const io::htslib::hts_idx_t* const idx, const s32 tid, const s64 i) {
  // If the next window exists in the linear index, return its start offset.
  if (i + 1 < idx->lidx[tid].n) {
    return idx->lidx[tid].offset[i + 1];
  }
  // Otherwise, return the maximum mapped offset for this reference sequence.
  return FindMaxMappedOffset(idx, tid);
}

/**
 * Calculate the estimated workload for a given range of virtual file offsets in terms of
 * bytes of BGZF compressed data.
 */
static u64 CalculateEstimatedWorkload(const VOffset start_offset, const VOffset end_offset) {
  if (end_offset < start_offset) {
    return 0;
  }
  const u64 end_file_offset = end_offset >> kUncompressedFileOffsetWidth;
  const u64 start_file_offset = start_offset >> kUncompressedFileOffsetWidth;
  return (end_file_offset > start_file_offset) ? (end_file_offset - start_file_offset) : 0;
}

/**
 * Find the optimal work size per partition for the given number of partitions and cost function.
 * This provides a better estimate than the naive approach of total_cost / target_partition_count
 * when the cost distribution is uneven.
 *
 * `target_partition_count` must be greater than zero.
 */
static u64 FindOptimalWorkSize(const vec<Region>& regions,
                               const u64 target_partition_count,
                               const RegionPartitionCostFunction auto& cost_function) {
  // Calculate the total work for all intervals
  u64 total_work = 0;
  for (const auto& interval : regions) {
    total_work += cost_function(interval);
  }

  if (total_work == 0 || target_partition_count == 0) {
    return 0;
  }

  // The goal is to find the largest possible work size for a partition
  // that still yields `target_partition_count` non-empty partitions.
  u64 low = 0;
  u64 high = total_work / target_partition_count + 1;
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

vec<Region> ParseBamLinearIndex(const fs::path& bam_path, const bool include_unmapped) {
  const auto hdr = io::SamHdrRead(io::SamOpen(bam_path, "r").get());
  const auto idx = io::htslib::HtsIdxLoad(bam_path, HTS_FMT_BAI);

  vec<Region> intervals;

  VOffset end_offset = 0;
  const auto n_ref = io::SamHdrNRef(hdr.get());
  for (s32 tid = 0; tid < n_ref; ++tid) {
    // If there is no data for this tid, create an interval which contains
    // the entire genome length with no workload for this tid.
    if (tid >= idx->n || idx->lidx[tid].n == 0) {
      intervals.emplace_back(
          Region{.tid = tid, .start = 0, .end = io::SamHdrTid2Length(hdr.get(), tid), .estimated_workload = 0});
      continue;
    }

    const auto& lidx = idx->lidx[tid];
    const auto min_offset = FindMinMappedOffset(idx.get(), tid);
    for (s64 i = 0; i < lidx.n; ++i) {
      // If the start offset is less than the minimum mapped offset for this tid,
      // the start offset is not valid, and we use the minimum mapped offset instead.
      const auto start_offset = std::max(min_offset, lidx.offset[i]);
      end_offset = FindNextLinearIndexOffset(idx.get(), tid, i);
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
    if (idx->n_no_coor > 0) {
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
  if (regions.empty() || bam_intervals_with_estimated_workload.empty()) {
    return;
  }
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
  return EstimateWorkloadForRegions(regions, ParseBamLinearIndex(bam_path, include_unmapped));
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
