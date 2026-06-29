#include "core/region.h"

#include <algorithm>
#include <cstddef>
#include <utility>

#include <csv.hpp>

#include <xoos/error/error.h>
#include <xoos/io/bed-region.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/io/malloc-ptr.h>
#include <xoos/types/float.h>

#include "core/region-partition.h"

namespace xoos::read_collapser {

u64 Region::Width() const {
  return static_cast<u64>(end > start ? end - start : 0);
}

bool Region::Overlaps(const Region& other) const {
  if (tid != other.tid) {
    return false;
  }
  return !(end <= other.start || other.end <= start);
}

vec<Region> ReadBamTargetRegions(const fs::path& bam_filename) {
  const io::HtsFilePtr bam = io::HtsOpen(bam_filename, "rb");
  const io::SamHdrPtr header = io::SamHdrRead(bam.get());

  vec<Region> regions = vec<Region>();
  for (s32 tid = 0; tid < header->n_targets; ++tid) {
    regions.emplace_back(tid, 0, header->target_len[tid]);
  }
  return regions;
}

vec<Region> ReadBedTargetRegions(const fs::path& bed_filename, const fs::path& bam_filename, const s32 padding) {
  vec<io::BedRegion> bed_regions = io::ParseBedFile(bed_filename);
  if (bed_regions.empty()) {
    throw error::Error("BED file '{}', contains no regions", bed_filename);
  }
  // Ensure that the BED regions are sorted and non-overlapping
  bed_regions = io::MergeOverlappingRegions(std::move(bed_regions));
  const vec<io::BedRegion> padded_regions = io::PadBedRegions(bed_regions, padding);

  vec<Region> regions;
  regions.reserve(padded_regions.size());

  const io::HtsFilePtr bam = io::HtsOpen(bam_filename, "rb");
  const io::SamHdrPtr header = io::SamHdrRead(bam.get());

  for (const auto& padded_region : padded_regions) {
    regions.emplace_back(
        io::SamHdrName2Tid(header.get(), padded_region.chromosome), padded_region.start, padded_region.end);
  }
  // Although io::MergeOverlappingRegions sort the regions, they are sorted by chromosome name
  // which may not correspond to the tid order in the BAM file.
  std::ranges::sort(regions);
  return regions;
}

vec<Region> SplitLargeRegions(const vec<Region>& regions, const u64 region_size) {
  vec<Region> result_regions;
  for (const auto& region : regions) {
    // Determine the full length of the target region and the number of regions it will be partitioned into
    const u64 target_len = region.Width();
    const u64 target_region_count = (target_len + region_size - 1) / region_size;

    // Generate a region for each partitioned region, the last region may be smaller than region_size
    // if the target region is not evenly divisible by region_size
    for (s64 j = 0; std::cmp_less(j, target_region_count); ++j) {
      const s64 beg = region.start + j * ToSigned(region_size);
      const s64 end = std::min(region.end, region.start + (j + 1) * ToSigned(region_size));
      result_regions.emplace_back(region.tid, beg, end);
    }
  }
  return result_regions;
}

vec<SuperRegion> DetermineSuperRegions(const ReadCollapserOptions& options) {
  const auto regions_initial = options.bed_input
                                   ? ReadBedTargetRegions(options.bed_input.value(), options.bam_input, options.padding)
                                   : ReadBamTargetRegions(options.bam_input);
  const auto regions = SplitLargeRegions(regions_initial, options.region_size);
  auto partitioned_regions = PartitionRegionsByWidthWithoutCountLimit(regions, options.region_size);
  // Manually append region for unmapped reads
  partitioned_regions.emplace_back(vec<Region>{Region{HTS_IDX_NOCOOR, 0, 0}});
  return partitioned_regions;
}

vec<SuperRegion> DetermineSuperRegions(const ReadCollapserOptions& options, const size_t max_super_region_count) {
  const auto regions_initial = options.bed_input
                                   ? ReadBedTargetRegions(options.bed_input.value(), options.bam_input, options.padding)
                                   : ReadBamTargetRegions(options.bam_input);
  auto regions = SplitLargeRegions(regions_initial, options.region_size);
  // Estimate workload and append unmapped region (done by the function)
  EstimateWorkloadForRegions(regions, options.bam_input, true);
  return PartitionRegionsByEstimatedWorkload(regions, max_super_region_count);
}

}  // namespace xoos::read_collapser
