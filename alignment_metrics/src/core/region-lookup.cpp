#include "core/region-lookup.h"

#include <cstddef>
#include <numeric>
#include <optional>
#include <utility>

#include <xoos/error/error.h>
#include <xoos/io/fasta-reader.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "core/interval.h"
#include "core/super-region.h"
#include "io/alignment-reader.h"
#include "io/bed-reader.h"
#include "util/region-util.h"

namespace xoos::alignment_metrics {

/**
 * @brief Constructs a RegionLookupTable from BAM and optional BED files.
 *
 * This constructor builds a lookup table for genomic regions of interest by:
 * 1. Loading regions from a BED file (if provided) or extracting them from the BAM file header
 * 2. Merging overlapping regions to avoid redundant processing
 * 3. Optionally partitioning large regions into smaller chunks for memory efficiency
 * 4. Querying the BAM file to identify regions with/without coverage
 * 5. Grouping nearby regions into "super regions" for efficient parallel processing
 *
 * @param bam_path Path to the BAM file to process
 * @param bed_path Optional path to BED file specifying regions of interest
 * @param reference_path Optional path to reference FASTA (required for HP-aware partitioning)
 * @param region_size Maximum size for partitioned regions and super region grouping distance
 * @param hp_aware If true, partition boundaries avoid splitting homopolymers
 * @param disable_region_partitioning If true, skip partitioning and use merged regions as-is
 */
RegionLookupTable::RegionLookupTable(const fs::path& bam_path,
                                     const std::optional<fs::path>& bed_path,
                                     const std::optional<fs::path>& reference_path,
                                     const u32 region_size,
                                     const bool hp_aware,
                                     const bool disable_region_partitioning) {
  auto regions_by_chromosome = bed_path.has_value() && fs::exists(bed_path.value())
                                   ? GetRegionsFromBed(bed_path.value(), bam_path)
                                   : GetRegionsFromBam(bam_path);
  auto merged_regions = MergeOverlappingRegions(regions_by_chromosome);
  RegionsByChromosome partitioned_regions;
  if (!disable_region_partitioning) {
    partitioned_regions = PartitionRegions(merged_regions, region_size, reference_path, hp_aware);
  } else {
    // If region partitioning is disabled, we use the merged regions as is
    // and because merged regions are no longer used beyond this point, we can
    // safely move it to partitioned regions to avoid unnecessary copying.
    partitioned_regions = std::move(merged_regions);
  }
  RegionsByChromosome regions_with_coverage_by_chr;
  // Move the partitioned regions into the lookup table's internal data structure
  // We query the BAM file for chromosomes with no reads mapped to them. For these chromosomes,
  // we add them to `_regions_without_coverage_by_chr` to avoid processing them in the metrics calculation phase.
  const AlignmentReader bam_reader = OpenAlignmentFile(bam_path);
  for (const auto& [chromosome, regions] : partitioned_regions) {
    // Query the BAM using an iterator to determine if there are any reads mapped to this chromosome
    const auto alignment_itr = io::SamItrQueryS(bam_reader.idx.get(), bam_reader.header.get(), chromosome);
    if (alignment_itr.get() == nullptr) {
      _regions_without_coverage_by_chr[chromosome] = regions;
    } else {
      const io::Bam1Ptr bam1_ptr(bam_init1());
      if (sam_itr_next(bam_reader.bam.get(), alignment_itr.get(), bam1_ptr.get()) < 0) {
        _regions_without_coverage_by_chr[chromosome] = regions;
      } else {
        regions_with_coverage_by_chr[chromosome] = regions;
      }
      _total_region_count += regions.size();
    }
  }

  // Now we have populated both `regions_with_coverage_by_chr` and `_regions_without_coverage_by_chr` based on the BED
  // file and the BAM file, we can proceed to create super regions for the regions with coverage. A super region is a
  // group of regions that are processed together. These are regions that are close to each other in the genome and can
  // be processed in a single pass to improve performance.
  vec<SuperRegion> super_regions;
  for (const auto& [chromosome, regions] : regions_with_coverage_by_chr) {
    vec<Interval> sub_regions;
    // If the current super region has no sub-regions, add this region and move on to the next region
    for (const auto& region : regions) {
      if (!sub_regions.empty()) {
        // If this region is too far away from the first region in the super region,
        // then complete this super region and start a new one.
        const auto first_region = sub_regions.front();
        if (std::cmp_greater_equal(region.start - first_region.start, region_size)) {
          // We complete the current super region which contains the sub-regions
          // from `regions[i - sub_regions.size()]` to `regions[i - 1]`.
          super_regions.emplace_back(CreateSuperRegion(chromosome, sub_regions));
          sub_regions.clear();
        }
      }
      sub_regions.emplace_back(region);
    }
    // Handle the remaining super region if there is one
    if (!sub_regions.empty()) {
      super_regions.emplace_back(CreateSuperRegion(chromosome, sub_regions));
    }
  }
  _super_regions = std::move(super_regions);
}

RegionLookupTable::RegionLookupTable(const vec<SuperRegion>& super_regions) : _super_regions(super_regions) {
  _total_region_count =
      std::accumulate(_super_regions.begin(), _super_regions.end(), 0u, [](const size_t sum, const SuperRegion& sr) {
        return sum + sr.subregions.size();
      });
}

const vec<SuperRegion>& RegionLookupTable::GetSuperRegions() const {
  return _super_regions;
}

const RegionsByChromosome& RegionLookupTable::GetRegionsWithoutCoverage() const {
  return _regions_without_coverage_by_chr;
}

// Validates that the trimmed alignment belongs in some region in the specified super region
// This is used as the bed regions will return the original alignments, but the trimmed alignments may not be part of
// the bam
bool RegionLookupTable::IsAlignmentInAnyRegion(const s64 rpos,
                                               const s64 reference_length,
                                               const size_t super_region_index,
                                               const size_t sub_region_index) const {
  const auto& super_region = _super_regions.at(super_region_index);
  // let's check the starting one first for efficiency
  const auto current_sub_region = super_region.subregions.at(sub_region_index);
  // if the sub_region overlaps with the alignment, we return it
  const Interval current_alignment_span{rpos, rpos + reference_length};
  if (current_sub_region.Overlaps(current_alignment_span)) {
    return true;
  }

  return std::ranges::any_of(super_region.subregions, [&current_alignment_span](const auto& sub_region) {
    return sub_region.Overlaps(current_alignment_span);
  });
}

bool RegionLookupTable::IsSafeToCountAlignment(const std::string_view chromosome,
                                               const Interval& alignment_span,
                                               const size_t super_region_index,
                                               const size_t sub_region_index) const {
  if (!alignment_span.Overlaps(
          {_super_regions.at(super_region_index).start, _super_regions.at(super_region_index).end})) {
    // If the alignment does not overlap with the current super region at all, it is not safe to count
    return false;
  }
  // If the alignment start position lies within the current subregion, then it is safe to count the alignment
  const auto& super_region = _super_regions.at(super_region_index);
  const auto& current_sub_region = super_region.subregions.at(sub_region_index);
  if (current_sub_region.Contains(alignment_span.start)) {
    return true;
  }
  // If the interval overlaps with a previous super region, we need to check the previous super region for overlapping
  // regions
  if (super_region_index > 0) {
    const auto& previous_super_region = _super_regions.at(super_region_index - 1);
    if (previous_super_region.chromosome == chromosome && previous_super_region.end > alignment_span.start) {
      for (const auto& sub_region : previous_super_region.subregions) {
        if (sub_region.Overlaps(alignment_span)) {
          // We have found a subregion that overlaps with the alignment in a previous super region
          // so we don't count this alignment for the current subregion to avoid double counting
          return false;
        }
      }
    }
  }
  // No sub-region in prior super region overlaps. It is safe to count this alignment for the current region.
  return true;
}

size_t RegionLookupTable::TotalRegionCount() const {
  return _total_region_count;
}

size_t RegionLookupTable::SuperRegionCount() const {
  return _super_regions.size();
}

}  // namespace xoos::alignment_metrics
