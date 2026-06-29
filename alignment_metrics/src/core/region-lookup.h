#pragma once

#include <optional>

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "alignment-metrics-options.h"
#include "core/interval.h"
#include "core/super-region.h"
#include "io/bed-reader.h"

namespace xoos::alignment_metrics {

/**
 * Each pileup entry is about 500 bytes in size. To limit memory usage to less than 2 GB per thread,
 * we use an estimated maximum region size of 2 million bases and warn the user if any region exceeds this size.
 */
constexpr size_t kRecommendedMaxRegionSize = 2'000'000;

/**
 * A lookup table that stores regions of interest grouped by chromosome,
 * along with super regions that group together nearby regions for more efficient processing.
 *
 * The lookup table also provides methods to determine if an alignment can be safely counted
 * towards the metrics of a given super region, avoiding double counting of alignments that
 * overlap with regions in previous super regions.
 */
class RegionLookupTable {
 public:
  /**
   * Returns the list of super regions that are stored in this lookup table.
   */
  const vec<SuperRegion>& GetSuperRegions() const;

  const RegionsByChromosome& GetRegionsWithoutCoverage() const;

  // returns the first sub-region that the trimmed alignment overlaps with
  bool IsAlignmentInAnyRegion(s64 rpos, s64 reference_length, size_t super_region_index, size_t sub_region_index) const;

  /**
   * Constructs a RegionLookupTable from a vector of super regions.
   */
  explicit RegionLookupTable(const fs::path& bam_path,
                             const std::optional<fs::path>& bed_path = std::nullopt,
                             const std::optional<fs::path>& reference_path = std::nullopt,
                             u32 region_size = kDefaultRegionSize,
                             bool hp_aware = true,
                             bool disable_region_partitioning = false);

  /**
   * Constructs a RegionLookupTable from a vector of super regions.
   * This is useful for manually creating a RegionLookupTable for testing purposes.
   * TODO: Move this to a test utility file.
   */
  explicit RegionLookupTable(const vec<SuperRegion>& super_regions);

  /**
   * Given an alignment spanning the interval `alignment_span` on chromosome `chromosome`,
   * this function determines if it is safe to count the alignment towards the metrics of the super region
   * identified by `super_region_index`.
   *
   * In particular, it checks if the alignment overlaps with a subregion in a previous super region. If it does, the
   * alignment should have already been counted by the metrics worker processing that previous super region, and thus we
   * should not count it again for the current super region to avoid double counting.

   * @param chromosome The chromosome on which the alignment is located.
   * @param alignment_span The genomic interval spanned by the alignment.
   * @param super_region_index The index of the super region for which we want to determine if it's safe to count the
   alignment.
   * @param sub_region_index The index of the subregion within the super region where the alignment is located.
   *
   * NOTE: This function assumes that the `alignment_span` overlaps with the subregion identified by
   * `super_region_index` and `sub_region_index`. It does not check for this condition.
   *
   * @return true if it is safe to count the alignment for the specified super region; false otherwise.
   */
  bool IsSafeToCountAlignment(std::string_view chromosome,
                              const Interval& alignment_span,
                              size_t super_region_index,
                              size_t sub_region_index) const;

  // TODO: Move this to a test utility as only used by tests.
  size_t TotalRegionCount() const;

  // TODO: Move this to a test utility as only used by tests.
  size_t SuperRegionCount() const;

 private:
  /**
   * Sub-regions without coverage grouped by chromosome.
   */
  RegionsByChromosome _regions_without_coverage_by_chr;

  /**
   * Super region is a group of regions that are processed together.
   * It is used to group regions that are close to each other in the genome.
   */
  vec<SuperRegion> _super_regions;
  /**
   * Total number of regions across all chromosomes.
   */
  size_t _total_region_count{};
};

}  // namespace xoos::alignment_metrics
