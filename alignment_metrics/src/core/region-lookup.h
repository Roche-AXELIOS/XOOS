#pragma once

#include <optional>

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "alignment-metrics-options.h"
#include "core/interval.h"
#include "core/region-labels-set.h"
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
                             bool disable_region_partitioning = false,
                             bool stratify_by_label = false);

  /**
   * Constructs a RegionLookupTable from a vector of super regions.
   * This is useful for manually creating a RegionLookupTable for testing purposes.
   * A non-empty vector is treated as if built from a BED file (HasBedFile() returns true); an empty
   * vector is treated as if no BED file was supplied (HasBedFile() returns false).
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

  /**
   * Determines whether the genomic interval `interval` on chromosome `chromosome` overlaps any
   * region of interest.
   *
   * Unlike IsSafeToCountAlignment, this query is independent of any super-region/sub-region index,
   * so it can test an interval located anywhere in the genome (e.g. a supplementary alignment
   * described by a read's SA tag, which may map to a different position or chromosome than the
   * primary alignment).
   *
   * When no BED file is supplied, the regions of interest span whole contigs, so any interval on a
   * known chromosome overlaps.
   *
   * @param chromosome The chromosome name the interval is located on (matched against SA tag rname).
   * @param interval The 0-based, half-open genomic interval to test.
   * @return true if the interval overlaps at least one region of interest; false otherwise.
   */
  bool OverlapsAnyRegionOfInterest(std::string_view chromosome, const Interval& interval) const;

  /**
   * Returns the union (bitwise-OR) of the label bit vectors of every region of interest on
   * `chromosome` that overlaps `interval`. Used to attribute read-level (non-pileup) metrics to
   * every label class a read spans (union attribution).
   *
   * Returns 0 when the table is unlabeled, when the chromosome is unknown, or when the interval
   * overlaps no labeled region.
   *
   * @param chromosome The chromosome the interval is located on.
   * @param interval The 0-based, half-open genomic interval to test.
   * @return The OR of the labels of all overlapping regions of interest.
   */
  u16 GetLabelsForInterval(std::string_view chromosome, const Interval& interval) const;

  /**
   * The region labels set describing all region labels seen while building this table. Empty when the
   * table is unlabeled (non-stratified mode).
   */
  const RegionLabelsSet& GetLabels() const;

  /**
   * Returns true when this table was built from a BED file and therefore restricts metrics to the
   * target regions it describes. Returns false when the table was built directly from the contigs in
   * the BAM header (no BED file), in which case every position on a known chromosome is a region of
   * interest and OverlapsAnyRegionOfInterest need not be consulted.
   */
  bool HasBedFile() const;

  // TODO: Move this to a test utility as only used by tests.
  size_t TotalRegionCount() const;

  // TODO: Move this to a test utility as only used by tests.
  size_t SuperRegionCount() const;

 private:
  /**
   * Builds `_regions_of_interest_by_chr` from the current super regions and the regions without
   * coverage, so all regions of interest can be queried by chromosome name.
   */
  void BuildRegionOfInterestIndex();

  /**
   * Sub-regions without coverage grouped by chromosome.
   */
  RegionsByChromosome _regions_without_coverage_by_chr;

  /**
   * All regions of interest (with and without coverage) grouped by chromosome name, sorted by start
   * position. Used by OverlapsAnyRegionOfInterest to test intervals located anywhere in the genome.
   */
  RegionsByChromosome _regions_of_interest_by_chr;

  /**
   * Super region is a group of regions that are processed together.
   * It is used to group regions that are close to each other in the genome.
   */
  vec<SuperRegion> _super_regions;
  /**
   * Total number of regions across all chromosomes.
   */
  size_t _total_region_count{};

  /**
   * True when the table was built from a BED file (restricted to target regions); false when built
   * directly from the BAM header contigs. See HasBedFile().
   */
  bool _has_bed_file{true};

  /**
   * Registry of region labels parsed from the BED file. Empty in non-stratified mode.
   */
  RegionLabelsSet _labels;
};

}  // namespace xoos::alignment_metrics
