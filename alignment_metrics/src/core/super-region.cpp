#include "core/super-region.h"

#include <algorithm>
#include <string>
#include <utility>

#include <csv.hpp>

#include <xoos/error/error.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/io/malloc-ptr.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "core/interval.h"

namespace xoos::alignment_metrics {

/**
 * @brief Creates a SuperRegion from a chromosome name and a vector of subregions.
 *
 * Constructs a SuperRegion by validating, sorting, and grouping subregions. The function ensures
 * subregions are non-overlapping and sorted by start position, then creates a super region spanning
 * from the earliest start to the latest end. SuperRegions enable processing multiple nearby genomic
 * regions in a single BAM iteration, reducing I/O overhead and improving performance while preventing
 * double counting through overlap validation.
 *
 * @param chromosome The chromosome/contig name where the subregions are located
 * @param sub_regions Vector of genomic intervals to group into a super region
 * @return A SuperRegion containing the sorted, validated subregions
 * @throws error::Error if sub_regions is empty or contains overlapping intervals
 */
SuperRegion CreateSuperRegion(const std::string& chromosome, const vec<Interval>& sub_regions) {
  if (sub_regions.empty()) {
    throw error::Error("Super region must have at least 1 sub-region");
  }
  // Validate that the sub-regions are sorted and non-overlapping
  auto sorted_sub_regions = sub_regions;
  std::ranges::sort(sorted_sub_regions, [](const Interval& a, const Interval& b) { return a.start < b.start; });
  const auto sub_regions_overlap =
      std::ranges::adjacent_find(sorted_sub_regions, [](const Interval& a, const Interval& b) {
        return a.end > b.start;
      }) != sorted_sub_regions.end();
  if (sub_regions_overlap) {
    throw error::Error("Super region must not have overlapping sub-regions");
  }
  return SuperRegion{.chromosome = chromosome,
                     .start = sub_regions.front().start,
                     .end = sub_regions.back().end,
                     .subregions = std::move(sorted_sub_regions)};
}

/**
 * @brief Creates an htslib region list from a SuperRegion for BAM file querying.
 *
 * Converts a SuperRegion into htslib's native hts_reglist_t format for efficient multi-region
 * BAM queries. Uses calloc/malloc (not new) because htslib's destructor uses free(), and mixing
 * allocation methods causes undefined behavior. The function creates an array of intervals and
 * transfers ownership to the reglist, enabling batch querying of all subregions in a single
 * sam_itr_regions() call.
 *
 * @param super_region The SuperRegion to convert to htslib format
 * @return A smart pointer to hts_reglist_t containing the region intervals
 */
io::HtsRegListPtr CreateHtsRegList(const SuperRegion& super_region) {
  const auto& sub_regions = super_region.subregions;
  // with htslib must use calloc / malloc to allocate the hts_pair because using free for destructor
  io::MallocPtr<hts_pair_pos_t[]> intervals{
      static_cast<hts_pair_pos_t*>(std::calloc(sub_regions.size(), sizeof(hts_pair_pos_t)))};
  for (size_t i = 0; i < sub_regions.size(); ++i) {
    intervals[i] = hts_pair_pos_t{sub_regions.at(i).start, sub_regions.at(i).end};
  }
  // with htslib must use calloc / malloc to allocate the reglist  because using free for destructor
  auto hts_reg_list_ptr =
      io::HtsRegListPtr{static_cast<hts_reglist_t*>(std::calloc(1, sizeof(hts_reglist_t))), io::HtsRegListDeleter{1}};
  hts_reg_list_ptr->reg = super_region.chromosome.c_str();
  // the intervals are owned by the reglist
  hts_reg_list_ptr->intervals = intervals.release();
  hts_reg_list_ptr->count = static_cast<u32>(sub_regions.size());
  return hts_reg_list_ptr;
}

}  // namespace xoos::alignment_metrics
