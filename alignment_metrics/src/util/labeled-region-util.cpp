#include "util/labeled-region-util.h"

#include <string>
#include <string_view>
#include <utility>

#include "core/interval.h"
#include "io/bed-reader.h"

namespace xoos::alignment_metrics {

RegionsByChromosome MergeRegionsByLabel(RegionsByChromosome regions) {
  RegionsByChromosome merged_regions;
  MergeRegionsByLabel(std::move(regions),
                      [&merged_regions](const std::string_view chromosome, const Interval& interval) {
                        merged_regions[std::string{chromosome}].emplace_back(interval);
                      });
  return merged_regions;
}

}  // namespace xoos::alignment_metrics
