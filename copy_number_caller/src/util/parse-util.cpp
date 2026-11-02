#include "util/parse-util.h"

#include <cstddef>
#include <string>
#include <tuple>

#include <xoos/genome/region.h>

namespace xoos::cnc {

bool IsMissingValue(const std::string_view value) {
  return util::IsNA(value) || value == kLegacyMissingValue;
}

/**
 * @brief parses a 1-based closed "contig:start-end" region string into a
 * (contig, start, end) tuple with 0-based half-open coordinates (start
 * inclusive, end exclusive, as in a BED file). Only the explicit range form is
 * accepted; whole-chromosome and single-position inputs are rejected.
 * @param r  region string in the form "contig:start-end" (1-based closed)
 * @return tuple - contig, 0-based start, 0-based exclusive end
 */
std::tuple<std::string, size_t, size_t> ParseRegionString(const std::string& r) {
  // This caller has always required an explicit "chr:start-end" range; restrict
  // the parser to that form so malformed input fails fast.
  const auto region = genome::ParseSamtoolsLikeRegion(r, genome::RegionForm::kRange);
  return std::make_tuple(region.chromosome, static_cast<size_t>(region.start), static_cast<size_t>(region.end));
}

}  // namespace xoos::cnc
