#pragma once

#include <optional>
#include <string>

#include <xoos/types/int.h>

#include "util/region-util.h"

namespace xoos::svc {

struct Region {
  std::string chrom;
  // 0-based inclusive start position on the reference
  u64 start{};
  // 0-based exclusive end position on the reference
  u64 end{};

  // The previous interval is intended for checking whether an overlapping deletion should be processed in this region
  // or the previous region.
  std::optional<Interval> prev_interval = std::nullopt;

  // Returns the 0-based inclusive start position, cast to s64.
  s64 Start() const;

  // Returns end + 1 to pad the query range by one base beyond the region boundary.
  // Used when querying BAM reads to include reads starting at the region boundary
  // for correct insertion depth computation at the end of the region.
  s64 EndPadded() const;
};

/**
 * @brief Extracts variants from a VCF file that are within the specified region.
 * @param vcf_file Path to the VCF file
 * @param region Region to extract variants from
 * @param region_padding Padding to add to the region
 * @return Set of keys for variants within the specified region
 */
// Default padding (in bases) applied to a region when extracting overlapping variants.
inline constexpr u64 kDefaultRegionPadding{1000};

StrUnorderedSet ExtractVariantKeySet(const fs::path& vcf_file,
                                     const Region& region,
                                     u64 region_padding = kDefaultRegionPadding);

}  // namespace xoos::svc
