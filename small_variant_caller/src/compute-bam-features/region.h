#pragma once

#include <optional>
#include <string>

#include <xoos/types/int.h>

#include "util/region-util.h"

namespace xoos::svc {

struct Region {
  std::string chrom;  // chromosome
  u64 start{};        // 0-based inclusive start position on the reference
  u64 end{};          // 0-based exclusive end position on the reference

  // The previous interval is intended for checking whether an overlapping deletion should be processed in this region
  // or the previous region.
  std::optional<Interval> prev_interval = std::nullopt;  // previous interval on the same chromosome

  s64 Start() const;

  s64 EndInclusive() const;
};

/**
 * @brief Extracts variants from a VCF file that are within the specified region.
 * @param vcf_file Path to the VCF file
 * @param region Region to extract variants from
 * @param region_padding Padding to add to the region
 * @return Set of keys for variants within the specified region
 */
StrUnorderedSet ExtractVariantKeySet(const fs::path& vcf_file, const Region& region, u64 region_padding = 1000);

}  // namespace xoos::svc
