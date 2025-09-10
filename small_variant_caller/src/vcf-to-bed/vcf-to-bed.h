#pragma once

#include <optional>

#include <xoos/types/fs.h>

#include "util/region-util.h"

namespace xoos::svc {

// CLI parameters for VCF to BED conversion
struct VcfToBedParam {
  fs::path vcf_file{};         // path to input VCF file
  fs::path output_bed_file{};  // path to output BED file
  u32 left_pad{};              // left-padding to variant start position
  u32 right_pad{};             // right-padding to variant end position
  u32 collapse_dist{};         // distance to collapse nearby variants after padding is applied
  std::optional<ChromIntervalsMap> chrom_intervals_map{};  // map of chromosome to vector of target region intervals
  size_t threads{};                                        // number of threads to use for parallel processing
};

ChromIntervalsMap ExtractVariantIntervalsSingleThreaded(const fs::path& vcf_path,
                                                        const ChromIntervalsMap& target_intervals,
                                                        u32 left_pad = 0,
                                                        u32 right_pad = 0,
                                                        u32 collapse_dist = 0);
ChromIntervalsMap ExtractVariantIntervalsParallelized(const fs::path& vcf_path,
                                                      const ChromIntervalsMap& target_regions,
                                                      u32 threads,
                                                      u32 left_pad = 0,
                                                      u32 right_pad = 0,
                                                      u32 collapse_dist = 0);
void ConvertVcfToBed(const VcfToBedParam& param);
}  // namespace xoos::svc
