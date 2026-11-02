#pragma once

#include <variant>

#include "allele-counter.h"
#include "cascade-depth-profile.h"
#include "sex-predictor-cascade.h"

namespace xoos::sex_predict {

/// `count` evenly spaced windows of `size` bases within [start, end). When the
/// span cannot accommodate that many disjoint windows, the whole span is
/// returned as a single window instead.
vec<ScanWindow> SampleWindows(s32 tid, s64 start, s64 end, u32 count, s64 size);

/// Cascade step 4: chrX heterozygosity vs. the autosome rate, measured
/// against `depth_profile`'s sampled contigs. Returns the updated
/// measurements (with het fields populated), or a fully-formed SexPrediction
/// if too few sites were evaluable to trust a rate. See the README's "chrX
/// heterozygosity" section.
std::variant<SexMeasurements, SexPrediction> MeasureHeterozygosity(const fs::path& bam_file_path,
                                                                   const DepthProfile& depth_profile,
                                                                   const RegionInfo& region_x,
                                                                   SexMeasurements measurements,
                                                                   u8 min_mapq,
                                                                   u32 num_threads,
                                                                   Sex ratio_sex);

}  // namespace xoos::sex_predict
