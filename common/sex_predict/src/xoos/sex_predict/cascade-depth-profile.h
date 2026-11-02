#pragma once

#include <cstddef>
#include <string>
#include <variant>

#include <xoos/io/alignment-reader.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "sex-predictor.h"

namespace xoos::sex_predict {

/// A recognized autosome contig sampled for the diploid depth baseline.
struct AutosomeContig {
  s32 tid;
  s64 len;
};

/// A contig included in the windowed depth scan (autosome or allosome).
struct ScanContig {
  s32 tid;
  s64 len;
  std::string name;
};

/// Sampled contigs plus the corroborating chrX/chrY depth-copy estimates and
/// the autosome dispersion measurement, built by BuildDepthProfile().
struct DepthProfile {
  vec<AutosomeContig> autosomes;
  vec<ScanContig> contigs;
  size_t x_index = 0;
  size_t y_index = 0;
  SexMeasurements measurements;
};

/// Cascade step 3: windowed depth profile across the sampled autosomes plus
/// both allosomes, and the autosome-baseline trust gate. Returns either the
/// profile, or a fully-formed SexPrediction for one of this step's own
/// fallback cases (no autosome baseline, zero median, or the autosomes
/// disagreeing beyond a dispersion threshold). See the README's "Depth
/// profile" section for the full derivation.
std::variant<DepthProfile, SexPrediction> BuildDepthProfile(const fs::path& bam_file_path,
                                                            const io::AlignmentReader& reader,
                                                            const RegionInfo& region_x,
                                                            const RegionInfo& region_y,
                                                            u8 min_mapq,
                                                            u32 num_threads,
                                                            Sex ratio_sex);

}  // namespace xoos::sex_predict
