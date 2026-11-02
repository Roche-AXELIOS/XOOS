#pragma once

#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

namespace xoos::cnc {
const s32 kVcfParsingOptionsDefaultNThreads = 1;
const s32 kVcfParsingOptionsDefaultSomaticNormalSampleMinDepth = 15;
const f64 kVcfParsingOptionsDefaultSomaticNormalSampleMinBAF = 0.367;
const f64 kVcfParsingOptionsDefaultSomaticNormalSampleMaxBAF = 0.633;
const s32 kVcfParsingOptionsDefaultSomaticTumorSampleMinDepth = 15;
const f64 kVcfParsingOptionsDefaultSomaticTumorSampleMinBAF = 0.367;
const f64 kVcfParsingOptionsDefaultSomaticTumorSampleMaxBAF = 0.633;
const s32 kVcfParsingOptionsDefaultGermlineNormalSampleMinDepth = 0;
const f64 kVcfParsingOptionsDefaultGermlineNormalSampleMinBAF = 0;
const f64 kVcfParsingOptionsDefaultGermlineNormalSampleMaxBAF = 1.0;
const s32 kVcfParsingOptionsDefaultTumorOnlyMinDepth = 10;
const f64 kVcfParsingOptionsDefaultTumorOnlyMinBAF = 0.1;
const f64 kVcfParsingOptionsDefaultTumorOnlyMaxBAF = 0.9;
const f64 kVcfParsingOptionsDefaultMinPopulationAF = 0.01;

/**
 * @brief Numeric thresholds for filtering variants read from a VCF.
 *
 */
struct BafFilterOptions {
  s32 normal_sample_min_depth = 0;
  f64 normal_sample_min_baf = 0;
  f64 normal_sample_max_baf = 0;
  s32 tumor_sample_min_depth = 0;
  f64 tumor_sample_min_baf = 0;
  f64 tumor_sample_max_baf = 0;

  bool force_enable_somatic_variant_parsing = false;
};

}  // namespace xoos::cnc
