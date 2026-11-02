#pragma once

#include <xoos/types/float.h>
#include <xoos/types/fs.h>

namespace xoos::cnc {

const bool kAugmentIntervalsDefaultWholeGenome = false;
const size_t kAugmentIntervalsDefaultWholeGenomeIntervalSize = 1000;
const f64 kAugmentIntervalsDefaultMinGcContent = 0.15;
const f64 kAugmentIntervalsDefaultMaxGcContent = 0.75;
const f64 kAugmentIntervalsDefaultWholeGenomeMinMappability = 0.6;
const size_t kAugmentIntervalsDefaultOnTargetIntervalSize = 400;
const size_t kAugmentIntervalsDefaultOffTargetMinWidth = 20000;
const size_t kAugmentIntervalsDefaultOffTargetTrimLength = 500;
const size_t kAugmentIntervalsDefaultOffTargetIntervalSize = 200000;
const f64 kAugmentIntervalsDefaultOnTargetMinMappability = 0.6;
const f64 kAugmentIntervalsDefaultOffTargetMinMappability = 0.1;
const f64 kAugmentIntervalsDefaultSexChromosomeTargetMinMappability = 0.9;
const bool kAugmentIntervalsDefaultNoOffTargets = false;
const bool kAugmentIntervalsDefaultIncludeAltContigs = false;
const size_t kAugmentIntervalsGermlineWGSDefaultIntervalSize = 500;
/// Whole-genome interval size used for germline WGS on SBX-S (simplex) input.
/// Simplex reads have different coverage/quality characteristics, so a larger
/// interval size than the SBX-D default (::kAugmentIntervalsGermlineWGSDefaultIntervalSize)
/// is used when the input BAM is detected as SBX-S.
const size_t kAugmentIntervalsGermlineWGSSbxSIntervalSize = 1000;
const f64 kAugmentIntervalsGermlineWGSDefaultWholeGenomeMinMappability = 0.8;

struct AugmentIntervalsOptions {
  size_t whole_genome_interval_size{};
  f64 whole_genome_min_mappability{};
  size_t on_target_interval_size = kAugmentIntervalsDefaultOnTargetIntervalSize;
  size_t off_target_min_width = kAugmentIntervalsDefaultOffTargetMinWidth;
  size_t off_target_trim_length = kAugmentIntervalsDefaultOffTargetTrimLength;
  size_t off_target_interval_size = kAugmentIntervalsDefaultOffTargetIntervalSize;
  f64 min_gc_content = kAugmentIntervalsDefaultMinGcContent;
  f64 max_gc_content = kAugmentIntervalsDefaultMaxGcContent;
  f64 on_target_min_mappability = kAugmentIntervalsDefaultOnTargetMinMappability;
  f64 off_target_min_mappability = kAugmentIntervalsDefaultOffTargetMinMappability;
  f64 sex_chromosome_min_mappability = kAugmentIntervalsDefaultSexChromosomeTargetMinMappability;
  bool generate_whole_genome_intervals = kAugmentIntervalsDefaultWholeGenome;
  bool no_off_targets = kAugmentIntervalsDefaultNoOffTargets;
  bool include_alt_contigs = kAugmentIntervalsDefaultIncludeAltContigs;
};
}  // namespace xoos::cnc
