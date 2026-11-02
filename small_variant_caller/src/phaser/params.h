#pragma once

#include <memory>
#include <optional>
#include <string>

#include <xoos/io/metadata-util.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "core/config.h"
#include "core/workflow.h"

namespace xoos::svc::phaser {

// Default maximum reference distance (bp) between variants considered for phasing.
inline constexpr s32 kDefaultMaxPhaseDistance = 10;
// Default minimum fraction of shared reads that must be in cis to phase a pair.
inline constexpr f32 kDefaultMinCisFraction = 0.8F;
// Match the SVC tumor-normal-wgs minimum base quality.
inline constexpr s32 kDefaultMinBaseQuality = static_cast<s32>(kConcordantMinBaseQuality);

struct Params {
  fs::path bam_path;
  fs::path vcf_path;
  fs::path output_dir{"."};
  fs::path vcf_output{"output.phased.vcf.gz"};
  fs::path output_path;
  // Optional reference FASTA for merging with gaps.
  fs::path ref_path;
  std::string tumor_rg;
  s32 max_phase_distance = kDefaultMaxPhaseDistance;
  s32 min_shared_reads = 2;
  s32 min_cis_reads = 2;
  f32 min_cis_fraction = kDefaultMinCisFraction;
  bool report_mnv = false;
  bool merge_mnv = true;
  bool merge_complex = true;
  bool keep_components = true;
  bool include_filtered = false;
  // Phase PASS + non-PASS pairs; rescue non-PASS if phased.
  bool rescue_filtered = false;
  // Skip phasing insertion + SNV pairs (reduces FPs near insertions).
  bool skip_ins_snv = false;
  // Minimum base quality at the variant position to count a read.
  s32 min_base_quality = kDefaultMinBaseQuality;
  size_t threads = 1;
  std::optional<io::CommandLineInfo> command_line;
  Workflow workflow = Workflow::kTumorNormalWgs;
};

using ParamsPtr = std::shared_ptr<Params>;

}  // namespace xoos::svc::phaser
