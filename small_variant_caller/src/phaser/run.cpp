#include "run.h"

#include <cstddef>
#include <utility>

#include <fmt/format.h>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "phase-engine.h"
#include "phaser/params.h"
#include "phaser/variant.h"
#include "vcf-io.h"

namespace xoos::svc::phaser {

namespace {

// A multi-sample VCF needs an explicit tumor read group so reads from the normal
// sample are not included in the phasing evidence.
void ValidateTumorReadGroup(const VcfData& vcf_data, const Params& params) {
  const s32 sample_count = vcf_data.header->GetNumSamples();
  if (params.workflow == Workflow::kTumorOnlyWgs && sample_count != 1) {
    throw error::Error(
        fmt::format("tumor-only-wgs phasing requires a single-sample VCF; found {} samples", sample_count));
  }
  if (params.workflow == Workflow::kTumorNormalWgs && sample_count > 1 && params.tumor_rg.empty()) {
    throw error::Error(
        "tumor-normal-wgs phasing requires --tumor-sample-name; specify --tumor-sample-name for a "
        "multi-sample VCF");
  }
}

// Drop phase blocks where no variant is PASS, keeping only blocks that rescue at least one
// non-PASS variant via a phased PASS variant. Logs a summary of the filtering result.
vec<PhaseBlock> ApplyRescueFiltering(vec<PhaseBlock>& blocks, const vec<Variant>& variants) {
  s32 rescued = 0;
  s32 dropped_blocks = 0;
  vec<PhaseBlock> kept_blocks;
  for (auto& block : blocks) {
    bool has_pass = false;
    for (const s32 vi : block.variant_indices) {
      if (variants[static_cast<size_t>(vi)].is_pass) {
        has_pass = true;
        break;
      }
    }
    if (!has_pass) {
      ++dropped_blocks;
      continue;
    }
    for (const s32 vi : block.variant_indices) {
      if (!variants[static_cast<size_t>(vi)].is_pass) {
        ++rescued;
      }
    }
    kept_blocks.push_back(std::move(block));
  }
  Logging::Info("Rescue filtering: kept {} blocks, dropped {} (all non-PASS), rescued {} filtered variants",
                kept_blocks.size(),
                dropped_blocks,
                rescued);
  return kept_blocks;
}

}  // namespace

void Run(const Params& params) {
  Logging::Info("Reading VCF: {}", params.vcf_path.string());
  // ReadVcf infers the tumor sample column from the VCF header metadata (`tumor_sample`).
  const auto vcf_data = ReadVcf(params.vcf_path, params);
  ValidateTumorReadGroup(vcf_data, params);

  Logging::Info("Building phase groups (max distance: {}bp)", params.max_phase_distance);
  const auto groups = BuildPhaseGroups(vcf_data.variants, params.max_phase_distance);

  Logging::Info("Running phasing on {} groups with {} threads", groups.size(), params.threads);
  vec<PhaseBlock> blocks;
  if (params.threads > 1) {
    blocks = RunPhasingParallel(params.bam_path, params.tumor_rg, vcf_data.variants, groups, params);
  } else {
    blocks = RunPhasing(params.bam_path, params.tumor_rg, vcf_data.variants, groups, params);
  }

  // --rescue-filtered: drop phase blocks where no variant is PASS.
  if (params.rescue_filtered) {
    blocks = ApplyRescueFiltering(blocks, vcf_data.variants);
  }

  Logging::Info("Writing output VCF: {}", params.output_path.string());
  WriteVcf(params.output_path, vcf_data, blocks, vcf_data.variants, params);

  s32 phased_variants = 0;
  for (const auto& block : blocks) {
    phased_variants += static_cast<s32>(block.variant_indices.size());
  }

  Logging::Info("Phasing complete: {} phase blocks, {} phased variants out of {} total",
                blocks.size(),
                phased_variants,
                vcf_data.variants.size());
}

}  // namespace xoos::svc::phaser
