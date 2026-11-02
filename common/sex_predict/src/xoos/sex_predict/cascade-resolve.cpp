#include "cascade-resolve.h"

#include <utility>

#include <xoos/log/logging.h>

namespace xoos::sex_predict {

SexConfidence ConfidenceFromMargin(const f64 value, const f64 gate) {
  return value >= gate * kHighConfidenceMarginMultiplier ? SexConfidence::kHigh : SexConfidence::kMedium;
}

SexPrediction ResolveFromHeterozygosity(const SexMeasurements& measurements,
                                        const Sex ratio_sex,
                                        const SexPredictionPolicy policy) {
  using enum SexEvidence;
  using enum SexUncertaintyReason;

  if (measurements.het_ratio >= kMinHetRatioForTwoCopies) {
    if (measurements.chr_x_depth_copies < kTwoCopyDepthThreshold) {
      // Heterozygosity wins: it cannot be faked by a depth artifact, whereas a
      // ~2x depth ratio is readily distorted by GC and mappability bias.
      Logging::Debug(
          "sex_predict: chrX heterozygosity indicates two or more copies but depth indicates {:.3f}; trusting "
          "heterozygosity",
          measurements.chr_x_depth_copies);
    }
    return SexPrediction{.sex = Sex::kFemale,
                         .evidence = kChrXHeterozygosity,
                         .confidence = ConfidenceFromMargin(measurements.het_ratio, kMinHetRatioForTwoCopies),
                         .reason = kNone,
                         .measurements = measurements};
  }

  SexUncertaintyReason reason = kAmbiguousZygosity;
  if (measurements.het_ratio <= kMaxHetRatioForOneCopy) {
    // A single copy of chrX with no chrY. LOY male, Turner/LoX female, tumor
    // chrX-LOH, and cultured X loss are mutually indistinguishable here.
    reason = measurements.chr_x_depth_copies >= kTwoCopyDepthThreshold ? kCopyNeutralLoh : kOneCopyXNoY;
  }

  // kOneCopyXNoY is the one irreducible case (see SexPredictionPolicy): under
  // kStrict, report Unknown instead of guessing. Every other reason keeps the
  // ratio's call regardless of policy -- they are measurement gaps, not
  // biological ambiguity, so there is no "wrong guess" to protect against.
  const Sex resolved_sex =
      (reason == kOneCopyXNoY && policy == SexPredictionPolicy::kStrict) ? Sex::kUnknown : ratio_sex;
  Logging::Debug(
      "sex_predict: chrX heterozygosity ({:.4f}) did not confirm two copies; resolved to {} with reason {} "
      "(policy={})",
      measurements.het_ratio,
      GetDescriptionForSex(resolved_sex),
      static_cast<u32>(std::to_underlying(reason)),
      static_cast<u32>(std::to_underlying(policy)));
  // Always kLow: het_ratio landed in or below the ambiguous band here, so
  // whatever `resolved_sex` says (the ratio's answer, or Unknown under
  // kStrict), it rests on unconfirmed evidence by construction.
  return SexPrediction{.sex = resolved_sex,
                       .evidence = kChrXHeterozygosity,
                       .confidence = SexConfidence::kLow,
                       .reason = reason,
                       .measurements = measurements};
}

}  // namespace xoos::sex_predict
