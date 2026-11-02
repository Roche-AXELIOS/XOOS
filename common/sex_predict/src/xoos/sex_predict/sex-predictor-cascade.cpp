#include "sex-predictor-cascade.h"

#include <optional>
#include <variant>

#include <xoos/io/alignment-reader.h>
#include <xoos/log/logging.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>

#include "cascade-depth-profile.h"
#include "cascade-heterozygosity.h"
#include "cascade-residual-chry.h"
#include "cascade-resolve.h"
#include "sex-predictor.h"

namespace xoos::sex_predict {

SexPrediction PredictSexDetailed(const fs::path& bam_file_path,
                                 const u8 min_mapq,
                                 const u64 min_reads,
                                 const u32 num_threads,
                                 const SexPredictionMode mode,
                                 const SexPredictionPolicy policy) {
  using enum Sex;
  using enum SexEvidence;
  using enum SexUncertaintyReason;

  // `policy` has no effect in kNormal mode regardless of which step below
  // ends up resolving the call -- kOneCopyXNoY (the only reason policy acts
  // on) is only reachable via the heterozygosity step kNormal mode always
  // skips. Warn once, up front, rather than only on the specific path that
  // happens to reach the mode branch below.
  if (mode == SexPredictionMode::kNormal && policy != kDefaultSexPredictionPolicy) {
    Logging::Debug(
        "sex_predict: SexPredictionPolicy has no effect in kNormal mode -- kOneCopyXNoY is only reachable via the "
        "heterozygosity step this mode skips");
  }

  // Step 1: the cheap chrX/chrY ratio scan.
  const auto counts = ScanChrXYReadCounts(bam_file_path, min_mapq, min_reads, num_threads);
  if (!counts.has_value()) {
    return SexPrediction{.sex = kUnknown,
                         .evidence = kNoSexChromosomes,
                         .confidence = SexConfidence::kLow,
                         .reason = kSexChromosomesAbsent,
                         .measurements = std::nullopt};
  }
  const f64 xy_ratio = static_cast<f64>(counts->chr_x) / static_cast<f64>(counts->chr_y);
  const Sex ratio_sex = ClassifyByRatio(xy_ratio, kMinRatioForSexDet, kMinRatioForSexNA);

  // Step 2: chrY present -> Male, returned immediately. A female cannot gain a
  // Y, and LOY only ever pushes a male toward Female/Unknown, so this call is
  // already LOY-safe and never pays for the evidence below.
  if (ratio_sex == kMale) {
    return SexPrediction{.sex = kMale,
                         .evidence = kChrYPresent,
                         .confidence = SexConfidence::kHigh,
                         .reason = kNone,
                         .measurements = std::nullopt};
  }

  // Step 3: no chrY. Genuine female, LOY male, or a tumor that lost chrX
  // material -- indistinguishable from chrX/chrY alone. Profile windowed depth
  // across the sampled autosomes plus both allosomes in one taskflow pass.
  const io::AlignmentReader reader = io::OpenAlignmentReader(bam_file_path);
  const auto region_x = GetXChromosome(reader.hdr);
  const auto region_y = GetYChromosome(reader.hdr);
  if (!region_x.has_value() || !region_y.has_value()) {
    return SexPrediction{.sex = ratio_sex,
                         .evidence = kChrXYRatioFallback,
                         .confidence = SexConfidence::kLow,
                         .reason = kSexChromosomesAbsent,
                         .measurements = std::nullopt};
  }
  auto depth_profile_result =
      BuildDepthProfile(bam_file_path, reader, *region_x, *region_y, min_mapq, num_threads, ratio_sex);
  if (std::holds_alternative<SexPrediction>(depth_profile_result)) {
    return std::get<SexPrediction>(depth_profile_result);
  }
  const DepthProfile& depth_profile = std::get<DepthProfile>(depth_profile_result);
  SexMeasurements measurements = depth_profile.measurements;

  // Step 3b: residual chrY -- see cascade-residual-chry.h.
  if (const auto residual_call = TryResidualChrY(measurements); residual_call.has_value()) {
    return *residual_call;
  }

  // kNormal mode: no chrY by ratio or residual depth. chrX LOH and
  // culture-induced X-loss -- the two confounds heterozygosity exists to
  // catch -- do not occur in germline tissue, so there is nothing left for
  // heterozygosity to disambiguate (see SexPredictionMode's doc comment).
  // Resolve directly from the ratio, skipping the heterozygosity pileup
  // entirely -- `measurements`'s het fields stay zero-initialized, which is
  // the load-bearing proof (not just a matching answer) that heterozygosity
  // was never measured.
  if (mode == SexPredictionMode::kNormal) {
    // A Female call gets a real, margin-based confidence (this mode's whole
    // thesis is that the ratio alone is trustworthy here); an Unknown call is
    // already the ratio's own hedge, so it stays kLow rather than a
    // synthesized "high confidence because we didn't check" value.
    const SexConfidence confidence =
        ratio_sex == kFemale ? ConfidenceFromMargin(xy_ratio, kMinRatioForSexDet) : SexConfidence::kLow;
    const SexUncertaintyReason reason = ratio_sex == kUnknown ? kAmbiguousChrXYRatio : kNone;
    return SexPrediction{.sex = ratio_sex,
                         .evidence = kNormalModeRatio,
                         .confidence = confidence,
                         .reason = reason,
                         .measurements = measurements};
  }

  // Step 4: chrX heterozygosity, the decisive signal. Measured against the
  // autosome rate so that ancestry, coverage, contamination, and genome-wide
  // LOH burden cancel out of the comparison.
  auto het_result =
      MeasureHeterozygosity(bam_file_path, depth_profile, *region_x, measurements, min_mapq, num_threads, ratio_sex);
  if (std::holds_alternative<SexPrediction>(het_result)) {
    return std::get<SexPrediction>(het_result);
  }
  measurements = std::get<SexMeasurements>(het_result);

  // Step 5: resolve. chrY is already known absent at this point; see
  // ResolveFromHeterozygosity's own doc comment for the decision logic.
  return ResolveFromHeterozygosity(measurements, ratio_sex, policy);
}

SexPrediction PredictSexDetailed(const fs::path& bam_file_path,
                                 const u8 min_mapq,
                                 const u64 min_reads,
                                 const u32 num_threads,
                                 const SexPredictionPolicy policy) {
  return PredictSexDetailed(bam_file_path, min_mapq, min_reads, num_threads, kDefaultSexPredictionMode, policy);
}

SexPrediction PredictSexDetailed(const fs::path& bam_file_path,
                                 const u8 min_mapq,
                                 const u64 min_reads,
                                 const u32 num_threads,
                                 const SexPredictionMode mode) {
  return PredictSexDetailed(bam_file_path, min_mapq, min_reads, num_threads, mode, kDefaultSexPredictionPolicy);
}

SexPrediction PredictSexDetailed(const fs::path& bam_file_path,
                                 const u8 min_mapq,
                                 const u64 min_reads,
                                 const u32 num_threads) {
  return PredictSexDetailed(bam_file_path, min_mapq, min_reads, num_threads, kDefaultSexPredictionPolicy);
}

SexPrediction PredictSexDetailed(const fs::path& bam_file_path, const u8 min_mapq, const u64 min_reads) {
  return PredictSexDetailed(bam_file_path, min_mapq, min_reads, kDefaultNumThreads);
}

SexPrediction PredictSexDetailed(const fs::path& bam_file_path, const u8 min_mapq) {
  return PredictSexDetailed(bam_file_path, min_mapq, kDefaultMinReads);
}

SexPrediction PredictSexDetailed(const fs::path& bam_file_path) {
  return PredictSexDetailed(bam_file_path, kDefaultMinMapQ);
}

}  // namespace xoos::sex_predict
