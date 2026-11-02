#pragma once

#include <filesystem>
#include <optional>

#include <xoos/types/float.h>
#include <xoos/types/int.h>

#include "sex.h"

namespace xoos::sex_predict {

namespace fs = std::filesystem;

// --- chrX heterozygosity ---------------------------------------------------
//
// The primary discriminator between one and two copies of chrX -- immune to
// the GC/mappability bias that rules out a depth-only estimate. See the
// README's "Why heterozygosity, not depth" and "chrX heterozygosity" sections.

// Number of sampled windows, and their width, per chromosome group (chrX, and
// the autosomes collectively).
constexpr u32 kNumHetWindowsPerGroup = 200;
constexpr s64 kHetWindowSize = 100'000;

// Minimum evaluated sites in *each* group before the heterozygosity ratio is
// trusted; below this the call falls back to the chrX/chrY ratio.
constexpr u64 kMinEvaluatedHetSites = 2'000;

// Gates on het_rate(chrX) / het_rate(autosomes) -- two copies run ~0.7, one
// copy runs 0.01-0.03. See the README's "Threshold derivation" section for
// the population-genetics basis and calibration status.
constexpr f64 kMinHetRatioForTwoCopies = 0.35;
constexpr f64 kMaxHetRatioForOneCopy = 0.10;

// Estimated chrX copies at or above which the (corroborating) depth signal is
// read as "two or more copies".
constexpr f64 kTwoCopyDepthThreshold = 1.5;

// --- Residual chrY (partial loss-of-Y) --------------------------------------
//
// A tumor with LOY is rarely 100% pure -- real samples retain some fraction of
// Y-positive stromal/immune cells, leaving a measurable, if reduced, chrY
// depth. See the README's "Residual chrY" section for the full derivation.

// Estimated chrY copies at or above which residual coverage is trusted as real
// signal rather than cross-mapping noise. Calibrated against an in-silico
// COLO829 dilution series (see README); UNVALIDATED FOR MOSAIC LOSS-OF-Y
// (mLOY) -- also used in `SexPredictionMode::kNormal`, where the analogous
// confound is mLOY clonal fraction rather than tumor purity. See the README's
// "Known limitations" for the caveat and the concrete next validation step.
constexpr f64 kMinResidualChrYCopiesForMale = 0.05;

// How far a measurement must clear its gate to be read as `SexConfidence::kHigh`
// rather than `kMedium`. A heuristic margin, not itself calibrated -- see the
// README's "SexConfidence" section for the mapping and real-data basis.
constexpr f64 kHighConfidenceMarginMultiplier = 1.5;

/// Which parts of the cascade PredictSexDetailed() runs. Steps 1-3 (ratio,
/// depth profile, residual-chrY) measure allosome dosage and always run;
/// only chrX heterozygosity (step 4) is mode-gated. Named around which BAM is
/// being scanned, not the caller's pipeline label -- see the README's
/// "SexPredictionMode" section for the full rationale and the naming story.
enum class SexPredictionMode : u8 {
  /// Skip chrX heterozygosity: the confounds it exists to catch (tumor chrX
  /// LOH, culture-induced X-loss) cannot occur in germline tissue, so there
  /// is nothing left for it to disambiguate. `SexPredictionPolicy` still
  /// parses but has no effect (see its own doc comment). Residual chrY depth
  /// is still checked -- see `kMinResidualChrYCopiesForMale`'s mLOY caveat.
  kNormal,
  /// Today's full cascade, unchanged.
  kCascade
};

/// Default mode: preserves the cascade's original (and only, until this
/// parameter existed) behavior.
constexpr SexPredictionMode kDefaultSexPredictionMode = SexPredictionMode::kCascade;

/// Which evidence PredictSexDetailed() reached its call on. See the README's
/// "Cascade" and "Resolving to a call" sections for the full derivation of
/// each.
enum class SexEvidence : u8 {
  /// chrY present -- Male. LOY only ever pushes toward Female/Unknown, so
  /// this call is already LOY-safe.
  kChrYPresent,
  /// Ratio didn't clear the Male threshold, but depth found chrY coverage
  /// above `kMinResidualChrYCopiesForMale` -- typically partial (not
  /// complete) LOY.
  kResidualChrYDepth,
  /// `SexPredictionMode::kNormal` only: neither the initial ratio nor residual
  /// chrY resolved Male, and heterozygosity was deliberately not measured.
  /// The result therefore comes directly from the chrX/chrY ratio. Distinct
  /// from `kChrXYRatioFallback`: this is a mode choice, not a measurement
  /// failure.
  kNormalModeRatio,
  /// No chrY; the call comes from chrX heterozygosity vs. the autosome rate.
  kChrXHeterozygosity,
  /// Heterozygosity could not be measured well enough to use, so the call
  /// falls back to the plain chrX/chrY ratio.
  kChrXYRatioFallback,
  /// chrX and/or chrY absent from the BAM header; no call possible.
  kNoSexChromosomes
};

/// Controls how PredictSexDetailed() resolves the one case that is genuinely
/// irreducible from a BAM alone (`SexUncertaintyReason::kOneCopyXNoY`). See
/// the README's "SexPredictionPolicy" section for the full rationale.
enum class SexPredictionPolicy : u8 {
  /// Keep the chrX/chrY ratio's call (default). Never regresses a call the
  /// ratio already gets right (e.g. HCC1395). `reason` still reports
  /// `kOneCopyXNoY` so a caller that wants to act on the ambiguity can.
  kReportBestCall,
  /// Report `Unknown` for `kOneCopyXNoY` instead of guessing -- for callers
  /// where a confidently wrong Female is worse than an Unknown.
  kStrict
};

/// Default policy: keep the chrX/chrY ratio's call rather than demoting an
/// unconfirmed one-copy-X result to Unknown.
constexpr SexPredictionPolicy kDefaultSexPredictionPolicy = SexPredictionPolicy::kReportBestCall;

/// Any uncertainty attached to a call, including retained Female ratio calls.
/// Machine-readable so callers can distinguish "genuinely ambiguous biology"
/// from "we could not measure it" -- see the README's "Resolving to a call"
/// section.
enum class SexUncertaintyReason : u8 {
  /// No uncertainty was recorded.
  kNone,
  /// One copy of chrX and no chrY -- the honest, irreducible core of the
  /// design (see `SexPredictionPolicy`).
  kOneCopyXNoY,
  /// Depth says two or more copies of chrX but heterozygosity says one:
  /// copy-neutral loss of heterozygosity, common in tumors.
  kCopyNeutralLoh,
  /// The heterozygosity ratio landed between the one-copy and two-copy gates.
  kAmbiguousZygosity,
  /// `SexPredictionMode::kNormal` stopped after the chrX/chrY ratio landed in
  /// its explicit Unknown band.
  kAmbiguousChrXYRatio,
  /// Too few sites passed the depth and quality filters to trust a rate.
  kInsufficientHetSites,
  /// Fewer than two autosomes carried usable coverage, so their consistency
  /// and the depth baseline could not be established.
  kNoAutosomeBaseline,
  /// The sampled autosomes' depths disagreed beyond `kMaxAutosomeDepthCv`
  /// (e.g. an aneuploid tumor), so the baseline was not trustworthy.
  kAutosomeDepthDispersed,
  /// chrX and/or chrY absent from the BAM header.
  kSexChromosomesAbsent
};

/// The measured quantities behind a call, for logging and threshold
/// calibration. Populated whenever the depth profile ran.
struct SexMeasurements {
  /// het_rate(chrX) / het_rate(autosomes) -- the primary discriminator.
  f64 het_ratio = 0.0;
  f64 chr_x_het_rate = 0.0;
  f64 autosome_het_rate = 0.0;
  u64 chr_x_evaluated_sites = 0;
  u64 autosome_evaluated_sites = 0;
  /// Continuous (unrounded) copy estimates from median window depth.
  f64 chr_x_depth_copies = 0.0;
  f64 chr_y_depth_copies = 0.0;
  /// Dispersion of the sampled autosomes' depths; the baseline trust gate.
  f64 autosome_depth_cv = 0.0;
};

/// How much to trust a call, independent of what the call is or which
/// `SexPredictionPolicy` produced it. Categorical rather than a numeric score
/// deliberately -- a precise float would imply a calibration precision this
/// module does not have. Derived entirely from measurements already
/// computed; no extra scanning. See the README's "SexConfidence" section for
/// the full evidence -> confidence mapping and its real-data basis.
enum class SexConfidence : u8 {
  kHigh,
  kMedium,
  kLow
};

/// Result of the predictor: the call, the evidence it rests on, any uncertainty
/// attached to it, how much to trust it, and the underlying measurements.
struct SexPrediction {
  Sex sex = Sex::kUnknown;
  SexEvidence evidence = SexEvidence::kNoSexChromosomes;
  SexConfidence confidence = SexConfidence::kLow;
  SexUncertaintyReason reason = SexUncertaintyReason::kSexChromosomesAbsent;
  std::optional<SexMeasurements> measurements;
};

/// Robust sex prediction from a BAM alone, safe for tumor samples. Runs a
/// 4-step cascade (chrY present -> residual chrY depth -> chrX heterozygosity
/// -> `SexPredictionPolicy`), strongest signal first; `mode` gates only the
/// heterozygosity step. See the README's "Cascade" and "Resolving to a call"
/// sections for the full derivation and the table of which state resolves to
/// what.
///
/// WGS only; the depth baseline assumes roughly uniform genome-wide coverage.
/// An explicit sex supplied by the caller should always take precedence over
/// this prediction.
SexPrediction PredictSexDetailed(const fs::path& bam_file_path,
                                 u8 min_mapq,
                                 u64 min_reads,
                                 u32 num_threads,
                                 SexPredictionMode mode,
                                 SexPredictionPolicy policy);

/// Overload using `kDefaultSexPredictionPolicy`, with `mode` still explicit.
SexPrediction PredictSexDetailed(  // NOSONAR cpp:S5536 -- used in sex-predictor-tests.cpp (5-arg calls with a mode)
    const fs::path& bam_file_path,
    u8 min_mapq,
    u64 min_reads,
    u32 num_threads,
    SexPredictionMode mode);

/// Overload using `kDefaultSexPredictionMode` (`kCascade`).
SexPrediction PredictSexDetailed(
    const fs::path& bam_file_path, u8 min_mapq, u64 min_reads, u32 num_threads, SexPredictionPolicy policy);

/// Overload using `kDefaultSexPredictionPolicy`.
SexPrediction PredictSexDetailed(const fs::path& bam_file_path, u8 min_mapq, u64 min_reads, u32 num_threads);

/// Overload using `kDefaultNumThreads`.
SexPrediction PredictSexDetailed(const fs::path& bam_file_path, u8 min_mapq, u64 min_reads);

/// Overload using `kDefaultMinReads` (0, no early stopping).
SexPrediction PredictSexDetailed(const fs::path& bam_file_path, u8 min_mapq);

/// Overload using `kDefaultMinMapQ`.
SexPrediction PredictSexDetailed(  // NOSONAR cpp:S5536 -- widely used in sex-predictor-tests.cpp
    const fs::path& bam_file_path);

}  // namespace xoos::sex_predict
