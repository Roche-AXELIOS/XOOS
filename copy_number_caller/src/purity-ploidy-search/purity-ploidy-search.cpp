#include "purity-ploidy-search/purity-ploidy-search.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <utility>

#include <xoos/enum/enum-util.h>
#include <xoos/error/error.h>
#include <xoos/io/metadata-util.h>
#include <xoos/log/logging.h>
#include <xoos/sex_predict/sex.h>
#include <xoos/stats/copy-number-stats.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/util/file-functions.h>
#include <xoos/util/missing-value.h>

#include "cli/copy-number-caller-options.h"
#include "io/copy-number-caller-default-filenames.h"
#include "io/write-segment-likelihoods.h"
#include "likelihood/baf-tools.h"
#include "likelihood/likelihood-mode.h"
#include "likelihood/likelihood.h"
#include "likelihood/total-copy-number-prior.h"
#include "observations/observations.h"
#include "segmentation/genomic-segments.h"
#include "segmentation/interval-trees.h"
#include "segmentation/read-segments.h"
#include "segmentation/segment-type.h"
#include "util/eigen-math.h"
#include "util/parse-util.h"

namespace xoos::cnc {
// Grid axis constants (kLikelihoodRo* / kLikelihoodPhi*) live in purity-ploidy-search.h so the search and the
// per-segment likelihood diagnostic share one definition of the grid.
const f64 kMinIntervalDensity = 0.5;
constexpr s32 kPercentMultiplier = 100;
const size_t kMinCandidateHetDelSegments = 5;
const f64 kHetDelExpectedTcnMin = 0.9;
const f64 kHetDelExpectedTcnMax = 1.1;
const f64 kHetDelMaxMBAF = 0.45;
const f64 kWGDMaxPeakDifference = 0.1;
const f64 kMinPurityToEvaluateWgdQc = 0.35;
// Lower purity bound for the low-purity LOH fallback. Below this purity even the LOH signal is too
// unreliable to distinguish WGD from noise, so the diploid-bias solution is accepted without running the LOH check.
const f64 kMinPurityForLohFallback = 0.2;
const size_t kNumTopPeaks = 2;
constexpr s32 kPurityPloidyOutputPrecision = 2;

/**
 * @brief Use global tumor purity/ploidy estimates (user-provided or optionally estimate those here) to calculate the
 * joint logr-BAF likelihood of each segment. In calculating this joint likelihood we also return the most likely total
 * copy number of the segment as well as the the most likely minor copy number. See docs/likelihood.md for more
 * information
 * @param params
 */
void PurityPloidySearchMain(const CopyNumberCallerOptions& options) {
  if (!options.logrs_fname.has_value()) {
    throw error::Error("PurityPloidySearch requires a logRs file (--logrs)");
  }
  if (!options.segments_fname.has_value()) {
    throw error::Error("PurityPloidySearch requires a segments file (--segments)");
  }
  if (!options.sample_metadata_options.sex.has_value()) {
    throw error::Error("PurityPloidySearch requires the sample sex (--sex)");
  }
  const sex_predict::Sex sex = options.sample_metadata_options.sex.value();
  std::ifstream logr_ifstream(options.logrs_fname.value());
  const Observations logrs = ReadObservations(logr_ifstream);

  // create output files
  const auto purity_ploidy_out_fname =
      file::PrefixedOutputPath(options.output_dir, options.output_prefix, kDefaultPurityPloidyOutput);
  const auto purity_ploidy_grid_out_fname =
      file::PrefixedOutputPath(options.output_dir, options.output_prefix, kDefaultPurityPloidyGridOutput);
  vec<fs::path> output_paths{purity_ploidy_out_fname, purity_ploidy_grid_out_fname};
  if (options.purity_ploidy_search_options.write_segment_likelihood) {
    output_paths.push_back(file::PrefixedOutputPath(
        options.output_dir, options.output_prefix, kDefaultSegmentLikelihoodDiploidBiasOutput));
    output_paths.push_back(
        file::PrefixedOutputPath(options.output_dir, options.output_prefix, kDefaultSegmentLikelihoodUniformOutput));
  }
  file::CheckOutputFilePathsForWriting(output_paths, options.overwrite);

  std::optional<Observations> ref_ads;
  std::optional<Observations> alt_ads;
  SegmentType segment_type = SegmentType::kUnknown;
  if (!options.bafs_fname.has_value()) {
    segment_type = SegmentType::kLogROnly;
  } else {
    auto ref_alts = GetDepthsFromFile(options.bafs_fname.value());
    ref_ads = std::move(ref_alts.ref_obvs);
    alt_ads = std::move(ref_alts.alt_obvs);
    segment_type = SegmentType::kBaf;
  }
  auto segments = ReadSegments(options.segments_fname.value(), segment_type);
  // Recompute mean_logr and num_obs from raw per-locus observations so that all downstream
  // steps (ShouldSkipSegment, EvaluatePurityPloidySolution, grid search) use consistent
  // full-precision values instead of the rounded MeanLogRatio from the SEG file.
  PopulateGenomicSegmentOptionalFields(segments, logrs);
  const PurityPloidySearchOut purity_ploidy_out = PurityPloidySearchWithWGDQC(segments,
                                                                              logrs,
                                                                              ref_ads,
                                                                              alt_ads,
                                                                              options.purity_ploidy_search_options,
                                                                              purity_ploidy_grid_out_fname,
                                                                              options.output_prefix,
                                                                              sex,
                                                                              options.command_line_info);
  WritePurityPloidy(purity_ploidy_out.purity,
                    purity_ploidy_out.ploidy,
                    purity_ploidy_out_fname,
                    options.command_line_info,
                    purity_ploidy_out.total_copy_number_prior);
  if (ref_ads.has_value() && alt_ads.has_value()) {
    CheckExtremeBafProportionAndWarn(ref_ads.value(), alt_ads.value(), purity_ploidy_out.purity);
  }
}

PurityPloidySearchOut PurityPloidySearch(const vec<GenomicSegment>& segments,
                                         const Observations& logrs,
                                         const std::optional<Observations>& ref_obvs,
                                         const std::optional<Observations>& alt_obvs,
                                         const PurityPloidyGridParams& params,
                                         const fs::path& grid_out_fname,
                                         const io::CommandLineInfo& command_line_info) {
  f64 purity = NAN;
  f64 ploidy = NAN;
  Eigen::MatrixXd ppgrid_joint;
  Eigen::MatrixXd ppgrid_logr;
  Eigen::MatrixXd ppgrid_baf;
  const TotalCopyNumberPrior prior =
      params.use_diploid_bias_prior ? TotalCopyNumberPrior::kDiploidBias : TotalCopyNumberPrior::kUniform;
  if (ref_obvs.has_value() && alt_obvs.has_value()) {
    std::tie(ppgrid_joint, ppgrid_logr, ppgrid_baf) =
        GetPurityPloidyGridSerialLogRBAF(segments, logrs, ref_obvs.value(), alt_obvs.value(), params);
    std::tie(purity, ploidy) = EstimatePurityPloidyByTotalMaximum(ppgrid_joint);
    WritePurityPloidyGrid(ppgrid_joint, ppgrid_logr, ppgrid_baf, grid_out_fname, command_line_info, prior);
  } else {
    ppgrid_joint = GetPurityPloidyGridLogROnly(segments, logrs, params);
    std::tie(purity, ploidy) = EstimatePurityPloidyByTotalMaximum(ppgrid_joint);
    WritePurityPloidyGrid(ppgrid_joint, {}, {}, grid_out_fname, command_line_info, prior);
  }
  return {.purity = purity,
          .ploidy = ploidy,
          .joint_ll_grid = ppgrid_joint,
          .logr_ll_grid = ppgrid_logr,
          .baf_ll_grid = ppgrid_baf,
          .total_copy_number_prior = prior};
}

/**
 * @brief Checks if a segment should be skipped based on mean_logr, num_obs, and optionally num_snps.
 * @param seg The segment to check.
 * @param seg_max_logr Maximum allowed mean_logr.
 * @param seg_min_num_logrs Minimum required number of logRs.
 * @param seg_min_num_snps Minimum required number of SNPs.
 * @param check_seg_min_snps Whether to check for minimum required number of SNPs.
 * @return bool True if the segment should be skipped, false otherwise.
 */
bool ShouldSkipSegment(const GenomicSegment& seg,
                       f64 seg_max_logr,
                       size_t seg_min_num_logrs,
                       std::optional<size_t> seg_min_num_snps,
                       const bool check_seg_min_snps) {
  if (seg.mean_logr > seg_max_logr) {
    Logging::Info("Skipping segment {}:{}-{} because mean logr value is greater than {}",
                  seg.contig,
                  seg.start,
                  seg.end,
                  seg_max_logr);
    return true;
  }

  if (!seg.num_obs.has_value()) {
    throw error::Error("Segment {}:{}-{} is missing its observation count (num_obs)", seg.contig, seg.start, seg.end);
  }
  if (seg.num_obs.value() < seg_min_num_logrs) {
    Logging::Info(
        "Skipping segment {}:{}-{} because it has < {} logr", seg.contig, seg.start, seg.end, seg_min_num_logrs);
    return true;
  }
  if (check_seg_min_snps) {
    if (!seg_min_num_snps.has_value()) {
      throw error::Error("ShouldSkipSegment requires seg_min_num_snps when check_seg_min_snps is set");
    }
    if (!seg.num_snps.has_value()) {
      throw error::Error("Segment {}:{}-{} is missing its SNP count (num_snps)", seg.contig, seg.start, seg.end);
    }
    if (seg.num_snps.value() < seg_min_num_snps.value()) {
      Logging::Info("Skipping segment {}:{}-{} because it has < {} SNPs",
                    seg.contig,
                    seg.start,
                    seg.end,
                    seg_min_num_snps.value());
      return true;
    }
  }
  return false;
}

/**
 * @brief Checks whether a segment has low interval density relative to its length.
 * A segment passes (returns false) when the total base coverage of the logR intervals that overlap the segment is at
 * least kMinIntervalDensity * segment_length. Segments that fail this check are large regions with very sparse interval
 * coverage and should be excluded from purity/ploidy estimation.
 * @param logr_idxs Indices into logrs for the intervals that overlap the segment.
 * @param logrs All logR observations (provides genomic positions via starts/ends).
 * @param seg The segment being evaluated.
 * @return true if the interval density is too low (segment should be skipped), false otherwise.
 */
bool HasLowIntervalDensity(const vec<size_t>& logr_idxs, const Observations& logrs, const GenomicSegment& seg) {
  if (seg.end == seg.start) {
    return false;
  }
  size_t total_covered_bases = 0;
  for (const size_t idx : logr_idxs) {
    const auto eidx = static_cast<Eigen::Index>(idx);
    total_covered_bases += logrs.ends[eidx] - logrs.starts[eidx];
  }
  const auto segment_length = static_cast<f64>(seg.end - seg.start);
  if (static_cast<f64>(total_covered_bases) < kMinIntervalDensity * segment_length) {
    Logging::Info("Skipping segment {}:{}-{} because interval coverage ({} bp) is below {}% of segment length ({} bp)",
                  seg.contig,
                  seg.start,
                  seg.end,
                  total_covered_bases,
                  static_cast<s32>(kMinIntervalDensity * kPercentMultiplier),
                  static_cast<size_t>(segment_length));
    return true;
  }
  return false;
}

namespace {

// Shared inputs for the logR-only grid build, bundled to keep the helper parameter counts within limits.
// classify_germline is derived once from tumor_only && use_diploid_bias_prior.
struct LogROnlyGridContext {
  vec<GenomicSegment> segments;
  Observations obvs;
  IntervalTrees trees;
  f64 seg_max_logr = 0.0;
  size_t seg_min_num_logrs = 0;
  sex_predict::Sex sex = sex_predict::Sex::kUnknown;
  bool use_diploid_bias_prior = false;
  bool tumor_only = false;
  bool classify_germline = false;
  f64 mean_segment_logr_sd = 0.0;
  SegmentLikelihoodWriter* likelihood_writer = nullptr;
};

// Per-segment inputs for the logR-only grid-cell loop. The logr_obvs reference binds to the segment's filtered logR
// subset, which outlives the cell loop. germline_logr_likelihood/germline_tcn are computed once per segment and reused
// for every cell; only meaningful when classify_germline is true.
struct LogROnlySegmentInputs {
  GenomicSegment seg;
  Eigen::VectorXd logr_obvs;
  f64 mean_segment_logr_sd = 0.0;
  bool use_diploid_bias_prior = false;
  bool expect_haploid = false;
  bool classify_germline = false;
  f64 germline_logr_likelihood = kNegativeInf;
  size_t germline_tcn = 0;
};

// Scratch likelihood buffers reused across grid cells and segments, bundled so the per-segment helpers take a single
// named argument. Members bind to buffers that outlive the grid build. cell_logr is the CN-loop output for the
// somatic-only branch; somatic_logr/somatic_baf and germline_logr are scratch for the germline-classified branch.
struct LogROnlyScratch {
  // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members) - write-through: caller-owned buffers reused across
  // segments and written into per cell; value members would reintroduce a per-segment allocation.
  Eigen::VectorXd& cell_logr;
  Eigen::VectorXd& somatic_logr;
  Eigen::VectorXd& somatic_baf;
  Eigen::VectorXd& germline_logr;
  // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)
};

// Evaluate one grid cell in the tumor-only diploid-bias (germline-classified) pass. Computes the somatic log-ratio
// likelihood with the same additive diploid-bias prior used by the germline helper so the two are comparable, then
// returns the cell value as the max of the germline and somatic logR likelihoods. There is no BAF data in the
// log-ratio-only path, so a somatic cell contributes only its (prior-weighted) logR likelihood. Optionally emits a
// diagnostic row.
f64 EvalLogROnlyGermlineClassifiedCell(const LogROnlySegmentInputs& in,
                                       const f64 purity,
                                       const f64 ploidy,
                                       const LogROnlyScratch& scratch,
                                       SegmentLikelihoodWriter* const likelihood_writer) {
  GenomicSegment somatic_segment;
  somatic_segment = SegmentSerialSummarizedLogRBAFLikelihood(somatic_segment,
                                                             in.logr_obvs,
                                                             purity,
                                                             ploidy,
                                                             in.mean_segment_logr_sd,
                                                             /*mean_baf_sd=*/0.0,
                                                             in.use_diploid_bias_prior,
                                                             in.expect_haploid,
                                                             scratch.somatic_logr,
                                                             scratch.somatic_baf,
                                                             /*baf_likelihood_weight=*/1.0);
  const f64 cell_value = std::max(in.germline_logr_likelihood, somatic_segment.logr_likelihood.value());
  if (likelihood_writer != nullptr) {
    const bool is_germline =
        in.germline_logr_likelihood >= (somatic_segment.logr_likelihood.value() - kGermlineTieMargin);
    likelihood_writer->WriteRow(
        SegmentLikelihoodRowKey{
            .contig = in.seg.contig, .start = in.seg.start, .end = in.seg.end, .purity = purity, .ploidy = ploidy},
        SomaticLikelihoodRow{.logr_likelihood = somatic_segment.logr_likelihood.value(),
                             .total_copy_number = somatic_segment.total_copy_number.value(),
                             .baf_likelihood = std::nullopt,
                             .joint_likelihood = somatic_segment.logr_likelihood,
                             .minor_copy_number = std::nullopt},
        GermlineLikelihoodRow{.logr_likelihood = in.germline_logr_likelihood, .total_copy_number = in.germline_tcn},
        is_germline);
  }
  return cell_value;
}

// Evaluate one grid cell in the somatic-only path (tumor-normal mode or the WGD re-run): take the max over total copy
// number of the (optionally diploid-bias-weighted) logR likelihood. Returns the cell value. There is no BAF in the
// log-ratio-only path, so the somatic BAF/joint/minor and germline/classification columns are emitted as N/A.
f64 EvalLogROnlySomaticCell(const LogROnlySegmentInputs& in,
                            const f64 purity,
                            const f64 ploidy,
                            const LogROnlyScratch& scratch,
                            SegmentLikelihoodWriter* const likelihood_writer) {
  Eigen::VectorXd& logr_likelihoods = scratch.cell_logr;
  for (size_t i = 0; i < static_cast<size_t>(logr_likelihoods.size()); ++i) {
    const auto tcn = static_cast<f64>(i);
    const auto n_cns = static_cast<f64>(logr_likelihoods.size());
    const f64 weight = (in.use_diploid_bias_prior && i == 2) ? 2.0 / n_cns : 1.0 / n_cns;
    logr_likelihoods[static_cast<Eigen::Index>(i)] =
        weight *
        SegmentSummarizedLogRLikelihood(in.logr_obvs, in.mean_segment_logr_sd, purity, ploidy, tcn, in.expect_haploid);
  }
  const size_t somatic_tcn = eigen::IndexMax(logr_likelihoods);
  const f64 cell_value = logr_likelihoods(static_cast<Eigen::Index>(somatic_tcn));
  if (likelihood_writer != nullptr) {
    likelihood_writer->WriteRow(
        SegmentLikelihoodRowKey{
            .contig = in.seg.contig, .start = in.seg.start, .end = in.seg.end, .purity = purity, .ploidy = ploidy},
        SomaticLikelihoodRow{.logr_likelihood = cell_value,
                             .total_copy_number = somatic_tcn,
                             .baf_likelihood = std::nullopt,
                             .joint_likelihood = cell_value,
                             .minor_copy_number = std::nullopt},
        std::nullopt,
        std::nullopt);
  }
  logr_likelihoods.setZero();
  return cell_value;
}

// Fill every (purity, ploidy) cell of one segment's grid, dispatching to the germline-classified or somatic-only cell
// evaluation. Iterates over a stable copy of the grid so the source is neither read nor advanced as the body writes
// into it (cpp:S892, cpp:S5316).
void FillLogROnlyGridCells(const LogROnlySegmentInputs& in,
                           Eigen::MatrixXd& grid,
                           const LogROnlyScratch& scratch,
                           SegmentLikelihoodWriter* const likelihood_writer) {
  // Iterate cells in column-major order (column outer, row inner) so the likelihood writer emits rows in a stable
  // order.
  for (Eigen::Index col_ind = 0; col_ind < grid.cols(); ++col_ind) {
    const f64 ploidy = kLikelihoodPhiMin + (static_cast<f64>(col_ind) * kLikelihoodPhiStep);
    for (Eigen::Index row_ind = 0; row_ind < grid.rows(); ++row_ind) {
      const f64 purity = kLikelihoodRoMin + (static_cast<f64>(row_ind) * kLikelihoodRoStep);
      grid(row_ind, col_ind) = in.classify_germline
                                   ? EvalLogROnlyGermlineClassifiedCell(in, purity, ploidy, scratch, likelihood_writer)
                                   : EvalLogROnlySomaticCell(in, purity, ploidy, scratch, likelihood_writer);
    }
  }
}

// Build the per-segment logR-only likelihood grid for one segment. Segments excluded by mode/sex, allosome, density, or
// quality checks leave their grid unmodified and produce no diagnostic rows. The scratch buffers are reused across
// segments.
void BuildLogROnlySegmentGrid(const LogROnlyGridContext& ctx,
                              const size_t seg_idx,
                              Eigen::MatrixXd& grid,
                              const LogROnlyScratch& scratch) {
  const GenomicSegment& seg = ctx.segments[seg_idx];
  Logging::Info("building grid for segment #{}", seg_idx);
  if (ctx.tumor_only && sex_predict::IsInAllosome(seg.contig)) {
    Logging::Info("skipping segment #{} {}:{}-{} because allosomes are excluded in tumor-only mode",
                  seg_idx,
                  seg.contig,
                  seg.start,
                  seg.end);
    return;
  }
  if (SkipBasedOnModeAndSex(seg, LikelihoodMode::kGermline, ctx.sex)) {
    Logging::Info(
        "skipping segment #{} {}:{}-{} because due to sex of sample and mode", seg_idx, seg.contig, seg.start, seg.end);
    return;
  }
  const bool expect_haploid = ExpectHaploid(seg, LikelihoodMode::kGermline, ctx.sex);
  const auto seg_obvs_idxs = ctx.trees.LookUp(seg.contig, seg.start, seg.end);
  if (seg_obvs_idxs.empty()) {
    Logging::Info("skipping segment #{} {}:{}-{} because it has no observations. Probably a seed segment",
                  seg_idx,
                  seg.contig,
                  seg.start,
                  seg.end);
    return;
  }
  if (HasLowIntervalDensity(seg_obvs_idxs, ctx.obvs, seg)) {
    return;
  }
  if (ShouldSkipSegment(seg, ctx.seg_max_logr, ctx.seg_min_num_logrs, std::nullopt, false)) {
    return;
  }
  const Observations seg_obvs = ctx.obvs.FilterByIdxs(seg_obvs_idxs);
  // Reset the CN-loop scratch for this segment (every cell fully overwrites it before use).
  scratch.cell_logr.setConstant(kNegativeInf);
  // The germline log-ratio likelihood is evaluated at fixed purity/ploidy, so it is computed once per segment and
  // reused for every grid cell. Only needed when germline classification is active.
  f64 germline_logr_likelihood = kNegativeInf;
  // Germline total copy number at the germline-state MLE (CN = 2 excluded). Computed once per segment alongside the
  // germline likelihood and reused for every grid cell; only meaningful when germline classification is active.
  size_t germline_tcn = 0;
  if (ctx.classify_germline) {
    germline_logr_likelihood = SegmentGermlineSummarizedLogRLikelihood(
        seg_obvs.obvs, ctx.mean_segment_logr_sd, expect_haploid, scratch.germline_logr);
    germline_tcn = eigen::IndexMax(scratch.germline_logr);
  }
  const LogROnlySegmentInputs in{.seg = seg,
                                 .logr_obvs = seg_obvs.obvs,
                                 .mean_segment_logr_sd = ctx.mean_segment_logr_sd,
                                 .use_diploid_bias_prior = ctx.use_diploid_bias_prior,
                                 .expect_haploid = expect_haploid,
                                 .classify_germline = ctx.classify_germline,
                                 .germline_logr_likelihood = germline_logr_likelihood,
                                 .germline_tcn = germline_tcn};
  FillLogROnlyGridCells(in, grid, scratch, ctx.likelihood_writer);
}

}  // namespace

/**
 * @brief Builds a grid of LogR likelihoods over each combination of purity and ploidy, where purity and ploidy steps
 * are pre-defined
 * @param segments
 * @param obvs
 * @return Eigen::MatrixXd - max LogR likelihood over all purity-ploidy steps
 */
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
Eigen::MatrixXd GetPurityPloidyGridLogROnly(const vec<GenomicSegment>& segments,
                                            const Observations& obvs,
                                            const PurityPloidyGridParams& params) {
  // iniitialize a 2D grid - purity, ploidy for each segment. Zero-initialized (Eigen::MatrixXd(rows, cols) leaves
  // contents uninitialized) because segments skipped by mode/sex/density/quality leave their grid untouched and
  // must contribute zero to the element-wise sum below.
  vec<Eigen::MatrixXd> grids;
  grids.reserve(segments.size());
  for (size_t i = 0; i < segments.size(); ++i) {
    grids.emplace_back(Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(kLikelihoodRoNumSteps),
                                             static_cast<Eigen::Index>(kLikelihoodPhiNumSteps)));
  }
  // Germline-vs-somatic classification is only active on the diploid-bias (first) pass (see
  // GetPurityPloidyGridSerialLogRBAF for why use_diploid_bias_prior doubles as the WGD-pass signal).
  const bool classify_germline = params.tumor_only && params.use_diploid_bias_prior;
  const IntervalTrees trees(obvs);
  const f64 sd = GetMeanSegmentLogRSD(segments, obvs, trees);
  // Scratch likelihood buffers reused across cells and segments. The germline buffer is only used when germline
  // classification is active; the germline likelihood is independent of the grid cell.
  Eigen::VectorXd cell_logr_likelihoods(kLikelihoodMaxCN1 + 1);
  Eigen::VectorXd somatic_logr_likelihoods(kLikelihoodMaxCN1 + 1);
  Eigen::VectorXd somatic_baf_likelihoods(kLikelihoodMaxCN1 + 1);
  Eigen::VectorXd germline_logr_likelihoods(kLikelihoodMaxCN1 + 1);
  const LogROnlyScratch scratch{.cell_logr = cell_logr_likelihoods,
                                .somatic_logr = somatic_logr_likelihoods,
                                .somatic_baf = somatic_baf_likelihoods,
                                .germline_logr = germline_logr_likelihoods};
  const LogROnlyGridContext ctx{.segments = segments,
                                .obvs = obvs,
                                .trees = trees,
                                .seg_max_logr = params.seg_max_logr,
                                .seg_min_num_logrs = params.seg_min_num_logrs,
                                .sex = params.sex,
                                .use_diploid_bias_prior = params.use_diploid_bias_prior,
                                .tumor_only = params.tumor_only,
                                .classify_germline = classify_germline,
                                .mean_segment_logr_sd = sd,
                                .likelihood_writer = params.likelihood_writer};
  for (size_t seg_idx = 0; seg_idx < segments.size(); ++seg_idx) {
    BuildLogROnlySegmentGrid(ctx, seg_idx, grids[seg_idx], scratch);
  }
  // Sum the per-segment grids element-wise into one grid.
  Eigen::MatrixXd sum_of_grids = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(kLikelihoodRoNumSteps),
                                                       static_cast<Eigen::Index>(kLikelihoodPhiNumSteps));
  for (const auto& grid : grids) {
    sum_of_grids += grid;
  }
  return sum_of_grids;
}

namespace {

// Per-segment inputs for the grid-cell evaluation loop, bundled to keep the helper parameter counts within limits.
struct SegmentGridInputs {
  GenomicSegment seg;
  Eigen::VectorXd logr_obvs;
  std::optional<f64> seg_mbaf;
  f64 mean_segment_logr_sd = 0.0;
  f64 mean_segment_mbaf_sd = 0.0;
  bool use_diploid_bias_prior = false;
  bool expect_haploid = false;
  bool classify_germline = false;
  f64 germline_logr_likelihood = 0.0;
  size_t germline_tcn = 0;
  f64 baf_likelihood_weight = 1.0;
};

// Mutable output grids for one segment (joint, logR, BAF), bundled to keep the helper parameter counts within limits.
struct SegmentGrids {
  // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members) - write-through: output grids that the helper writes
  // into caller-owned per-segment storage; value members would break the write-back.
  Eigen::MatrixXd& joint;
  Eigen::MatrixXd& logr;
  Eigen::MatrixXd& baf;
  // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)
};

// Per-segment likelihood grids for every metric (one grid per segment, per metric). Bundled so the sum helper takes a
// single named argument instead of three same-typed vectors (cpp:S5419).
struct PerSegmentGridSet {
  vec<Eigen::MatrixXd> joint;
  vec<Eigen::MatrixXd> logr;
  vec<Eigen::MatrixXd> baf;
};

// Scratch likelihood buffers reused across grid cells and segments, bundled so the per-segment helpers take a single
// named argument instead of several same-typed Eigen::VectorXd references (cpp:S5419). Members bind to buffers that
// outlive the grid build. germline_logr is only populated when germline classification is active.
struct GridScratchBuffers {
  // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members) - write-through: caller-owned buffers reused across
  // segments and written into per cell; value members would reintroduce a per-segment allocation.
  Eigen::VectorXd& logr;
  Eigen::VectorXd& baf;
  Eigen::VectorXd& germline_logr;
  // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)
};

// Emit one diagnostic row for the current grid cell. The germline columns are populated only on the diploid-bias pass
// (classify_germline); on the uniform pass they are emitted as N/A. The somatic BAF and minor columns are N/A when the
// segment has no usable BAF; the joint column equals the somatic logR likelihood in that case, since only the log-ratio
// is available to the search.
void WriteSegmentGridCellRow(SegmentLikelihoodWriter& writer,
                             const SegmentGridInputs& in,
                             const GenomicSegment& segment,
                             const f64 purity,
                             const f64 ploidy,
                             const bool is_germline) {
  const std::optional<GermlineLikelihoodRow> germline_row =
      in.classify_germline ? std::optional<GermlineLikelihoodRow>(GermlineLikelihoodRow{
                                 .logr_likelihood = in.germline_logr_likelihood, .total_copy_number = in.germline_tcn})
                           : std::nullopt;
  const std::optional<bool> classification = in.classify_germline ? std::optional<bool>(is_germline) : std::nullopt;
  writer.WriteRow(
      SegmentLikelihoodRowKey{
          .contig = in.seg.contig, .start = in.seg.start, .end = in.seg.end, .purity = purity, .ploidy = ploidy},
      SomaticLikelihoodRow{.logr_likelihood = segment.logr_likelihood.value(),
                           .total_copy_number = segment.total_copy_number.value(),
                           .baf_likelihood = segment.baf_likelihood,
                           .joint_likelihood = segment.joint_likelihood,
                           .minor_copy_number = segment.minor_copy_number},
      germline_row,
      classification);
}

// Loop-invariant inputs shared across every segment of the serial logR-BAF grid build, bundled to keep the per-segment
// helper parameter count within limits.
struct GridBuildContext {
  vec<GenomicSegment> segments;
  Observations logrs;
  Observations ref_ads;
  Observations alt_ads;
  IntervalTrees logr_trees;
  IntervalTrees variant_trees;
  f64 seg_max_logr;
  size_t seg_min_num_logrs;
  size_t seg_min_num_snps;
  sex_predict::Sex sex;
  bool use_diploid_bias_prior;
  bool tumor_only;
  bool classify_germline;
  f64 mean_segment_logr_sd;
  f64 mean_segment_mbaf_sd;
  f64 baf_likelihood_weight;
  SegmentLikelihoodWriter* likelihood_writer;
};

// Sum the per-segment likelihood grids element-wise into a single grid per metric. NaN BAF entries (segments with no
// usable BAF) are treated as a zero contribution before summation.
std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> SumPerSegmentGrids(PerSegmentGridSet& grids) {
  const auto rows = static_cast<Eigen::Index>(kLikelihoodRoNumSteps);
  const auto cols = static_cast<Eigen::Index>(kLikelihoodPhiNumSteps);
  Eigen::MatrixXd sum_joint_ll_grids = Eigen::MatrixXd::Zero(rows, cols);
  Eigen::MatrixXd sum_logr_ll_grids = Eigen::MatrixXd::Zero(rows, cols);
  Eigen::MatrixXd sum_baf_ll_grids = Eigen::MatrixXd::Zero(rows, cols);
  for (const auto& joint_ll_grid : grids.joint) {
    sum_joint_ll_grids += joint_ll_grid;
  }
  for (const auto& logr_ll_grid : grids.logr) {
    sum_logr_ll_grids += logr_ll_grid;
  }
  for (auto& baf_ll_grid : grids.baf) {
    eigen::ReplaceNaN(baf_ll_grid, 0);
    sum_baf_ll_grids += baf_ll_grid;
  }
  return {sum_joint_ll_grids, sum_logr_ll_grids, sum_baf_ll_grids};
}

// Fill every (purity, ploidy) cell of one segment's grids with its joint/logR/BAF likelihoods, classifying each cell as
// germline or somatic and optionally emitting a diagnostic row. The scratch likelihood buffers are reused across cells.
void FillSegmentGridCells(const SegmentGridInputs& in,
                          SegmentGrids grids,
                          const GridScratchBuffers& scratch,
                          SegmentLikelihoodWriter* const likelihood_writer) {
  // Iterate cells in column-major order (column outer, row inner) so the likelihood writer emits rows in a stable
  // order. The grid dimensions are fixed (kLikelihoodRoNumSteps x kLikelihoodPhiNumSteps), so reading grids.joint's
  // extents while the loop body writes into it is safe and needs no defensive copy.
  for (Eigen::Index col_ind = 0; col_ind < grids.joint.cols(); ++col_ind) {
    const f64 ploidy = kLikelihoodPhiMin + (static_cast<f64>(col_ind) * kLikelihoodPhiStep);
    for (Eigen::Index row_ind = 0; row_ind < grids.joint.rows(); ++row_ind) {
      const f64 purity = kLikelihoodRoMin + (static_cast<f64>(row_ind) * kLikelihoodRoStep);
      GenomicSegment segment;
      segment.mbaf = in.seg_mbaf;
      segment = SegmentSerialSummarizedLogRBAFLikelihood(segment,
                                                         in.logr_obvs,
                                                         purity,
                                                         ploidy,
                                                         in.mean_segment_logr_sd,
                                                         in.mean_segment_mbaf_sd,
                                                         in.use_diploid_bias_prior,
                                                         in.expect_haploid,
                                                         scratch.logr,
                                                         scratch.baf,
                                                         in.baf_likelihood_weight);
      // Classify the segment as germline or somatic using a log-ratio-only comparison (BAF is uninformative for
      // germline events). A germline segment contributes only its germline logR likelihood (BAF set to 0, joint =
      // germline logR); a somatic segment contributes its somatic logR, BAF, and joint likelihoods.
      const bool is_germline = in.classify_germline &&
                               (in.germline_logr_likelihood >= (segment.logr_likelihood.value() - kGermlineTieMargin));
      if (is_germline) {
        grids.logr(row_ind, col_ind) = in.germline_logr_likelihood;
        grids.baf(row_ind, col_ind) = 0.0;
        grids.joint(row_ind, col_ind) = in.germline_logr_likelihood;
      } else {
        grids.joint(row_ind, col_ind) = segment.joint_likelihood.value();
        grids.logr(row_ind, col_ind) = segment.logr_likelihood.value();
        grids.baf(row_ind, col_ind) = segment.baf_likelihood.value_or(NAN);
      }
      if (likelihood_writer != nullptr) {
        WriteSegmentGridCellRow(*likelihood_writer, in, segment, purity, ploidy, is_germline);
      }
    }
  }
}

// Build the per-segment likelihood grids for one segment, writing into the segment's slot of the output grids. Segments
// excluded by mode/sex, allosome, density, or quality checks leave their grids unmodified and produce no diagnostic
// rows. The scratch likelihood buffers are reused across segments.
void BuildSegmentGrid(const GridBuildContext& ctx,
                      const size_t seg_idx,
                      SegmentGrids grids,
                      const GridScratchBuffers& scratch) {
  const GenomicSegment& seg = ctx.segments[seg_idx];
  Logging::Info("building grid for segment #{}", seg_idx);
  // In tumor-only mode the allosome log ratios come from a median normalization anchored on the autosomes, so the
  // somatic ploidy policy (SkipBasedOnModeAndSex with kSomatic) governs which allosome segments are usable: chrX is
  // retained for female samples and chrX pseudo-autosomal regions for male samples, while chrY and single-copy non-PAR
  // chrX are excluded. Allosomes are dropped entirely when the sex is unknown.
  if (SkipBasedOnModeAndSex(seg, LikelihoodMode::kSomatic, ctx.sex)) {
    Logging::Info(
        "skipping segment #{} {}:{}-{} because due to sex of sample and mode", seg_idx, seg.contig, seg.start, seg.end);
    return;
  }
  const bool expect_haploid = ExpectHaploid(seg, LikelihoodMode::kSomatic, ctx.sex);
  const vec<size_t> logr_idxs = ctx.logr_trees.LookUp(seg.contig, seg.start, seg.end);
  if (logr_idxs.empty()) {
    Logging::Info(
        "skipping segment #{} {}:{}-{} because it has no logr observations", seg_idx, seg.contig, seg.start, seg.end);
    return;
  }
  if (HasLowIntervalDensity(logr_idxs, ctx.logrs, seg)) {
    return;
  }
  if (ShouldSkipSegment(seg, ctx.seg_max_logr, ctx.seg_min_num_logrs, ctx.seg_min_num_snps)) {
    return;
  }

  const Observations logr_subset = ctx.logrs.FilterByIdxs(logr_idxs);
  const vec<size_t> variant_idxs = ctx.variant_trees.LookUp(seg.contig, seg.start, seg.end);
  const Observations ref_ad_subset = ctx.ref_ads.FilterByIdxs(variant_idxs);
  const Observations alt_ad_subset = ctx.alt_ads.FilterByIdxs(variant_idxs);
  std::optional<f64> seg_mbaf = std::nullopt;
  if (ref_ad_subset.obvs.size() != 0) {
    const f64 mbaf = GetPeakOfBAFDistributionThenMirror(ref_ad_subset.obvs, alt_ad_subset.obvs);
    seg_mbaf = std::isnan(mbaf) ? std::nullopt : std::optional<f64>(mbaf);
  }
  // The germline log-ratio likelihood is evaluated at fixed purity/ploidy, so it is computed once per segment and
  // reused for every grid cell. Only needed when germline classification is active.
  f64 germline_logr_likelihood = kNegativeInf;
  // Germline total copy number at the germline-state MLE (CN = 2 excluded). Computed once per segment alongside the
  // germline likelihood and reused for every grid cell; only meaningful when germline classification is active.
  size_t germline_tcn = 0;
  if (ctx.classify_germline) {
    germline_logr_likelihood = SegmentGermlineSummarizedLogRLikelihood(
        logr_subset.obvs, ctx.mean_segment_logr_sd, expect_haploid, scratch.germline_logr);
    germline_tcn = eigen::IndexMax(scratch.germline_logr);
  }
  const SegmentGridInputs grid_inputs{.seg = seg,
                                      .logr_obvs = logr_subset.obvs,
                                      .seg_mbaf = seg_mbaf,
                                      .mean_segment_logr_sd = ctx.mean_segment_logr_sd,
                                      .mean_segment_mbaf_sd = ctx.mean_segment_mbaf_sd,
                                      .use_diploid_bias_prior = ctx.use_diploid_bias_prior,
                                      .expect_haploid = expect_haploid,
                                      .classify_germline = ctx.classify_germline,
                                      .germline_logr_likelihood = germline_logr_likelihood,
                                      .germline_tcn = germline_tcn,
                                      .baf_likelihood_weight = ctx.baf_likelihood_weight};
  FillSegmentGridCells(grid_inputs, grids, scratch, ctx.likelihood_writer);
}

}  // namespace

/**
 * @brief Builds a grid of joint LogR-BAF likelihoods over each combination of purity and ploidy, where purity and
 * ploidy steps are pre-defined
 * @param segments
 * @param obvs
 * @return matrices of joint-likelihoods, logrs, bafs for each purity-ploidy combination
 */
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> GetPurityPloidyGridSerialLogRBAF(
    const vec<GenomicSegment>& segments,
    const Observations& logrs,  // NOSONAR(cpp:S5419) adjacent Observations params
    const Observations& ref_ads,
    const Observations& alt_ads,
    const PurityPloidyGridParams& params) {
  // Germline-vs-somatic classification is only active on the diploid-bias (first) pass of the purity/ploidy search.
  // use_diploid_bias_prior doubles as the "which WGD pass" signal: the first pass runs with the prior on, the WGD
  // re-run runs with it off. ANDing on it therefore disables germline classification during the re-run, where the
  // germline ploidy = 2 assumption is invalid, without threading a separate flag through the call stack.
  const bool classify_germline = params.tumor_only && params.use_diploid_bias_prior;
  // iniitialize a 2D grid - purity, ploidy for each segment
  PerSegmentGridSet per_segment_grids;
  per_segment_grids.joint.reserve(segments.size());
  per_segment_grids.logr.reserve(segments.size());
  per_segment_grids.baf.reserve(segments.size());
  // Zero-initialize each per-segment grid (Eigen::MatrixXd(rows, cols) leaves contents uninitialized) because segments
  // excluded by mode/sex/density/quality leave their grids untouched and must contribute zero to the element-wise sums
  // in SumPerSegmentGrids.
  const auto grid_rows = static_cast<Eigen::Index>(kLikelihoodRoNumSteps);
  const auto grid_cols = static_cast<Eigen::Index>(kLikelihoodPhiNumSteps);
  for (size_t i = 0; i < segments.size(); ++i) {
    per_segment_grids.joint.emplace_back(Eigen::MatrixXd::Zero(grid_rows, grid_cols));
    per_segment_grids.logr.emplace_back(Eigen::MatrixXd::Zero(grid_rows, grid_cols));
    per_segment_grids.baf.emplace_back(Eigen::MatrixXd::Zero(grid_rows, grid_cols));
  }
  const IntervalTrees logr_trees(logrs);
  IntervalTrees variant_trees(ref_ads);
  Eigen::VectorXd logr_likelihoods(kLikelihoodMaxCN1 + 1);
  Eigen::VectorXd baf_likelihoods(kLikelihoodMaxCN1 + 1);
  // Germline log-ratio likelihood is evaluated at fixed purity/ploidy and is therefore independent of the grid cell.
  // It is reused for every (purity, ploidy) cell of a segment when germline classification is active.
  Eigen::VectorXd germline_logr_likelihoods(kLikelihoodMaxCN1 + 1);
  const f64 mean_segment_logr_sd = GetMeanSegmentLogRSD(segments, logrs, logr_trees);
  // The mirror process reduces the variances of the BAFs. We f64 the mean segment BAF SD to account for this
  // we also w
  const f64 mean_segment_mbaf_sd = kMirrorMBAFSDScale * GetMeanSegmentMBAFSD(segments, ref_ads, alt_ads, variant_trees);
  const GridBuildContext ctx{.segments = segments,
                             .logrs = logrs,
                             .ref_ads = ref_ads,
                             .alt_ads = alt_ads,
                             .logr_trees = logr_trees,
                             .variant_trees = variant_trees,
                             .seg_max_logr = params.seg_max_logr,
                             .seg_min_num_logrs = params.seg_min_num_logrs,
                             .seg_min_num_snps = params.seg_min_num_snps,
                             .sex = params.sex,
                             .use_diploid_bias_prior = params.use_diploid_bias_prior,
                             .tumor_only = params.tumor_only,
                             .classify_germline = classify_germline,
                             .mean_segment_logr_sd = mean_segment_logr_sd,
                             .mean_segment_mbaf_sd = mean_segment_mbaf_sd,
                             .baf_likelihood_weight = params.baf_likelihood_weight,
                             .likelihood_writer = params.likelihood_writer};
  const GridScratchBuffers scratch{
      .logr = logr_likelihoods, .baf = baf_likelihoods, .germline_logr = germline_logr_likelihoods};
  for (size_t seg_idx = 0; seg_idx < segments.size(); ++seg_idx) {
    BuildSegmentGrid(ctx,
                     seg_idx,
                     SegmentGrids{.joint = per_segment_grids.joint[seg_idx],
                                  .logr = per_segment_grids.logr[seg_idx],
                                  .baf = per_segment_grids.baf[seg_idx]},
                     scratch);
  }
  // Sum the per-segment grids element-wise into one grid per metric.
  return SumPerSegmentGrids(per_segment_grids);
}

/**
 * @brief given a purity ploidy grid, return the purity ploidy associated with the absolute maximum of the
 * grid
 * @param purity_ploidy_grid
 * @return
 */
std::tuple<f64, f64> EstimatePurityPloidyByTotalMaximum(const Eigen::MatrixXd& purity_ploidy_grid) {
  const auto [row_max_ll, col_max_ll] = eigen::ArgMaxRowCol(purity_ploidy_grid);
  const f64 purity = kLikelihoodRoMin + (kLikelihoodRoStep * static_cast<f64>(row_max_ll));
  const f64 ploidy = kLikelihoodPhiMin + (kLikelihoodPhiStep * static_cast<f64>(col_max_ll));
  return {purity, ploidy};
}

/**
 * @brief write joint, logr and BAF likelihoods from the purity ploidy grid to disk
 * @param purity_ploidy_grid
 * @param out_fname
 */
void WritePurityPloidyGrid(const Eigen::MatrixXd& purity_ploidy_grid_joint_lls,
                           const std::optional<Eigen::MatrixXd>& purity_ploidy_grid_logr_lls,
                           const std::optional<Eigen::MatrixXd>& purity_ploidy_grid_baf_lls,
                           const fs::path& out_fname,
                           const io::CommandLineInfo& command_line_info,
                           std::optional<TotalCopyNumberPrior> prior) {
  std::ofstream ofs(out_fname);
  if (!ofs.is_open()) {
    Logging::Error("Failed to open purity ploidy grid output file: {}", out_fname.string());
    throw std::runtime_error("Failed to open purity ploidy grid output file: " + out_fname.string());
  }
  if (!command_line_info.version.empty()) {
    io::WriteTsvMetadata(ofs, command_line_info);
  }
  if (prior.has_value()) {
    ofs << "##TotalCopyNumberPrior=" << enum_util::FormatEnumName(prior.value()) << "\n";
  }
  ofs << kPurityPloidyGridOutPurity << "\t" << kPurityPloidyGridOutPloidy << "\t"
      << kPurityPloidyGridOutTotalJointLikelihood << "\t" << kPurityPloidyGridOutTotalLogRatioLikelihood << "\t"
      << kPurityPloidyGridOutTotalBAlleleFrequencyLikelihood << '\n';
  f64 purity_val = NAN;
  f64 ploidy_val = NAN;
  for (Eigen::Index row_idx = 0; row_idx < purity_ploidy_grid_joint_lls.rows(); ++row_idx) {
    // Determine corresponding purity value of the row index
    purity_val = kLikelihoodRoMin + (static_cast<f64>(row_idx) * kLikelihoodRoStep);
    for (Eigen::Index col_idx = 0; col_idx < purity_ploidy_grid_joint_lls.cols(); ++col_idx) {
      // Determine corresponding ploidy value of the column index
      ploidy_val = kLikelihoodPhiMin + static_cast<f64>(col_idx) * kLikelihoodPhiStep;
      ofs << FloatAsString(purity_val, kPurityPloidyOutputPrecision) << "\t"
          << FloatAsString(ploidy_val, kPurityPloidyOutputPrecision) << "\t"
          << FloatAsString(purity_ploidy_grid_joint_lls(row_idx, col_idx));
      if (purity_ploidy_grid_logr_lls.has_value()) {
        ofs << "\t" << FloatAsString(purity_ploidy_grid_logr_lls.value()(row_idx, col_idx));
      } else {
        ofs << "\t" << util::kNA;
      }
      if (purity_ploidy_grid_baf_lls.has_value()) {
        ofs << "\t" << FloatAsString(purity_ploidy_grid_baf_lls.value()(row_idx, col_idx));
      } else {
        ofs << "\t" << util::kNA;
      }
      ofs << '\n';
    }
  }
}

/**
 * @brief write purity and ploidy to disk
 * @param purity
 * @param ploidy
 * @param out_fname
 */
void WritePurityPloidy(const f64 purity,
                       const f64 ploidy,
                       const fs::path& out_fname,
                       const io::CommandLineInfo& command_line_info,
                       std::optional<TotalCopyNumberPrior> prior) {
  std::ofstream ofs(out_fname);
  if (!ofs.is_open()) {
    Logging::Error("Failed to open purity ploidy output file: {}", out_fname.string());
    throw std::runtime_error("Failed to open purity ploidy output file: " + out_fname.string());
  }
  if (!command_line_info.version.empty()) {
    io::WriteTsvMetadata(ofs, command_line_info);
  }
  if (prior.has_value()) {
    ofs << "##TotalCopyNumberPrior=" << enum_util::FormatEnumName(prior.value()) << "\n";
  }
  ofs << kPurityPloidyOutPurity << "\t" << kPurityPloidyOutPloidy << '\n';
  // Match the purity/ploidy precision used by the segment, metrics, and VCF outputs.
  ofs << FloatAsString(purity, kPurityPloidyOutputPrecision) << "\t"
      << FloatAsString(ploidy, kPurityPloidyOutputPrecision) << "\n";
}

namespace {

/**
 * @brief Returns the subset of segments classified as somatic at a fixed purity/ploidy solution (tumor-only mode).
 *
 * Each autosomal segment is classified as germline or somatic using a log-ratio-only comparison: the germline
 * log-ratio likelihood (evaluated at the fixed germline purity/ploidy) is compared against the somatic log-ratio
 * likelihood evaluated at the supplied @p purity and @p ploidy. Segments whose germline likelihood exceeds their
 * somatic likelihood are dropped so that only somatic segments feed the downstream WGD QC step. Allosomes are
 * excluded entirely, regardless of sex, because WGD QC relies on the most reliable diploid segments.
 *
 * @param segments All genomic segments for the sample.
 * @param logrs Per-locus log-ratio observations.
 * @param purity Top-solution tumor purity from the diploid-bias search.
 * @param ploidy Top-solution tumor ploidy from the diploid-bias search.
 * @param sex Sex of the sample.
 * @return The somatic subset of @p segments.
 */
vec<GenomicSegment> FilterSomaticSegmentsForWGDQC(const vec<GenomicSegment>& segments,
                                                  const Observations& logrs,
                                                  const f64 purity,
                                                  const f64 ploidy,
                                                  const sex_predict::Sex sex) {
  const IntervalTrees logr_trees(logrs);
  const f64 mean_segment_logr_sd = GetMeanSegmentLogRSD(segments, logrs, logr_trees);
  Eigen::VectorXd germline_logr_likelihoods(kLikelihoodMaxCN1 + 1);
  Eigen::VectorXd somatic_logr_likelihoods(kLikelihoodMaxCN1 + 1);
  Eigen::VectorXd somatic_baf_likelihoods(kLikelihoodMaxCN1 + 1);
  vec<GenomicSegment> somatic_segments;
  somatic_segments.reserve(segments.size());
  // This helper only runs in tumor-only mode. WGD QC relies on the most reliable diploid segments, so all allosomes are
  // excluded regardless of sex. Remaining (autosomal) segments are then classified as germline or somatic using a
  // log-ratio-only likelihood comparison at the supplied purity/ploidy solution.
  for (const auto& seg : segments) {
    if (sex_predict::IsInAllosome(seg.contig)) {
      continue;
    }
    if (SkipBasedOnModeAndSex(seg, LikelihoodMode::kSomatic, sex)) {
      somatic_segments.push_back(seg);
      continue;
    }
    const vec<size_t> logr_idxs = logr_trees.LookUp(seg.contig, seg.start, seg.end);
    const Observations logr_subset = logrs.FilterByIdxs(logr_idxs);
    if (logr_subset.obvs.size() == 0) {
      somatic_segments.push_back(seg);
      continue;
    }
    const bool expect_haploid = ExpectHaploid(seg, LikelihoodMode::kSomatic, sex);
    const f64 germline_logr_likelihood = SegmentGermlineSummarizedLogRLikelihood(
        logr_subset.obvs, mean_segment_logr_sd, expect_haploid, germline_logr_likelihoods);
    GenomicSegment somatic_segment;
    somatic_segment = SegmentSerialSummarizedLogRBAFLikelihood(somatic_segment,
                                                               logr_subset.obvs,
                                                               purity,
                                                               ploidy,
                                                               mean_segment_logr_sd,
                                                               /*mean_baf_sd=*/0.0,
                                                               /*prior_2n=*/true,
                                                               expect_haploid,
                                                               somatic_logr_likelihoods,
                                                               somatic_baf_likelihoods,
                                                               /*baf_likelihood_weight=*/1.0);
    if (germline_logr_likelihood >= (somatic_segment.logr_likelihood.value() - kGermlineTieMargin)) {
      Logging::Info("WGD QC: dropping germline-classified segment {}:{}-{} from WGD QC (tumor-only mode)",
                    seg.contig,
                    seg.start,
                    seg.end);
      continue;
    }
    somatic_segments.push_back(seg);
  }
  return somatic_segments;
}

/**
 * @brief Identifies candidate clonal heterozygous deletion segments for WGD QC.
 *
 * Filters segments by: autosome membership, minimum logR/SNP observations, interval density,
 * and expected total copy number in the het-del range.
 */
vec<std::pair<const GenomicSegment*, f64>> FindCandidateHetDelSegments(const vec<GenomicSegment>& segments,
                                                                       const Observations& logrs,
                                                                       const RefAltObservations& ref_alt_depths,
                                                                       const f64 purity,
                                                                       const f64 ploidy,
                                                                       const PurityPloidySearchOptions& options) {
  const IntervalTrees logr_trees(logrs);
  const IntervalTrees variant_trees(ref_alt_depths.ref_obvs);
  vec<std::pair<const GenomicSegment*, f64>> candidates;
  for (const auto& seg : segments) {
    // WGD QC relies on the most reliable diploid segments, so all allosomes are excluded regardless of sex: their log
    // ratios come from an autosome-anchored normalization and are not dependable for the het-del mBAF analysis.
    if (sex_predict::IsInAllosome(seg.contig) || !seg.mean_logr.has_value()) {
      continue;
    }
    if (!seg.num_obs.has_value() || seg.num_obs.value() < options.seg_min_num_logrs) {
      continue;
    }
    if (!seg.num_snps.has_value() || seg.num_snps.value() < options.seg_min_num_snps) {
      continue;
    }
    const vec<size_t> logr_idxs = logr_trees.LookUp(seg.contig, seg.start, seg.end);
    if (HasLowIntervalDensity(logr_idxs, logrs, seg)) {
      continue;
    }
    const f64 expected_tcn = stats::ExpectedTotalCopyNumber(purity, ploidy, seg.mean_logr.value(), false);
    if (expected_tcn < kHetDelExpectedTcnMin || expected_tcn > kHetDelExpectedTcnMax) {
      continue;
    }
    const vec<size_t> variant_idxs = variant_trees.LookUp(seg.contig, seg.start, seg.end);
    if (variant_idxs.empty()) {
      continue;
    }
    const Observations ref_subset = ref_alt_depths.ref_obvs.FilterByIdxs(variant_idxs);
    const Observations alt_subset = ref_alt_depths.alt_obvs.FilterByIdxs(variant_idxs);
    const f64 seg_mbaf = GetPeakOfBAFDistributionThenMirror(ref_subset.obvs, alt_subset.obvs);
    if (std::isnan(seg_mbaf)) {
      continue;
    }
    if (seg_mbaf >= kHetDelMaxMBAF) {
      Logging::Info("WGD QC: excluding candidate segment {}:{}-{} because mBAF {:.4f} > {:.2f}",
                    seg.contig,
                    seg.start,
                    seg.end,
                    seg_mbaf,
                    kHetDelMaxMBAF);
      continue;
    }
    candidates.emplace_back(&seg, seg_mbaf);
  }
  return candidates;
}

/**
 * @brief Fraction of the assessable autosomal genome (by length) that is loss-of-heterozygosity at a fixed solution.
 *
 * At the given purity/ploidy every segment is assigned a total and minor copy number by the same production
 * call-assignment used by predict-somatic-cnv (AllSegmentsSerialLogRBAFLikelihoods in somatic mode): the full 0-20
 * copy-number range, sex-aware allosome handling, and per-segment haploid expectation. A segment is counted as
 * loss-of-heterozygosity (LOH) when its inferred minor copy number is 0. The returned value is the total length of LOH
 * segments divided by the total length of all assessable autosomal segments.
 *
 * The fraction is restricted to autosomes: allosome log ratios come from an autosome-anchored normalization and are not
 * dependable for this test. A segment is "assessable" only when the somatic model produced a joint call with a BAF
 * term (baf_likelihood set) — this excludes allosome segments the model skips (which carry a sentinel minor copy
 * number of 0) and segments with no overlapping BAF variants (which have no minor copy number), so neither inflates
 * the denominator or the LOH numerator.
 *
 * @param segments All genomic segments for the sample.
 * @param logrs Per-locus log-ratio observations.
 * @param ref_alt_depths Per-locus reference/alternate allelic depths (BAF support).
 * @param purity Top-solution tumor purity from the diploid-bias search.
 * @param ploidy Top-solution tumor ploidy from the diploid-bias search.
 * @param sex Sample sex (drives allosome handling in the somatic call assignment).
 * @return LOH fraction in [0, 1], or 0.0 when no segment can be assessed.
 */
f64 GenomeFractionWithLOH(const vec<GenomicSegment>& segments,
                          const Observations& logrs,
                          const RefAltObservations& ref_alt_depths,
                          const f64 purity,
                          const f64 ploidy,
                          const sex_predict::Sex sex) {
  vec<GenomicSegment> scored = segments;
  try {
    // Assign total/minor copy number to every segment exactly as predict-somatic-cnv does, so the LOH fraction reflects
    // the calls the sample would actually receive at this solution.
    AllSegmentsSerialLogRBAFLikelihoods(
        scored, logrs, ref_alt_depths.ref_obvs, ref_alt_depths.alt_obvs, purity, ploidy, LikelihoodMode::kSomatic, sex);
  } catch (const std::runtime_error& e) {  // NOSONAR(cpp:S1181)
    // GetMeanSegment{LogR,MBAF}SD signal insufficient log-ratio/BAF support via error::Error, which returns a plain
    // std::runtime_error, so std::runtime_error is the narrowest type available to catch here. A genome-wide LOH
    // fraction cannot be assessed in that case, so treat it as no detectable LOH.
    Logging::Info("WGD QC: cannot assess genome-wide LOH (insufficient log-ratio/BAF support): {}", e.what());
    return 0.0;
  }
  f64 loh_length = 0;
  f64 assessable_length = 0;
  for (const auto& seg : scored) {
    // Restrict to autosomes with a genuine joint call: a set BAF likelihood means the somatic model both kept the
    // segment (not a skipped allosome) and had BAF variants to call a minor copy number from.
    if (sex_predict::IsInAllosome(seg.contig) || !seg.baf_likelihood.has_value() ||
        !seg.minor_copy_number.has_value()) {
      continue;
    }
    const auto seg_length = static_cast<f64>(seg.end - seg.start);
    assessable_length += seg_length;
    if (seg.minor_copy_number.value() == 0) {
      loh_length += seg_length;
    }
  }
  if (assessable_length <= 0) {
    return 0.0;
  }
  return loh_length / assessable_length;
}

}  // namespace

/**
 * @brief Checks whether purity/ploidy solution is valid using a diploid bias prior
 *
 * The algorithm may converge on ploidy~2 even though the real ploidy is ~4 (whole genome doubling, WGD) when using a
 * diploid bias prior for total copy number. This function detects that scenario by inspecting the mirrored BAF
 * distribution of putative clonal heterozygous-deletion segments (expected total copy number ≈ 1).
 *
 * In a true WGD tumour those segments are actually copy-neutral LOH (TCN≈2), so their mBAF values cluster near 0 and
 * near 0.5 (assuming 100% purity) simultaneously, producing two distinct peaks in the KDE. A non-WGD tumour produces a
 * single peak for these deletions. If two peaks separated by more than @p kWGDMaxPeakDifference are detected the
 * function returns TotalCopyNumberPrior::kUniform, signalling that the search should be re-run with a uniform
 * copy-number prior.
 *
 * The mBAF-peak test needs enough tumor signal to resolve the two peaks, so it is skipped when purity is below
 * kMinPurityToEvaluateWgdQc. In that low-purity regime a fallback still detects WGD: if at least
 * kMinLohFractionForWgd of the assessable autosomal genome (by length) is in loss of heterozygosity, the
 * diploid-bias solution is rejected in favour of a uniform prior. Genome-wide LOH at this
 * scale is not expected in non-WGD tumours. The fallback itself has a lower purity floor (kMinPurityForLohFallback):
 * below it even the LOH signal is too unreliable to distinguish WGD from noise, so the diploid-bias solution is
 * accepted outright.
 *
 * Early-exit conditions that accept the diploid-bias solution without further testing:
 *   - purity <= kMinPurityForLohFallback
 *   - kMinPurityForLohFallback < purity < kMinPurityToEvaluateWgdQc AND genome-wide LOH fraction <
 *     kMinLohFractionForWgd
 *   - fewer than kMinCandidateHetDelSegments segments pass all filters
 *   - fewer than kMinCandidateHetDelSegments of those segments yield a valid mBAF estimate
 *   - the KDE of candidate-segment mBAFs has fewer than two local maxima
 * @param segments All genomic segments for the sample.
 * @param ref_alt_depths Per-locus reference and alternate allele depths (used to look up BAF values).
 * @param purity Estimated tumor purity from the purity/ploidy search.
 * @param ploidy Estimated tumor ploidy from the purity/ploidy search.
 * @param options Purity/ploidy search options (segment filter thresholds).
 * @param sex Sample sex (drives allosome handling in the genome-wide LOH assignment).
 * @return TotalCopyNumberPrior::kDiploidBias if the diploid-bias solution is accepted,
 *         TotalCopyNumberPrior::kUniform if it is rejected (WGD suspected).
 */
// NOLINTNEXTLINE(misc-use-internal-linkage): intentionally external, tested via forward declaration
TotalCopyNumberPrior EvaluatePurityPloidySolution(
    const vec<GenomicSegment>& segments,  // NOSONAR(S831): intentionally internal, tested via forward declaration
    const Observations& logrs,
    const RefAltObservations& ref_alt_depths,
    const f64 purity,
    const f64 ploidy,
    const PurityPloidySearchOptions& options,
    const sex_predict::Sex sex) {
  if (purity < kMinPurityToEvaluateWgdQc) {
    // Purity is too low for the mBAF-peak test to resolve the two peaks. Fall back to a genome-wide LOH check: a genome
    // dominated by loss-of-heterozygosity segments is a signature of whole-genome doubling.
    if (purity <= kMinPurityForLohFallback) {
      // Below this purity even the LOH signal is too unreliable to distinguish WGD from noise; accept the
      // diploid-bias solution without running the LOH check as a safety measure against spurious WGD calls.
      Logging::Info(
          "WGD QC: purity {:.2f} is below the LOH-fallback floor {:.2f}; accepting diploid-bias solution without the "
          "genome-wide LOH check",
          purity,
          kMinPurityForLohFallback);
      return TotalCopyNumberPrior::kDiploidBias;
    }
    const f64 loh_fraction = GenomeFractionWithLOH(segments, logrs, ref_alt_depths, purity, ploidy, sex);
    if (loh_fraction >= options.min_loh_fraction_for_wgd) {
      Logging::Info(
          "WGD QC: purity {:.2f} is below threshold {:.2f}, but genome-wide LOH fraction {:.2f} meets or exceeds "
          "{:.2f}; "
          "rejecting diploid-bias solution (WGD suspected), re-running with uniform prior (purity={:.2f}, "
          "ploidy={:.2f})",
          purity,
          kMinPurityToEvaluateWgdQc,
          loh_fraction,
          options.min_loh_fraction_for_wgd,
          purity,
          ploidy);
      return TotalCopyNumberPrior::kUniform;
    }
    Logging::Info(
        "WGD QC: purity {:.2f} is below threshold {:.2f} and genome-wide LOH fraction {:.2f} does not exceed {:.2f}, "
        "accepting diploid-bias solution",
        purity,
        kMinPurityToEvaluateWgdQc,
        loh_fraction,
        options.min_loh_fraction_for_wgd);
    return TotalCopyNumberPrior::kDiploidBias;
  }

  const auto candidate_segments = FindCandidateHetDelSegments(segments, logrs, ref_alt_depths, purity, ploidy, options);

  Logging::Info("WGD QC: {} candidate clonal het-del segments identified (need {})",
                candidate_segments.size(),
                kMinCandidateHetDelSegments);
  for (const auto& [seg, mbaf] : candidate_segments) {
    Logging::Info(
        "WGD QC:   candidate segment {}:{}-{} num_obs={} num_snps={} mean_logr={:.6f} mean_dh={:.4f} mbaf={:.4f}",
        seg->contig,
        seg->start,
        seg->end,
        seg->num_obs.value_or(0),
        seg->num_snps.value_or(0),
        seg->mean_logr.value(),
        seg->mean_dh.value_or(NAN),
        mbaf);
  }

  if (candidate_segments.size() < kMinCandidateHetDelSegments) {
    Logging::Info(
        "WGD QC: only {} candidate clonal het-del segments (need {}), accepting purity/ploidy solution "
        "(purity={:.2f}, ploidy={:.2f})",
        candidate_segments.size(),
        kMinCandidateHetDelSegments,
        purity,
        ploidy);
    return TotalCopyNumberPrior::kDiploidBias;
  }

  // Build mBAF vector from the already-computed values
  vec<f64> candidate_mbafs_vec;
  candidate_mbafs_vec.reserve(candidate_segments.size());
  for (const auto& [seg, mbaf] : candidate_segments) {
    candidate_mbafs_vec.push_back(mbaf);
  }
  const Eigen::VectorXd candidate_mbafs = eigen::FromStdVec(candidate_mbafs_vec);

  // Fit KDE on per-segment mBAF values (already in [0, 0.5])
  const Eigen::VectorXd weights = Eigen::VectorXd::Ones(candidate_mbafs.size());
  const f64 wgd_bandwidth_multiplier = 1.1;
  const auto kde = GaussianKernelDensityEstimate(candidate_mbafs, weights, wgd_bandwidth_multiplier);

  auto local_maxima = FindLocalMaxima(kde.y);
  if (local_maxima.size() < kNumTopPeaks) {
    Logging::Info(
        "WGD QC: fewer than 2 peaks in clonal heterozygous deletion mBAF distribution, accepting purity/ploidy "
        "solution (purity={:.2f}, ploidy={:.2f})",
        purity,
        ploidy);
    return TotalCopyNumberPrior::kDiploidBias;
  }

  // Find the top 2 peaks by density
  std::ranges::partial_sort(local_maxima, local_maxima.begin() + kNumTopPeaks, [&kde](const size_t a, const size_t b) {
    return kde.y(static_cast<Eigen::Index>(a)) > kde.y(static_cast<Eigen::Index>(b));
  });
  const f64 peak1 = kde.x(static_cast<Eigen::Index>(local_maxima[0]));
  const f64 peak2 = kde.x(static_cast<Eigen::Index>(local_maxima[1]));

  Logging::Info("WGD QC: top 2 clonal heterozygous deletion mBAF peaks at {:.3f} and {:.3f} (diff {:.3f})",
                peak1,
                peak2,
                std::abs(peak1 - peak2));

  if (std::abs(peak1 - peak2) > kWGDMaxPeakDifference) {
    Logging::Info(
        "WGD QC: peak difference exceeds {:.2f}, rejecting purity/ploidy solution (purity={:.2f}, ploidy={:.2f})",
        kWGDMaxPeakDifference,
        purity,
        ploidy);
    return TotalCopyNumberPrior::kUniform;
  }
  return TotalCopyNumberPrior::kDiploidBias;
}

/**
 * @brief Runs purity/ploidy search with an optional whole-genome-duplication (WGD) quality-control step.
 *
 * First performs a standard purity/ploidy search using a diploid-bias copy-number prior. If BAF data are available,
 * the result is then passed to EvaluatePurityPloidySolution. If that check determines the diploid-bias solution is
 * inconsistent with the observed mBAF distribution (i.e. WGD is suspected), the search is re-run with a uniform
 * copy-number prior and the updated purity/ploidy estimates are returned instead.
 *
 * When BAF data are absent the WGD QC step is skipped and the diploid-bias result is returned directly.
 *
 * @param segments Genomic segments to use in the search.
 * @param logrs Per-locus log-ratio observations.
 * @param ref_obvs Per-locus reference allele depths (optional; required for WGD QC).
 * @param alt_obvs Per-locus alternate allele depths (optional; required for WGD QC).
 * @param options Purity/ploidy search options (seg_max_logr, seg_min_num_logrs, seg_min_num_snps,
 *                baf_likelihood_weight). baf_likelihood_weight is applied to the BAF term only on the WGD re-run.
 * @param grid_out_fname Path to write the purity/ploidy likelihood grid (diploid-bias run).
 * @param sex Sex of the sample (used to handle allosome segments correctly).
 * @param command_line_info Metadata to prepend to output files.
 * @return PurityPloidySearchOut containing the selected purity, ploidy, likelihood grids, and the
 *         TotalCopyNumberPrior that was ultimately used (kDiploidBias or kUniform).
 */
PurityPloidySearchOut PurityPloidySearchWithWGDQC(const vec<GenomicSegment>& segments,
                                                  const Observations& logrs,
                                                  const std::optional<Observations>& ref_obvs,
                                                  const std::optional<Observations>& alt_obvs,
                                                  const PurityPloidySearchOptions& options,
                                                  const fs::path& grid_out_fname,
                                                  const std::string& output_prefix,
                                                  const sex_predict::Sex sex,
                                                  const io::CommandLineInfo& command_line_info) {
  // The per-segment likelihood diagnostic is opt-in (--write-segment-likelihood) and written one file per search pass,
  // alongside the grid output. A null writer pointer disables emission entirely.
  std::optional<SegmentLikelihoodWriter> diploid_bias_writer;
  if (options.write_segment_likelihood) {
    diploid_bias_writer.emplace(SegmentLikelihoodsParams{
        .sex = sex,
        .tumor_only = options.tumor_only,
        .prior = TotalCopyNumberPrior::kDiploidBias,
        .out_fname = file::PrefixedOutputPath(
            grid_out_fname.parent_path().parent_path(), output_prefix, kDefaultSegmentLikelihoodDiploidBiasOutput),
        .command_line_info = command_line_info});
  }
  const PurityPloidyGridParams diploid_bias_params{
      .seg_max_logr = options.seg_max_logr,
      .seg_min_num_logrs = options.seg_min_num_logrs,
      .seg_min_num_snps = options.seg_min_num_snps,
      .sex = sex,
      .use_diploid_bias_prior = true,
      .tumor_only = options.tumor_only,
      .likelihood_writer = diploid_bias_writer ? &diploid_bias_writer.value() : nullptr};
  PurityPloidySearchOut out =
      PurityPloidySearch(segments, logrs, ref_obvs, alt_obvs, diploid_bias_params, grid_out_fname, command_line_info);

  TotalCopyNumberPrior prior = TotalCopyNumberPrior::kDiploidBias;
  if (ref_obvs.has_value() && alt_obvs.has_value()) {
    // In tumor-only mode, classify segments at the top diploid-bias solution and drop germline-classified segments so
    // that only somatic segments contribute to the WGD QC mBAF analysis.
    const vec<GenomicSegment> wgd_qc_segments =
        options.tumor_only ? FilterSomaticSegmentsForWGDQC(segments, logrs, out.purity, out.ploidy, sex) : segments;
    prior = EvaluatePurityPloidySolution(wgd_qc_segments,
                                         logrs,
                                         RefAltObservations{.ref_obvs = ref_obvs.value(),
                                                            .alt_obvs = alt_obvs.value(),
                                                            .somatic_vafs = std::nullopt,
                                                            .vcf_check = std::nullopt},
                                         out.purity,
                                         out.ploidy,
                                         options,
                                         sex);
  }

  if (prior == TotalCopyNumberPrior::kUniform) {
    Logging::Info(
        "WGD QC: re-running purity/ploidy search with uniform copy number prior and BAF-likelihood weight {:.2f}",
        options.baf_likelihood_weight);
    std::optional<SegmentLikelihoodWriter> uniform_writer;
    if (options.write_segment_likelihood) {
      uniform_writer.emplace(SegmentLikelihoodsParams{
          .sex = sex,
          .tumor_only = options.tumor_only,
          .prior = TotalCopyNumberPrior::kUniform,
          .out_fname = file::PrefixedOutputPath(
              grid_out_fname.parent_path().parent_path(), output_prefix, kDefaultSegmentLikelihoodUniformOutput),
          .command_line_info = command_line_info});
    }
    // Re-run with uniform prior, overwriting grid_out_fname so the persisted grid matches the chosen solution
    const PurityPloidyGridParams uniform_params{
        .seg_max_logr = options.seg_max_logr,
        .seg_min_num_logrs = options.seg_min_num_logrs,
        .seg_min_num_snps = options.seg_min_num_snps,
        .sex = sex,
        .use_diploid_bias_prior = false,
        .tumor_only = options.tumor_only,
        // Trust BAF over log ratio on the re-run: the WGD QC step that brought us here rejected the initial solution
        // purely on BAF evidence, so up-weighting the BAF term lets that same signal drive the re-search.
        .baf_likelihood_weight = options.baf_likelihood_weight,
        .likelihood_writer = uniform_writer ? &uniform_writer.value() : nullptr};
    PurityPloidySearchOut uniform_out =
        PurityPloidySearch(segments, logrs, ref_obvs, alt_obvs, uniform_params, grid_out_fname, command_line_info);
    out.purity = uniform_out.purity;
    out.ploidy = uniform_out.ploidy;
    out.joint_ll_grid = std::move(uniform_out.joint_ll_grid);
    out.logr_ll_grid = std::move(uniform_out.logr_ll_grid);
    out.baf_ll_grid = std::move(uniform_out.baf_ll_grid);
    Logging::Info("WGD QC: uniform prior solution: purity={:.2f}, ploidy={:.2f}", out.purity, out.ploidy);
  }

  out.total_copy_number_prior = prior;
  return out;
}

}  // namespace xoos::cnc
