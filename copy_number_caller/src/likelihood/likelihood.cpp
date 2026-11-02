#include "likelihood/likelihood.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <fstream>

#include <fmt/core.h>

#include <xoos/enum/enum-util.h>
#include <xoos/error/error.h>
#include <xoos/io/metadata-util.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/log/logging.h>
#include <xoos/stats/copy-number-stats.h>
#include <xoos/stats/distributions.h>
#include <xoos/stats/stats.h>
#include <xoos/types/vec.h>
#include <xoos/util/file-functions.h>
#include <xoos/util/math.h>

#include "cli/copy-number-caller-options.h"
#include "cli/sample-metadata-options.h"
#include "io/copy-number-caller-default-filenames.h"
#include "io/write-segments.h"
#include "likelihood/baf-tools.h"
#include "likelihood/flag-copy-number-calls.h"
#include "seg-to-vcf/seg-to-vcf.h"
#include "segmentation/read-segments.h"
#include "segmentation/segments-header.h"
#include "util/eigen-math.h"
#include "util/mapq-util.h"

namespace xoos::cnc {
using std::log;
using xoos::file::PrefixedOutputPath;

// for CN likelihood estimation after purity-ploidy search
const size_t kLikelihoodMaxCN2 = 20;
constexpr size_t kNCN2 = kLikelihoodMaxCN2 + 1;

// simple helper function to resolve sample ID for VCF output from the input segments
namespace {
std::string ResolveSampleIdFromSegments(const vec<GenomicSegment>& segments, const std::string& fallback_sample_id) {
  if (segments.empty() || fallback_sample_id != kSampleMetadataOptionsDefaultSampleId) {
    return fallback_sample_id;
  }
  for (const auto& seg : segments) {
    if (seg.id.has_value() && !seg.id->empty()) {
      return seg.id.value();
    }
  }
  return fallback_sample_id;
}
}  // namespace

/**
 * @brief Use global tumor purity/ploidy estimates (user-provided or optionally estimate those here) to calculate the
 * joint logr-BAF likelihood of each segment. In calculating this joint likelihood we also return the most likely total
 * copy number of the segment as well as the the most likely minor copy number. See docs/likelihood.md for more
 * information
 * @param params
 */

void CalculateLikelihoodsMain(const CopyNumberCallerOptions& options) {
  using enum_util::FormatEnumName;
  const auto& dir = options.output_dir;
  const auto& prefix = options.output_prefix;
  const auto likelihood_out = (options.likelihood_options.mode == LikelihoodMode::kSomatic)
                                  ? PrefixedOutputPath(dir, prefix, kDefaultSomaticCNCallsetSegOutput)
                                  : PrefixedOutputPath(dir, prefix, kDefaultGermlineCNCallsetSegOutput);
  const auto vcf_out = PrefixedOutputPath(dir, prefix, kDefaultSomaticCNCallsetVcfOutput);
  vec<fs::path> output_paths = {likelihood_out};
  if (options.likelihood_options.mode == LikelihoodMode::kSomatic) {
    output_paths.push_back(vcf_out);
  }
  file::CheckOutputFilePathsForWriting(output_paths, options.overwrite);
  if (!options.segments_fname.has_value()) {
    throw error::Error("CalculateLikelihoods requires a segments file (--segments)");
  }
  if (!options.logrs_fname.has_value()) {
    throw error::Error("CalculateLikelihoods requires a logRs file (--logrs)");
  }
  vec<GenomicSegment> segments =
      ReadSegments(options.segments_fname.value(), options.likelihood_options.input_segments_type);
  const std::string sample_id = ResolveSampleIdFromSegments(segments, options.sample_metadata_options.sample_id);
  const SegmentsHeader header = ReadHeaderFromSegments(options.segments_fname.value());
  std::ifstream logr_ifstream(options.logrs_fname.value());
  const Observations logrs = ReadObservations(logr_ifstream);
  LikelihoodModel ll_model = options.likelihood_options.likelihood_model;
  if ((ll_model == LikelihoodModel::kSerialSummarized) && !options.bafs_fname.has_value()) {
    Logging::Debug("likelihood model {} specified, but no BAFs provided. Defaulting to {}",
                   FormatEnumName(ll_model),
                   FormatEnumName(LikelihoodModel::kLogrSummarizedOnly));
    ll_model = LikelihoodModel::kLogrSummarizedOnly;
  }
  std::optional<Observations> ref_ads;
  std::optional<Observations> alt_ads;
  if (options.bafs_fname.has_value()) {
    auto ref_alts = GetDepthsFromFile(options.bafs_fname.value());
    ref_ads = std::move(ref_alts.ref_obvs);
    alt_ads = std::move(ref_alts.alt_obvs);
  }
  std::optional<Observations> mapqs;
  if (options.mapping_qualities_fname.has_value()) {
    std::ifstream mapq_ifs(options.mapping_qualities_fname.value());
    mapqs = ReadObservations(mapq_ifs);
  }
  const auto ll_opts = options.likelihood_options;
  segments = CalculateLikelihoods(
      segments, logrs, ref_ads, alt_ads, mapqs, ll_model, ll_opts, options.sample_metadata_options);
  for (auto& seg : segments) {
    seg.sex = options.sample_metadata_options.sex;
    seg.id = sample_id;
  }
  // In tumor-only mode, drop the same allosomal segments the purity/ploidy search excluded so that CNV prediction and
  // output cover exactly the segments used to fit purity and ploidy. This mirrors the somatic ploidy policy: chrX
  // pseudo-autosomal regions (any sex) and female chrX are retained; chrY and male/unknown non-PAR chrX are dropped.
  if (options.likelihood_options.tumor_only) {
    const sex_predict::Sex sex = options.sample_metadata_options.sex.value();
    std::erase_if(segments, [sex](const GenomicSegment& seg) {
      return SkipBasedOnModeAndSex(seg, LikelihoodMode::kSomatic, sex);
    });
  }
  WriteSegments(likelihood_out, segments, options.likelihood_options.output_segments_type, options.command_line_info);
  if (options.likelihood_options.mode == LikelihoodMode::kSomatic) {
    const std::string reference_file =
        header.reference_file.empty() ? options.reference_genome_fname.string() : header.reference_file;
    const SeqLenMap seq_lengths = ParseFai(options.reference_genome_fai_fname);
    const StrUnorderedMap<size_t> seq_order = GetContigOrder(options.reference_genome_fai_fname);
    WriteSegmentsToVcf(vcf_out,
                       segments,
                       {.seq_lengths = seq_lengths,
                        .seq_order = seq_order,
                        .sample_id = sample_id,
                        .reference_file = reference_file,
                        .sex = options.sample_metadata_options.sex.value(),
                        .mode = options.likelihood_options.tumor_only ? CopyNumberCallerModes::kSomaticTumorOnlyWGS
                                                                      : CopyNumberCallerModes::kSomaticTumorNormalWGS,
                        .purity = options.sample_metadata_options.purity,
                        .ploidy = options.sample_metadata_options.ploidy,
                        .vcf_purity_source = VcfPuritySource::kUserInput,
                        .command_line_info = options.command_line_info});
  }
}

vec<GenomicSegment>& CalculateLikelihoods(vec<GenomicSegment>& segments,
                                          const Observations& logrs,
                                          const std::optional<Observations>& ref_ads,
                                          const std::optional<Observations>& alt_ads,
                                          const std::optional<Observations>& mapqs,
                                          LikelihoodModel ll_model,
                                          const LikelihoodOptions& options,
                                          const SampleMetadataOptions& sample_metadata_options) {
  // calculate likelihoods
  if (!sample_metadata_options.sex.has_value()) {
    throw error::Error("Must provide sex to CalculateLikelihoods");
  }
  const sex_predict::Sex sex = sample_metadata_options.sex.value();
  if (!sample_metadata_options.purity.has_value() || !sample_metadata_options.ploidy.has_value()) {
    Logging::Error("Must provide purity and ploidy to CalculateLikelihoods!");
    throw std::runtime_error("purity and ploidy not provided to CalculateLikelihoods");
  }
  const f64 purity = sample_metadata_options.purity.value();
  const f64 ploidy = sample_metadata_options.ploidy.value();
  switch (ll_model) {
    case LikelihoodModel::kSerialSummarized: {
      if (!ref_ads.has_value() || !alt_ads.has_value()) {
        throw std::runtime_error("Serial model requires both ref and alt allele depths");
      }
      segments = AllSegmentsSerialLogRBAFLikelihoods(
          segments, logrs, ref_ads.value(), alt_ads.value(), purity, ploidy, options.mode, sex);
      break;
    }
    case LikelihoodModel::kLogrSummarizedOnly: {
      segments = AllSegmentsLogRLikelihoods(segments, logrs, purity, ploidy, options.mode, sex);
      break;
    }
    default: {
      Logging::Error("Invalid likelihood model specified for CalculateLikelihoods");
      throw std::runtime_error("no likelihood model specified (this should never happen)");
    }
  }
  // add metadata to the segments
  // also initiate an empty flags std::vector
  for (auto& seg : segments) {
    seg.sex = sample_metadata_options.sex;
    seg.id = sample_metadata_options.sample_id;
    seg.flags.emplace(vec<std::string>{});
  }
  // flag calls here
  // assign mean MAPQs and flag if necessary
  if (mapqs.has_value()) {
    const Eigen::VectorXd avg_mean_mapqs = GetAvgMeanMapqPerSegment(segments, mapqs.value());
    for (size_t i = 0; i < segments.size(); ++i) {
      segments[i].avg_mean_mapq = avg_mean_mapqs[static_cast<Eigen::Index>(i)];
    }
    if (options.mapq_cutoff_for_calls > 0) {
      FlagCallsByMeanMapq(segments, options.mapq_cutoff_for_calls);
    }
  }
  // flag by expected LogR
  FlagCallsByExpectedTotalCopyNumber(segments, sex, options.non_integer_tcn_filter_options);
  FlagCallsByLength(segments, options.cnv_length_flag_min_size);
  FlagChrYCallsIfFemale(segments);
  // Tumor-only mode: classify each segment as germline or somatic. This runs last so the recomputed expected total copy
  // number on germline-classified segments is not overwritten by FlagCallsByExpectedTotalCopyNumber.
  if (options.tumor_only) {
    ClassifyGermlineVsSomaticSegments(segments, logrs, sex);
  }
  return segments;
}

bool SkipBasedOnModeAndSex(const GenomicSegment& seg, const LikelihoodMode mode, const sex_predict::Sex sex) {
  // skip if:
  // 1. sex is unknown and in allosome
  // 2. sex is female and in chrY (keep all of chrX regardless of mode)
  // 3. sex is male and mode is somatic, and we are in chrY or a non-PAR region
  // 4. sex is male and mode is germline and in PAR region
  using enum sex_predict::Sex;
  vec<bool> conditions{sex == kUnknown && seg.in_allosome.value_or(false),
                       sex == kFemale && seg.in_allosome.value_or(false) && sex_predict::IsInChromY(seg.contig),
                       sex == kMale && mode == LikelihoodMode::kSomatic && seg.in_allosome.value_or(false) &&
                           (!seg.in_pseudo_autosomal_region.value_or(false) || sex_predict::IsInChromY(seg.contig)),
                       sex == kMale && mode == LikelihoodMode::kGermline && seg.in_allosome.value_or(false) &&
                           sex_predict::IsInChromY(seg.contig) && seg.in_pseudo_autosomal_region.value_or(false)};
  return std::ranges::any_of(conditions, [](bool x) { return x; });
}

bool ExpectHaploid(const GenomicSegment& seg, const LikelihoodMode mode, const sex_predict::Sex sex) {
  // Only expect haploid in germline mode, in non-PAR region of an allosome, and if sex is male
  return (mode == LikelihoodMode::kGermline && seg.in_allosome.value_or(false) &&
          !seg.in_pseudo_autosomal_region.value_or(false) && sex == sex_predict::Sex::kMale);
}

namespace {
// A segment can be compared against the germline model only when the somatic model can call it and a copy-number
// solution with the required log-ratio inputs is present.
bool SegmentClassifiableForGermlineComparison(const GenomicSegment& seg, const sex_predict::Sex sex) {
  if (SkipBasedOnModeAndSex(seg, LikelihoodMode::kSomatic, sex)) {
    return false;
  }
  return seg.logr_likelihood.has_value() && seg.purity.has_value() && seg.ploidy.has_value() &&
         seg.mean_logr.has_value();
}
}  // namespace

/**
 * @brief Estimate total copy number per segment using only logr obserbations
 * @param segments
 * @param logrs
 * @param purity
 * @param ploidy
 * @return  std::vector of GenomicSegment
 */
vec<GenomicSegment>& AllSegmentsLogRLikelihoods(vec<GenomicSegment>& segments,
                                                const Observations& logrs,
                                                f64 purity,
                                                f64 ploidy,
                                                LikelihoodMode mode,
                                                const sex_predict::Sex sex) {
  const IntervalTrees logr_trees(logrs);
  // we need to conisder all multiplicities from 0 to kMaxLikelihood2 inclusive
  // idx=coverage
  Eigen::VectorXd logr_likelihoods(kNCN2);
  const f64 mean_segment_logr_sd = GetMeanSegmentLogRSD(segments, logrs, logr_trees);
  for (auto& seg : segments) {
    seg.ploidy = ploidy;
    seg.purity = purity;
    logr_likelihoods.setConstant(kNegativeInf);
    if (SkipBasedOnModeAndSex(seg, mode, sex)) {
      seg.total_copy_number = 0;
      seg.logr_likelihood = std::nullopt;
      seg.avg_mean_mapq = std::nullopt;
      seg.expected_total_copy_number = 0;
      seg.joint_likelihood = std::nullopt;
      continue;
    }
    const bool expect_haploid = ExpectHaploid(seg, mode, sex);
    // If there are no observations and mean logR has no value, then this should represent a seed segment with no
    // overlapping observations. We can skip likelihood model predictions on this segment.
    if ((!seg.num_obs.has_value()) || (seg.num_obs.value() == 0 && !seg.mean_logr.has_value())) {
      continue;
    }
    const vec<size_t> logr_idxs = logr_trees.LookUp(seg.contig, seg.start, seg.end);
    const Observations logr_obvs = logrs.FilterByIdxs(logr_idxs);
    SegmentAllSummarizedLogRLikelihoods(
        logr_obvs.obvs, mean_segment_logr_sd, purity, ploidy, expect_haploid, logr_likelihoods);
    // Choose best TCN state based on LogR likelihood alone
    size_t total_cn = eigen::IndexMax(logr_likelihoods);
    seg.total_copy_number = total_cn;
    seg.logr_likelihood = logr_likelihoods(static_cast<Eigen::Index>(total_cn));
    seg.joint_likelihood = logr_likelihoods(static_cast<Eigen::Index>(total_cn));
  }
  return segments;
}

/**
 * @brief Serial method for returning joint LogR-BAF likelihoods for all
 * segments. Here, we choose the total copy number state (TCN) solely by the
 * best LogR likelihood (as opposed to the joint method, which chooses the TCN
 * by the joint LogR-BAF likelihood).  The maximum BAF likelihood at this TCN
 * state is then added to the LogR likelihood at this state. The returned value
 * is this joint likelihood along with the corresponding TCN state and minor CN
 * state.
 * @param segments
 * @param logrs
 * @param ref_ads
 * @param alt_ads
 * @param purity
 * @param ploidy
 * @param mode (kGermline or kSomatic)
 * @param sex
 * @return  std::vector of GenomicSegment
 */
vec<GenomicSegment>& AllSegmentsSerialLogRBAFLikelihoods(vec<GenomicSegment>& segments,
                                                         const Observations& logrs,
                                                         const Observations& ref_ads,
                                                         const Observations& alt_ads,
                                                         f64 purity,
                                                         f64 ploidy,
                                                         LikelihoodMode mode,
                                                         const sex_predict::Sex sex) {
  const IntervalTrees logr_trees(logrs);
  IntervalTrees variant_trees(ref_ads);
  f64 mean_segment_baf_sd = 0;
  mean_segment_baf_sd = GetMeanSegmentMBAFSD(segments, ref_ads, alt_ads, variant_trees);
  // we need to conisder all multiplicities from 0 to kMaxLikelihood2 inclusive
  // idx=coverage
  Eigen::VectorXd logr_likelihoods(kNCN2);
  // row_idx=TCN, col_idx=multiplicity
  Eigen::VectorXd baf_likelihoods(kNCN2);
  const f64 mean_segment_logr_sd = GetMeanSegmentLogRSD(segments, logrs, logr_trees);
  for (auto& seg : segments) {
    seg.ploidy = ploidy;
    seg.purity = purity;
    if (SkipBasedOnModeAndSex(seg, mode, sex)) {
      seg.total_copy_number = 0;
      seg.minor_copy_number = 0;
      seg.major_copy_number = 0;
      seg.logr_likelihood = std::nullopt;
      seg.baf_likelihood = std::nullopt;
      seg.joint_likelihood = std::nullopt;
      seg.avg_mean_mapq = std::nullopt;
      seg.expected_total_copy_number = 0;
      continue;
    }
    const bool expect_haploid = ExpectHaploid(seg, mode, sex);
    // If there are no observations and mean logR has no value, then this should represent a seed segment with no
    // overlapping observations. We can skip likelihood model predictions on this segment.
    if ((!seg.num_obs.has_value()) || (seg.num_obs.value() == 0 && !seg.mean_logr.has_value())) {
      continue;
    }
    // gather data from the segments
    const vec<size_t> logr_idxs = logr_trees.LookUp(seg.contig, seg.start, seg.end);
    const Observations logr_obvs = logrs.FilterByIdxs(logr_idxs);
    const vec<size_t> variant_idxs = variant_trees.LookUp(seg.contig, seg.start, seg.end);
    const Observations ref_ad_subset = ref_ads.FilterByIdxs(variant_idxs);
    const Observations alt_ad_subset = alt_ads.FilterByIdxs(variant_idxs);
    if (logr_idxs.empty()) {
      Logging::Warn("Segment {}:{}-{} has 0 LogR observations. Probably a seed segment. Skipping...",
                    seg.contig,
                    seg.start + 1,
                    seg.end);
      continue;
    }
    if (variant_idxs.empty()) {
      Logging::Warn(
          "Segment {}:{}-{} has 0 BAF observations. Will set this segment's joint likelihood to its LogR likelihood, "
          "and set MajorCopyNumber and MinorCopyNumber to NA",
          seg.contig,
          seg.start + 1,
          seg.end);
      seg.num_snps = 0;
      seg.major_copy_number = std::nullopt;
      seg.minor_copy_number = std::nullopt;
      seg.mbaf = std::nullopt;
    } else {
      f64 mbaf = GetPeakOfBAFDistributionThenMirror(ref_ad_subset.obvs, alt_ad_subset.obvs);
      seg.mbaf = std::isnan(mbaf) ? std::nullopt : std::optional<f64>(mbaf);
    }
    // then, find all copy number likelihoods
    seg = SegmentSerialSummarizedLogRBAFLikelihood(seg,
                                                   logr_obvs.obvs,
                                                   purity,
                                                   ploidy,
                                                   mean_segment_logr_sd,
                                                   mean_segment_baf_sd,
                                                   false,
                                                   expect_haploid,
                                                   logr_likelihoods,
                                                   baf_likelihoods,
                                                   /*baf_likelihood_weight=*/1.0);
  }
  return segments;
}

vec<GenomicSegment>& ClassifyGermlineVsSomaticSegments(vec<GenomicSegment>& segments,
                                                       const Observations& logrs,
                                                       const sex_predict::Sex sex) {
  const IntervalTrees logr_trees(logrs);
  const f64 mean_segment_logr_sd = GetMeanSegmentLogRSD(segments, logrs, logr_trees);
  // Scratch buffers reused across segments. Both germline and somatic comparisons run over the same eight
  // copy-number states (0-7) with the diploid bias prior so the log-ratio likelihoods are directly comparable.
  Eigen::VectorXd germline_logr_likelihoods(static_cast<Eigen::Index>(kLikelihoodMaxCN1 + 1));
  Eigen::VectorXd somatic_logr_likelihoods(static_cast<Eigen::Index>(kLikelihoodMaxCN1 + 1));
  Eigen::VectorXd somatic_baf_likelihoods(static_cast<Eigen::Index>(kLikelihoodMaxCN1 + 1));
  for (auto& seg : segments) {
    // Allosome handling mirrors the purity/ploidy search grid: the somatic ploidy policy governs which allosomal
    // segments are usable. chrX (female) and chrX pseudo-autosomal regions are retained; chrY and male/unknown non-PAR
    // chrX are skipped. Segments the somatic model cannot call (or that lack a copy-number solution) are left
    // unclassified.
    if (!SegmentClassifiableForGermlineComparison(seg, sex)) {
      continue;
    }
    const vec<size_t> logr_idxs = logr_trees.LookUp(seg.contig, seg.start, seg.end);
    const Observations logr_subset = logrs.FilterByIdxs(logr_idxs);
    if (logr_subset.obvs.size() == 0) {
      continue;
    }
    const bool expect_haploid = ExpectHaploid(seg, LikelihoodMode::kSomatic, sex);
    const GermlineLikelihoodResult germline =
        ComputeGermlineLikelihood(logr_subset.obvs, mean_segment_logr_sd, expect_haploid, germline_logr_likelihoods);
    // Re-evaluate the somatic log-ratio likelihood under the diploid bias prior over the same eight states so the
    // comparison matches the germline model. This mirrors the classification used by the purity-ploidy grid search.
    GenomicSegment somatic_probe;
    somatic_probe = SegmentSerialSummarizedLogRBAFLikelihood(somatic_probe,
                                                             logr_subset.obvs,
                                                             seg.purity.value(),
                                                             seg.ploidy.value(),
                                                             mean_segment_logr_sd,
                                                             /*mean_baf_sd=*/0.0,
                                                             /*prior_2n=*/true,
                                                             expect_haploid,
                                                             somatic_logr_likelihoods,
                                                             somatic_baf_likelihoods,
                                                             /*baf_likelihood_weight=*/1.0);
    const f64 somatic_logr_likelihood = somatic_probe.logr_likelihood.value();
    // Break ties toward germline using the same margin as the purity/ploidy grid search.
    if (germline.likelihood >= (somatic_logr_likelihood - kGermlineTieMargin)) {
      // Germline CNV: override the total copy number with the germline best fit. The germline model is log-ratio-only,
      // so there is no allele-specific resolution and major/minor CN are cleared. The expected total copy number is
      // recomputed at the germline purity/ploidy while the reported Purity/Ploidy columns keep the somatic solution.
      seg.is_germline_cnv = true;
      seg.total_copy_number = germline.total_copy_number;
      seg.major_copy_number = std::nullopt;
      seg.minor_copy_number = std::nullopt;
      seg.expected_total_copy_number =
          stats::ExpectedTotalCopyNumber(kGermlinePurity, kDiploidPloidy, seg.mean_logr.value(), expect_haploid);
    } else {
      // Somatic: keep the allele-specific state computed during prediction.
      seg.is_germline_cnv = false;
    }
  }
  return segments;
}

const f64 kLikelihoodDefaultDF = 5;

f64 SegmentSummarizedLogRLikelihood(const Eigen::VectorXd& logrs,  // NOSONAR(cpp:S5419) adjacent f64 model params
                                    const f64 sd,
                                    const f64 purity,
                                    const f64 tumor_ploidy,
                                    const f64 total_copy_number,
                                    const bool expect_haploid) {
  const f64 expected_logr = stats::ExpectedLogR(purity, tumor_ploidy, total_copy_number, expect_haploid);
  const f64 mean_logr = eigen::Mean(logrs);
  f64 pdf = stats::NonStandardTPDF(mean_logr, kLikelihoodDefaultDF, expected_logr, sd);
  pdf = !math::IsCloseToZero(pdf) ? pdf : DBL_MIN;
  return log(pdf);
}

/**
 * @brief Calculates the BAF likelihood of a segment
 * @param mbaf - segment mirrored BAF estimate
 * @param sd - standard deviation of mirrored BAFs across all segments
 * @param purity - tumor purity
 * @param total_copy_number
 * @param multiplicity - minor copy number
 * @param expect_haploid
 * @return
 */
namespace {
f64 SegmentSummarizedBAFLikelihood(const f64 mbaf,
                                   const f64 sd,
                                   const f64 purity,
                                   const f64 total_copy_number,
                                   const f64 multiplicity,
                                   const bool expect_haploid) {
  if (expect_haploid) {
    return NAN;
  }
  const f64 expected_af = stats::ExpectedAF(purity, total_copy_number, multiplicity, expect_haploid);
  return log(stats::NonStandardTPDF(mbaf, kLikelihoodDefaultDF, expected_af, sd));
}
}  // namespace

f64 GetMeanSegmentLogRSD(const vec<GenomicSegment>& segments,
                         const Observations& logrs,
                         const IntervalTrees& logr_trees) {
  f64 total_sd = 0;
  f64 denom = 0;
  for (const auto& seg : segments) {
    if (sex_predict::IsInAllosome(seg.contig)) {
      continue;
    }
    const vec<size_t> logr_idxs = logr_trees.LookUp(seg.contig, seg.start, seg.end);
    const Observations logr_obvs = logrs.FilterByIdxs(logr_idxs);
    const f64 sd = eigen::Stddev(logr_obvs.obvs);
    if (std::isnan(sd)) {
      throw error::Error("Segment {}:{}-{} has NaN logR standard deviation", seg.contig, seg.start, seg.end);
    }
    total_sd += sd;
    denom += 1;
  }

  if (math::IsCloseToZero(denom)) {
    throw error::Error("No segments on autosomal chromosomes; cannot calculate mean segment logR standard deviation");
  }

  return total_sd / denom;
}

/**
 * @brief get the standard deviation of the weighted mBAF means of a set of segments
 * @param segments
 * @param ref_ads
 * @param alt_ads
 * @param variant_trees
 * @return
 */
f64 GetMeanSegmentMBAFSD(const vec<GenomicSegment>& segments,
                         const Observations& ref_ads,
                         const Observations& alt_ads,
                         IntervalTrees& variant_trees) {
  f64 num = 0;
  f64 denom = 0;
  for (const auto& seg : segments) {
    if (sex_predict::IsInAllosome(seg.contig)) {
      continue;
    }
    const vec<size_t> variant_idxs = variant_trees.LookUp(seg.contig, seg.start, seg.end);
    if (variant_idxs.size() > 1) {
      const Observations ref_ad_subset = ref_ads.FilterByIdxs(variant_idxs);
      const Observations alt_ad_subset = alt_ads.FilterByIdxs(variant_idxs);
      const f64 mbaf_mean = GetPeakOfBAFDistributionThenMirror(ref_ad_subset.obvs, alt_ad_subset.obvs);
      num += WeightedMBAFSD(ref_ad_subset.obvs, alt_ad_subset.obvs, mbaf_mean);
      denom += 1;
    }
  }

  if (math::IsCloseToZero(denom)) {
    throw error::Error("No autosomal segments with >1 variant; cannot calculate mean segment mBAF standard deviation");
  }

  return num / denom;
}

/**
 * @brief Get the summarized likelihoods (mean of point for all total CNs from 0 -> ret.size()
 * likelihood for TCN=c: likelihood(mean(logrs) | c)
 * @param logrs - std::vector of observations associated with the segment
 * @param sd - standard deviation of mean logr across all segments
 * @param purity - tumor purity
 * @param tumor_ploidy
 * @param ret - std::vector to fill. size of std::vector determines how many likelihoods to calculate [0, ret.size()-1]
 * @return reference to ret
 */
Eigen::VectorXd& SegmentAllSummarizedLogRLikelihoods(
    const Eigen::VectorXd& logrs,  // NOSONAR(cpp:S5419) adjacent f64 model params
    const f64 sd,
    const f64 purity,
    const f64 tumor_ploidy,
    const bool expect_haploid,
    Eigen::VectorXd& ret) {
  for (Eigen::Index c = 0; c < ret.size(); ++c) {
    ret(c) = SegmentSummarizedLogRLikelihood(logrs, sd, purity, tumor_ploidy, static_cast<f64>(c), expect_haploid);
  }
  return ret;
}

/**
 * @brief Get all the summarized BAF likelihoods associated with a segment. The number of likelihoods is each possible
 * minor copy number, ie. total_copy_number / 2 likelihood for TCN=c, MinorCN=m: likelihood(mean(baf) | c,m) for m in
 * [min(a,d) for a,d in zip(ref_ads, alt_ads)]
 * @param ref_ads - reference AD observations associated with the segment
 * @param alt_ads - alt AD observations associated with the segment
 * @param sd - standard deviation of mean mBAF across all segments
 * @param purity  - tumor purity
 * @param total_copy_number
 * @param ret - std::vector to fill. size of std::vector must be > total_copy_number/2
 * @return reference to ret
 */
SegmentAllSummarizedBAFLikelihoodsResults SegmentAllSummarizedBAFLikelihoods(
    const f64 mbaf, const f64 sd, const f64 purity, const f64 total_copy_number, const bool expect_haploid) {
  // calculate baf likelihood for each TCN/MinorCN combo
  auto tcn_int = static_cast<size_t>(total_copy_number);
  // row_idx=TCN, col_idx=multiplicity
  Eigen::VectorXd baf_likelihoods(kNCN2);
  baf_likelihoods.setConstant(kNegativeInf);
  for (size_t m = 0; m <= tcn_int / 2; ++m) {
    baf_likelihoods(static_cast<Eigen::Index>(m)) =
        SegmentSummarizedBAFLikelihood(mbaf, sd, purity, total_copy_number, static_cast<f64>(m), expect_haploid);
  }
  return {baf_likelihoods};
}

namespace {
f64 GetTcnPrior(const size_t n_cns, const size_t cn) {
  // Consider the copy neutral state (2) as 2 possible states. This makes the total number of states n_cns + 1.
  // This allows the copy neutral state to have f64 the prior.
  const auto n = static_cast<f64>(n_cns) + 1;
  constexpr f64 kCopyNeutralPriorWeight = 2.0;
  return (cn == kDiploidCopyNumber) ? log(kCopyNeutralPriorWeight / n) : log(1.0 / n);
}
}  // namespace

GermlineLikelihoodResult ComputeGermlineLikelihood(const Eigen::VectorXd& logrs,
                                                   const f64 mean_logr_sd,
                                                   const bool expect_haploid,
                                                   Eigen::VectorXd& logr_likelihoods) {
  // Germline events are assumed to sit in a normal diploid background, so the germline log-ratio
  // likelihood is evaluated at fixed purity = 0.99 and ploidy = 2.
  logr_likelihoods.setConstant(kNegativeInf);
  SegmentAllSummarizedLogRLikelihoods(
      logrs, mean_logr_sd, kGermlinePurity, kDiploidPloidy, expect_haploid, logr_likelihoods);
  // Apply the diploid bias prior across all copy-number states, exactly as in the somatic model, so the
  // non-diploid states are weighted identically. The copy-neutral state (CN = 2) is then forced to -Inf so a
  // germline segment can never be classified as copy-neutral.
  for (size_t i = 0; i < static_cast<size_t>(logr_likelihoods.size()); ++i) {
    logr_likelihoods[static_cast<Eigen::Index>(i)] += GetTcnPrior(static_cast<size_t>(logr_likelihoods.size()), i);
  }
  if (static_cast<size_t>(logr_likelihoods.size()) > kCopyNeutralState) {
    logr_likelihoods[static_cast<Eigen::Index>(kCopyNeutralState)] = kNegativeInf;
  }
  const u64 best_cn = eigen::IndexMax(logr_likelihoods);
  return {.likelihood = logr_likelihoods[static_cast<Eigen::Index>(best_cn)], .total_copy_number = best_cn};
}

f64 SegmentGermlineSummarizedLogRLikelihood(const Eigen::VectorXd& logrs,
                                            const f64 mean_logr_sd,
                                            const bool expect_haploid,
                                            Eigen::VectorXd& logr_likelihoods) {
  return ComputeGermlineLikelihood(logrs, mean_logr_sd, expect_haploid, logr_likelihoods).likelihood;
}

GenomicSegment& SegmentSerialSummarizedLogRBAFLikelihood(GenomicSegment& segment,
                                                         const Eigen::VectorXd& logrs,
                                                         f64 purity,
                                                         f64 ploidy,
                                                         const f64 mean_logr_sd,
                                                         const f64 mean_baf_sd,
                                                         const bool prior_2n,
                                                         const bool expect_haploid,
                                                         Eigen::VectorXd& logr_likelihoods,
                                                         Eigen::VectorXd& baf_likelihoods,
                                                         const f64 baf_likelihood_weight) {
  logr_likelihoods.setConstant(kNegativeInf);
  SegmentAllSummarizedLogRLikelihoods(logrs, mean_logr_sd, purity, ploidy, expect_haploid, logr_likelihoods);
  // Choose best TCN state based on LogR likelihood alone
  if (prior_2n) {
    for (size_t i = 0; i < static_cast<size_t>(logr_likelihoods.size()); ++i) {
      logr_likelihoods[static_cast<Eigen::Index>(i)] += GetTcnPrior(static_cast<size_t>(logr_likelihoods.size()), i);
    }
  }
  size_t total_cn = eigen::IndexMax(logr_likelihoods);
  f64 logr_likelihood = logr_likelihoods(static_cast<Eigen::Index>(total_cn));
  // calculate baf likelihoods for the total CN determined from the logr_likelihoods, over all minor CNs
  if (segment.mbaf.has_value() && !expect_haploid) {
    const auto baf_results = SegmentAllSummarizedBAFLikelihoods(
        segment.mbaf.value(), mean_baf_sd, purity, static_cast<f64>(total_cn), expect_haploid);
    baf_likelihoods = baf_results.baf_likelihoods;
    size_t minor_cn = eigen::IndexMax(baf_likelihoods);
    segment.total_copy_number = total_cn;
    segment.minor_copy_number = minor_cn;
    segment.major_copy_number = total_cn - minor_cn;
    segment.logr_likelihood = logr_likelihood;
    segment.baf_likelihood = baf_likelihoods(static_cast<Eigen::Index>(minor_cn));
    segment.joint_likelihood =
        segment.logr_likelihood.value() + (baf_likelihood_weight * segment.baf_likelihood.value());

  } else {
    segment.total_copy_number = total_cn;
    segment.logr_likelihood = logr_likelihood;
    segment.baf_likelihood = std::nullopt;
    segment.joint_likelihood = logr_likelihood;
  }
  segment.ploidy = ploidy;
  segment.purity = purity;
  return segment;
}

}  // namespace xoos::cnc
