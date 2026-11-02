#pragma once
#include <Eigen/Dense>
#include <limits>

#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "cli/copy-number-caller-options.h"
#include "likelihood/likelihood-mode.h"
#include "likelihood/likelihood-model.h"
#include "likelihood/likelihood-options.h"
#include "observations/observations.h"
#include "segmentation/genomic-segments.h"
#include "segmentation/interval-trees.h"

namespace xoos::cnc {

const f64 kNegativeInf = -std::numeric_limits<f64>::infinity();

// Maximum total copy-number state (minus one) used for CN likelihood estimation in the purity-ploidy search.
const size_t kLikelihoodMaxCN1 = 7;
// The mirror process halves the BAF variance, so the mean segment BAF SD is scaled back up by this factor.
const f64 kMirrorMBAFSDScale = 2.0;

// Fixed purity used to evaluate the germline log-ratio likelihood in tumor-only mode. Germline events are
// assumed to sit in a normal diploid background, hence purity ~ 1 and ploidy = kDiploidPloidy.
const f64 kGermlinePurity = 0.99;
// Copy-neutral (diploid) total copy-number state, excluded from germline classification.
const size_t kCopyNeutralState = 2;
// Tolerance for the germline-vs-somatic log-ratio likelihood comparison. A segment is classified germline when its
// germline log-ratio likelihood is within this margin of (or above) the somatic log-ratio likelihood. The somatic and
// germline expected log-ratios converge in the high-purity, near-diploid corner, producing near-ties (down to a single
// ULP) where a high-copy-number somatic state at low purity mimics a germline event. A bit-exact comparison cannot
// resolve these, so ties within this margin are broken toward germline to avoid depositing the somatic joint
// (log-ratio + BAF) into high-purity cells and inflating the purity estimate.
const f64 kGermlineTieMargin = 1e-3;

struct SegmentAllSummarizedBAFLikelihoodsResults {
  Eigen::VectorXd baf_likelihoods;
};

void CalculateLikelihoodsMain(const CopyNumberCallerOptions& options);
vec<GenomicSegment>& CalculateLikelihoods(vec<GenomicSegment>& segments,
                                          const Observations& logrs,
                                          const std::optional<Observations>& ref_ads,
                                          const std::optional<Observations>& alt_ads,
                                          const std::optional<Observations>& mapqs,
                                          LikelihoodModel ll_model,
                                          const LikelihoodOptions& options,
                                          const SampleMetadataOptions& sample_metadata_options);
vec<GenomicSegment>& AllSegmentsLogRLikelihoods(vec<GenomicSegment>& segments,
                                                const Observations& logrs,
                                                f64 purity,
                                                f64 ploidy,
                                                LikelihoodMode mode,
                                                sex_predict::Sex sex);
vec<GenomicSegment>& AllSegmentsSerialLogRBAFLikelihoods(vec<GenomicSegment>& segments,
                                                         const Observations& logrs,
                                                         const Observations& ref_ads,
                                                         const Observations& alt_ads,
                                                         f64 purity,
                                                         f64 ploidy,
                                                         LikelihoodMode mode,
                                                         sex_predict::Sex sex);
// Classify each predicted segment as a germline or somatic copy-number event by comparing the germline (fixed diploid
// purity/ploidy) and somatic log-ratio likelihoods. Runs as tumor-only post-processing. Germline-classified segments
// have their total and expected copy number recomputed under the germline model and their allele-specific copy numbers
// cleared; the segment's is_germline_cnv field records the outcome. Allosome handling mirrors the purity/ploidy search
// grid.
vec<GenomicSegment>& ClassifyGermlineVsSomaticSegments(vec<GenomicSegment>& segments,
                                                       const Observations& logrs,
                                                       sex_predict::Sex sex);
f64 SegmentSummarizedLogRLikelihood(
    const Eigen::VectorXd& logrs, f64 sd, f64 purity, f64 tumor_ploidy, f64 total_copy_number, bool expect_haploid);
f64 GetMeanSegmentLogRSD(const vec<GenomicSegment>& segments,
                         const Observations& logrs,
                         const IntervalTrees& logr_trees);
f64 GetMeanSegmentMBAFSD(const vec<GenomicSegment>& segments,
                         const Observations& ref_ads,
                         const Observations& alt_ads,
                         IntervalTrees& variant_trees);
Eigen::VectorXd& SegmentAllSummarizedLogRLikelihoods(
    const Eigen::VectorXd& logrs, f64 sd, f64 purity, f64 tumor_ploidy, bool expect_haploid, Eigen::VectorXd& ret);

// Best germline log-ratio fit for a segment: the maximum germline log-ratio likelihood and the copy-number state that
// achieves it.
struct GermlineLikelihoodResult {
  f64 likelihood;
  size_t total_copy_number;
};

// Germline log-ratio likelihood at fixed purity = kGermlinePurity, ploidy = kDiploidPloidy. The diploid bias prior
// is applied across all copy-number states and the copy-neutral state (CN = 2) is forced to -Inf so a germline
// segment can never be classified as copy-neutral. Returns the best germline log-ratio likelihood and its copy number.
GermlineLikelihoodResult ComputeGermlineLikelihood(const Eigen::VectorXd& logrs,
                                                   f64 mean_logr_sd,
                                                   bool expect_haploid,
                                                   Eigen::VectorXd& logr_likelihoods);
// Returns only the best (max) germline log-ratio likelihood; see ComputeGermlineLikelihood for details.
f64 SegmentGermlineSummarizedLogRLikelihood(const Eigen::VectorXd& logrs,
                                            f64 mean_logr_sd,
                                            bool expect_haploid,
                                            Eigen::VectorXd& logr_likelihoods);
SegmentAllSummarizedBAFLikelihoodsResults SegmentAllSummarizedBAFLikelihoods(
    f64 mbaf, f64 sd, f64 purity, f64 total_copy_number, bool expect_haploid);
GenomicSegment& SegmentSerialSummarizedLogRBAFLikelihood(GenomicSegment& seg,
                                                         const Eigen::VectorXd& logrs,
                                                         f64 purity,
                                                         f64 ploidy,
                                                         f64 mean_logr_sd,
                                                         f64 mean_baf_sd,
                                                         bool prior_2n,
                                                         bool expect_haploid,
                                                         Eigen::VectorXd& logr_likelihoods,
                                                         Eigen::VectorXd& baf_likelihoods,
                                                         f64 baf_likelihood_weight);
bool SkipBasedOnModeAndSex(const GenomicSegment& seg, LikelihoodMode mode, sex_predict::Sex sex);
bool ExpectHaploid(const GenomicSegment& seg, LikelihoodMode mode, sex_predict::Sex sex);
}  // namespace xoos::cnc
