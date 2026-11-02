#pragma once
#include <Eigen/Dense>
#include <optional>
#include <string>
#include <tuple>

#include <xoos/io/metadata-util.h>
#include <xoos/sex_predict/sex.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "cli/copy-number-caller-options.h"
#include "likelihood/total-copy-number-prior.h"
#include "observations/observations.h"
#include "segmentation/genomic-segments.h"

namespace xoos::cnc {

// Forward declaration: the per-segment likelihood diagnostic is emitted inline from the grid loops. A null pointer
// disables emission, so the diagnostic adds no work on the default (flag-off) path.
class SegmentLikelihoodWriter;

// Purity/ploidy grid axis definitions. Single source of truth for the grid the search iterates and for any consumer
// (e.g. the per-segment likelihood diagnostic) that must reconstruct the same purity/ploidy values from a cell index.
// Step size for purity
constexpr f64 kLikelihoodRoStep = 0.01;
// Minimum purity to consider
constexpr f64 kLikelihoodRoMin = 0.15;
// Maximum purity to consider
constexpr f64 kLikelihoodRoMax = 0.95;
constexpr f64 kLikelihoodRoRange = kLikelihoodRoMax - kLikelihoodRoMin;
constexpr size_t kLikelihoodRoNumSteps = static_cast<size_t>(kLikelihoodRoRange / kLikelihoodRoStep) + 1;
// Step size for ploidy
constexpr f64 kLikelihoodPhiStep = 0.1;
// Minimum ploidy to consider
constexpr f64 kLikelihoodPhiMin = 1.4;
// Maximum ploidy to consider
constexpr f64 kLikelihoodPhiMax = 6.0;
constexpr f64 kLikelihoodPhiRange = kLikelihoodPhiMax - kLikelihoodPhiMin;
constexpr size_t kLikelihoodPhiNumSteps = static_cast<size_t>(kLikelihoodPhiRange / kLikelihoodPhiStep) + 1;

struct PurityPloidySearchOut {
  f64 purity = 0;
  f64 ploidy = 0;
  Eigen::MatrixXd joint_ll_grid;
  Eigen::MatrixXd logr_ll_grid;
  Eigen::MatrixXd baf_ll_grid;
  TotalCopyNumberPrior total_copy_number_prior = TotalCopyNumberPrior::kDiploidBias;
};

// Shared scalar configuration for the purity/ploidy grid search and its per-metric sub-grids. Bundled to keep the
// grid/search function parameter counts within limits. seg_min_num_snps is unused by the logR-only path. The
// likelihood_writer pointer is optional: nullptr disables per-segment diagnostic emission.
struct PurityPloidyGridParams {
  f64 seg_max_logr = 0;
  size_t seg_min_num_logrs = 0;
  size_t seg_min_num_snps = 0;
  sex_predict::Sex sex = sex_predict::Sex::kUnknown;
  bool use_diploid_bias_prior = false;
  bool tumor_only = false;
  // Weight applied to the BAF term of the joint likelihood (joint = logr + weight * baf). Only the WGD re-run sets a
  // non-default value; all other grid builds keep 1.0 (no change to the joint likelihood).
  f64 baf_likelihood_weight = 1.0;
  SegmentLikelihoodWriter* likelihood_writer = nullptr;
};

const std::string kPurityPloidyGridOutPurity = "Purity";
const std::string kPurityPloidyGridOutPloidy = "Ploidy";
const std::string kPurityPloidyGridOutTotalJointLikelihood = "TotalJointLikelihood";
const std::string kPurityPloidyGridOutTotalLogRatioLikelihood = "TotalLogRatioLikelihood";
const std::string kPurityPloidyGridOutTotalBAlleleFrequencyLikelihood = "TotalBAlleleFrequencyLikelihood";
const std::string kPurityPloidyOutPurity = "Purity";
const std::string kPurityPloidyOutPloidy = "Ploidy";

void PurityPloidySearchMain(const CopyNumberCallerOptions& options);
PurityPloidySearchOut PurityPloidySearch(const vec<GenomicSegment>& segments,
                                         const Observations& logrs,
                                         const std::optional<Observations>& ref_obvs,
                                         const std::optional<Observations>& alt_obvs,
                                         const PurityPloidyGridParams& params,
                                         const fs::path& grid_out_fname,
                                         const io::CommandLineInfo& command_line_info);
PurityPloidySearchOut PurityPloidySearchWithWGDQC(const vec<GenomicSegment>& segments,
                                                  const Observations& logrs,
                                                  const std::optional<Observations>& ref_obvs,
                                                  const std::optional<Observations>& alt_obvs,
                                                  const PurityPloidySearchOptions& options,
                                                  const fs::path& grid_out_fname,
                                                  const std::string& output_prefix,
                                                  sex_predict::Sex sex,
                                                  const io::CommandLineInfo& command_line_info);
Eigen::MatrixXd GetPurityPloidyGridLogROnly(const vec<GenomicSegment>& segments,
                                            const Observations& obvs,
                                            const PurityPloidyGridParams& params);
std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd> GetPurityPloidyGridSerialLogRBAF(
    const vec<GenomicSegment>& segments,
    const Observations& logrs,
    const Observations& ref_ads,
    const Observations& alt_ads,
    const PurityPloidyGridParams& params);
std::tuple<f64, f64> EstimatePurityPloidyByTotalMaximum(const Eigen::MatrixXd& purity_ploidy_grid);
void WritePurityPloidyGrid(const Eigen::MatrixXd& purity_ploidy_grid_joint_lls,
                           const std::optional<Eigen::MatrixXd>& purity_ploidy_grid_logr_lls,
                           const std::optional<Eigen::MatrixXd>& purity_ploidy_grid_baf_lls,
                           const fs::path& out_fname,
                           const io::CommandLineInfo& command_line_info,
                           std::optional<TotalCopyNumberPrior> prior);
void WritePurityPloidy(f64 purity,
                       f64 ploidy,
                       const fs::path& out_fname,
                       const io::CommandLineInfo& command_line_info,
                       std::optional<TotalCopyNumberPrior> prior);
bool ShouldSkipSegment(const GenomicSegment& seg,
                       f64 seg_max_logr,
                       size_t seg_min_num_logrs,
                       std::optional<size_t> seg_min_num_snps,
                       bool check_seg_min_snps = true);
bool HasLowIntervalDensity(const vec<size_t>& logr_idxs, const Observations& logrs, const GenomicSegment& seg);

}  // namespace xoos::cnc
