#include "likelihood/baf-tools.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <stdexcept>

#include <xoos/log/logging.h>
#include <xoos/types/float.h>
#include <xoos/types/vec.h>
#include <xoos/util/math.h>

#include "util/eigen-math.h"

namespace xoos::cnc {

namespace {
constexpr f64 kLowerQuartile = 0.25;
constexpr f64 kUpperQuartile = 0.75;
constexpr f64 kSilvermanFactor = 0.9;
constexpr f64 kSilvermanExponentDenominator = 5.0;
constexpr f64 kWeightsSumTolerance = 0.00001;
constexpr std::size_t kNumKdeGridPoints = 512;
constexpr f64 kBafMirrorThreshold = 0.5;
}  // namespace

/**
 * @brief calculates bandwidth for kernel density estimation using Silverman's rule of thumb.
  // https://en.wikipedia.org/wiki/Kernel_density_estimation#A_rule-of-thumb_bandwidth_estimator
 * @param x  data for which kernel density estimation is to be calculated
 * @return bandwidth parameter to use for kernel density estimation
 */
f64 GetBandwidthSilvermanRuleOfThumb(const Eigen::VectorXd& x) {
  // Rule of thumb bandwidth estimator
  const auto n = static_cast<f64>(x.size());
  const f64 std_dev = eigen::Stddev(x);
  const Eigen::VectorXd lower_upper_quartiles = eigen::Quantile(x, eigen::FromStdVec({kLowerQuartile, kUpperQuartile}));
  const f64 iqr_coeff = (lower_upper_quartiles[1] - lower_upper_quartiles[0]) / 1.34;
  f64 a = std::min(std_dev, iqr_coeff);

  // Additional handling for a == 0
  if (math::IsCloseToZero(a)) {
    a = std_dev;
    if (math::IsCloseToZero(a)) {
      a = std::abs(x.size() > 0 ? x[0] : 0.0);
      if (math::IsCloseToZero(a)) {
        a = 1.0;
      }
    }
  }

  return kSilvermanFactor * a * std::pow(n, -1.0 / kSilvermanExponentDenominator);
}

/**
 * @brief return the Guassian Kernel Density Estimate of x, weigted by `weights`. Bandwidth is calculated using the
 * rule-of-thumb bandwidth estimator
 * https://en.wikipedia.org/wiki/Kernel_density_estimation#A_rule-of-thumb_bandwidth_estimator
 * Binned KDE performed using procedure described in https://www.sfu.ca/sasdoc/sashtml/stat/chap33/sect12.htm
 * @param x data for which to calculate kernel density estimates
 * @param weights weights for each data in x
 * @param bandwidth_multiplier scalar applied to the Silverman rule-of-thumb bandwidth; default 1.0; must be > 0
 * @return KernelDensityEstimateRet -
 *          .x BAF bins at which density was estiamted. The value of each BAF bin is the weighted average of all the BAF
 * values in the same bin .y the density at each bin
 */
KernelDensityEstimateRet GaussianKernelDensityEstimate(const Eigen::VectorXd& x, const Eigen::VectorXd& weights) {
  return GaussianKernelDensityEstimate(x, weights, 1.0);
}

KernelDensityEstimateRet GaussianKernelDensityEstimate(const Eigen::VectorXd& x,
                                                       const Eigen::VectorXd& weights,
                                                       const f64 bandwidth_multiplier) {
  const f64 weights_sum = eigen::Sum(weights);
  const f64 diff_from_expected = fabs(weights_sum - static_cast<f64>(x.size()));
  if (weights.size() != 0 && diff_from_expected > kWeightsSumTolerance) {
    Logging::Error("GaussianKernelDensityEstimate: weights are expected to sum to x.size()!!!");
    throw std::runtime_error("weights do not sum to 1");
  }
  // Sort the x values and weights
  const eigen::VectorXu64 sorted_indices = eigen::SortIndex(x);
  const Eigen::VectorXd sorted_x = eigen::Elem(x, sorted_indices);
  const Eigen::VectorXd sorted_weights =
      weights.size() == 0 ? (Eigen::VectorXd::Ones(sorted_x.size()) / static_cast<f64>(sorted_x.size())).eval()
                          : eigen::Elem(weights, sorted_indices);
  // Bandwidth for Gaussian kernel
  const f64 bandwidth = GetBandwidthSilvermanRuleOfThumb(x) * bandwidth_multiplier;
  const auto n = static_cast<f64>(x.size());
  const f64 coef = 1 / (n * bandwidth);

  // Calculate extension and new min/max for kde_x
  const f64 extension = 3.0 * bandwidth;
  const f64 min_x = sorted_x.minCoeff() - extension;
  const f64 max_x = sorted_x.maxCoeff() + extension;

  Eigen::VectorXd density = Eigen::VectorXd::Zero(kNumKdeGridPoints);
  Eigen::VectorXd kde_x = eigen::Linspace(min_x, max_x, kNumKdeGridPoints);
  for (Eigen::Index i = 0; i < kde_x.size(); ++i) {
    const Eigen::VectorXd arg = ((kde_x(i) - sorted_x.array()) / bandwidth).matrix();
    const Eigen::VectorXd weighted = sorted_weights.cwiseProduct(eigen::Normpdf(arg));
    density(i) = eigen::Sum(weighted);
  }
  density = density * coef;
  return KernelDensityEstimateRet{.x = kde_x, .y = density};
}

/**
 * @brief calculate BAF values given ref and alt allelic depths
 * @param ref_ads  ref depths for a list of alleles
 * @param alt_ads  alt depths for the same list of alleles
 * @return std::vector of B-allele fractions (alt depth over total depth) - unmirrored.
 */
Eigen::VectorXd GetBAFFromAD(const Eigen::VectorXd& ref_ads, const Eigen::VectorXd& alt_ads) {
  Eigen::VectorXd bafs(ref_ads.size());
  // calculate BAF from ref_ads and alt_ads
  for (Eigen::Index i = 0; i < ref_ads.size(); ++i) {
    const f64 total_depth = ref_ads(i) + alt_ads(i);
    if (math::IsCloseToZero(total_depth)) {
      bafs(i) = NAN;
    } else {
      bafs(i) = alt_ads(i) / total_depth;
    }
  }
  return bafs;
}

/**
 * @brief find indexes local maxima in a given array. If there's a "plateau" in the data where the left side of the
 * plateau is increasing and the right side is decreasing, count the median of the plateau as a local maxima.
 * @param x std::vector in which to find local maxima
 * @return std::vector of indices of local maxima
 */
vec<size_t> FindLocalMaxima(const Eigen::VectorXd& x) {
  // state 0 = increasing
  // state 1 = no change (plateau)
  // state 2 = decreasing
  size_t state = 0;
  vec<size_t> local_maxima_idxs;
  if (x[0] > x[1]) {
    local_maxima_idxs.push_back(0);
    state = 2;
  }
  size_t plateau_start_idx = 0;
  bool plateau_possibly_maximum = false;
  for (size_t i = 1; i < static_cast<size_t>(x.size()) - 1; ++i) {
    const auto ei = static_cast<Eigen::Index>(i);
    // previously increasing
    if (state == 0) {
      if (math::IsEqualWithTolerance(x[ei], x[ei - 1])) {
        state = 1;
        plateau_start_idx = i - 1;
        plateau_possibly_maximum = true;
        // immediately known to be a local maximum
      } else if (x[ei] < x[ei - 1]) {
        state = 2;
        local_maxima_idxs.push_back(i - 1);
      }  // otherwise keep state at 0 and do nothing
      // previously "equal"
    } else if (state == 1) {
      if (x[ei] < x[ei - 1]) {
        // now decreasing, we can set the middle of the plateau as a local maximum
        state = 2;
        if (plateau_possibly_maximum) {
          const size_t idx = plateau_start_idx + ((i - plateau_start_idx) / 2);
          local_maxima_idxs.push_back(idx);
        }
      } else if (x[ei] > x[ei - 1]) {
        // do nothing else since this is not a maximum
        state = 0;
      }
    } else if (state == 2) {
      if (x[ei] > x[ei - 1]) {
        // now increasing, but we don't care about minima here so do nothing
        state = 0;
      } else if (math::IsEqualWithTolerance(x[ei], x[ei - 1])) {
        // plateau, do nothing
        state = 1;
        plateau_possibly_maximum = false;
      }  // otherwise keep state at 2 and do nothing
    }
  }
  return local_maxima_idxs;
}

/**
 * @brief find the peak of a BAF distribution. This is meant as an alternative
 * to getting the mean BAF after mirroring. By finding the peak of the
 * distribution and then mirroring, we provide a more accurate estimation of the
 * mean mirrored baf
 * @param ref_ads
 * @param alt_ads
 * @return
 */
f64 GetPeakOfBAFDistributionThenMirror(const Eigen::VectorXd& ref_ads, const Eigen::VectorXd& alt_ads) {
  const Eigen::VectorXd bafs = GetBAFFromAD(ref_ads, alt_ads);
  // add in 1-baf in order to f64 the data for calculating density (pronounce the curves more)
  const Eigen::VectorXd aafs = (1 - bafs.array()).matrix();
  Eigen::VectorXd kde_input = eigen::JoinCols(bafs, aafs);
  const eigen::VectorXu64 sort_order = eigen::SortIndex(kde_input);
  // weight the aaf+baf observations by how many reads overlap that observation
  Eigen::VectorXd total_depth = ref_ads + alt_ads;
  total_depth = eigen::JoinCols(total_depth, total_depth);
  Eigen::VectorXd weights = total_depth / (eigen::Sum(total_depth)) * static_cast<f64>(total_depth.size());
  // sort all the input
  kde_input = eigen::Elem(kde_input, sort_order);
  weights = eigen::Elem(weights, sort_order);
  total_depth = eigen::Elem(total_depth, sort_order);
  auto kde = GaussianKernelDensityEstimate(kde_input, weights);
  // for now, just simply calculate all local maxima, and then take the max of those.
  // if needed, we could find "strong" maxima, but at this moment I think that is overkill
  auto local_maxima = FindLocalMaxima(kde.y);
  // If failed to find local maxima
  if (local_maxima.empty()) {
    local_maxima.push_back(eigen::IndexMax(kde.y));
  }
  // we should be suspicious of any segment that has more than 2 peaks - could indicate a segment that could be broken
  // down even further
  if (local_maxima.size() > 2) {
    Logging::Debug("Segment has more than 2 peaks in its BAF distribution");
  }
  size_t global_maximum_idx = 0;
  for (const auto local_maximum_idx : local_maxima) {
    if (kde.y[static_cast<Eigen::Index>(local_maximum_idx)] > kde.y[static_cast<Eigen::Index>(global_maximum_idx)]) {
      global_maximum_idx = local_maximum_idx;
    }
  }
  f64 max_mirrored_baf = kde.x[static_cast<Eigen::Index>(global_maximum_idx)];
  if (max_mirrored_baf > kBafMirrorThreshold) {
    // if the maximum is above 0.5, we should mirror it
    max_mirrored_baf = 1 - max_mirrored_baf;
  }
  return max_mirrored_baf;
}

f64 WeightedMBAFSD(const Eigen::VectorXd& ref_ads, const Eigen::VectorXd& alt_ads, const f64 mean_baf) {
  f64 num = 0;
  f64 denom = 0;
  for (Eigen::Index i = 0; i < ref_ads.size(); ++i) {
    const f64 minor_depth = std::min({ref_ads[i], alt_ads[i]});
    const f64 total_depth = ref_ads[i] + alt_ads[i];
    if (minor_depth < 0) {
      Logging::Error("Found a reference or alternative allele depth that is < 0");
      throw std::runtime_error("Found a reference or alternative allele depth that is < 0");
    }
    num += sqrt(total_depth) * pow((minor_depth / total_depth) - mean_baf, 2);
    denom += sqrt(total_depth);
  }
  return sqrt(num / (denom - 1));
}

}  // namespace xoos::cnc
