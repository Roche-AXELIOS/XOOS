#include "segmentation/cbs-util.h"

#include <Eigen/Dense>
#include <cmath>
#include <cstddef>
#include <utility>

#include <xoos/types/float.h>
#include <xoos/types/int.h>

#include "util/eigen-math.h"

namespace xoos::cnc {
/**
 * @brief Imagine segment as circle joined at both ends. This function calculates t-statistic for testing the hypothesis
 * that arc from i+1 to j and its complement have different means.
 * @param n_v1 size of minor arc
 * @param mean_v1 mean of minor arc
 * @param var_v1 variance of minor arc
 * @param n_v2 size of major arc
 * @param mean_v2 mean of major arc
 * @param var_v2 variance of major arc
 * @return t-statistic describing difference in means between minor and major arc
 */
f64 TStat(const f64 n_v1, const f64 mean_v1, const f64 var_v1, const f64 n_v2, const f64 mean_v2, const f64 var_v2) {
  // pooled variance of v1 and v2
  const f64 s_ij = std::pow(((n_v1 - 1) * var_v1 + (n_v2 - 1) * var_v2) / (n_v1 + n_v2 - 2), 0.5);
  // t stat is a function of the means of both std::vectors and their pooled variance
  const f64 t_ij = (mean_v1 - mean_v2) / (s_ij * std::pow(((1.0 / n_v1) + (1.0 / n_v2)), 0.5));
  return t_ij;
}

/**
 * @brief calculates t-statistic between two std::vectors that can be of unequal size and variance
 * @param v1 std::vector 1
 * @param v2 std::vector 2
 * @return t-statistic between v1 and v2
 */
f64 TStatFromVecs(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2) {
  const auto n_v1 = static_cast<f64>(v1.size());
  const auto n_v2 = static_cast<f64>(v2.size());
  const f64 mean_v1 = eigen::Mean(v1);
  const f64 mean_v2 = eigen::Mean(v2);
  const f64 var_v1 = eigen::Var(v1);
  const f64 var_v2 = eigen::Var(v2);
  return TStat(n_v1, mean_v1, var_v1, n_v2, mean_v2, var_v2);
}

/* warning - this is the biased estimator for variance. Do we need an unbiased estimator??? */
f64 Variance(const f64 m, const f64 mean, const f64 sq_sum) {
  return (sq_sum / m) - std::pow(mean, 2);
}

Eigen::VectorXd& FisherYatesShuffleInPlace(Eigen::VectorXd& arr, eigen::ArmaRng& rng) {
  for (s32 i = static_cast<s32>(arr.size()) - 1; i > 0; --i) {
    // get a random index in [0, i]
    const auto random_index = static_cast<Eigen::Index>(rng.Randi(0, i));
    std::swap(arr[i], arr[random_index]);
  }
  return arr;
}

/**
 * @brief return a t-stat that the difference in means between a "within" sub-segment (specified by i,j) against its
 * complementary "without" segment
 *
 * sub-segmnets. this generally works if we're looking for small sub-segments within very large segments
 * @param psums
 * @param i
 * @param j
 * @return
 */
f64 TStatSquaredMean(const Eigen::VectorXd& psums, const size_t i, const size_t j) {
  const auto inner_length = static_cast<f64>(j - i);
  const auto outer_length = static_cast<f64>(psums.size() - j + i);
  const auto ei = static_cast<Eigen::Index>(i);
  const auto ej = static_cast<Eigen::Index>(j);
  const f64 sum1 = i > 0 ? psums[ej - 1] - psums[ei - 1] : psums[ej - 1];
  const f64 sum2 = psums[psums.size() - 1] - sum1;
  return std::pow((sum1 / inner_length) - (sum2 / outer_length), 2) / (1 / inner_length + 1 / outer_length);
}

f64 TStatSquaredMean(const f64 total_sum, const f64 total_length, const f64 minor_sum, const f64 minor_length) {
  const f64 major_sum = total_sum - minor_sum;
  const f64 major_length = total_length - minor_length;
  return std::pow((major_sum / major_length) - (minor_sum / minor_length), 2) / (1 / major_length + 1 / minor_length);
}

f64 TStatSquaredMean(const f64 total_sum, const size_t total_length, const f64 minor_sum, const size_t minor_length) {
  return TStatSquaredMean(total_sum, static_cast<f64>(total_length), minor_sum, static_cast<f64>(minor_length));
}

}  // namespace xoos::cnc
