#include "distributions.h"

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>

#include <xoos/types/float.h>
#include <xoos/util/math.h>

#include "external/asa243/asa243.hpp"
#include "stats.h"

namespace xoos::stats {
using std::log;
using std::pow;
using std::tgamma;
using std::numbers::pi;

/**
 * @brief binomial PMF - probability of seeing k successes given n trials, with probability p of seeing a success per
 * trial
 * @param k
 * @param n
 * @param p
 * @return probability of seeing k successes given binomial distribution with paramerters (n, p)
 */
f64 BinomialPMF(f64 k, f64 n, f64 p) {
  return NCombinations(n, k) * pow(p, k) * pow(1 - p, n - k);
}

/**
 * @brief log binomial PMF - natural log of  probability of seeing k successes given n trials, with probability p of
 * seeing a success per trial.
 * @param k
 * @param n
 * @param p
 * @return natural log of probability of seeing k successes given binomial distribution with paramerters (n, p)
 */
f64 LogBinomialPMF(f64 k, f64 n, f64 p) {
  if (math::IsCloseToZero(p) || math::IsEqualWithTolerance(p, 1.0)) {
    throw std::runtime_error("p==" + std::to_string(p) + " not allowed for LogBinomialPMF");
  }
  return LogNCombinations(n, k) + (k * log(p)) + ((n - k) * log(1 - p));
}

f64 NonStandardTCDF(f64 x, f64 df, f64 ncp) {
  const f64 ret = stats::tnc(x, df, ncp);
  return ret;
}

f64 NonStandardTCDF(f64 x, f64 df, f64 mean, f64 sd) {
  const f64 ret = stats::tnc((x - mean) / sd, df, 0) / sd;
  return ret;
}

f64 NonStandardTPDF(f64 x, f64 df, f64 ncp) {
  // Threshold below which x is treated as zero to select the analytic limit form of the PDF.
  constexpr f64 kNearZeroThreshold = 1e-15;
  if (std::abs(x) > kNearZeroThreshold) {
    return df / x * (NonStandardTCDF(x * sqrt(1 + (2 / df)), df + 2, ncp) - NonStandardTCDF(x, df, ncp));
  }
  // Gaussian limit form of the non-central t PDF at x == 0.
  constexpr f64 kHalf = 0.5;
  return tgamma(df + kHalf) / (sqrt(pi * df) * tgamma(df / 2)) * exp(-kHalf * ncp * ncp);
}

f64 NonStandardTPDF(f64 x, f64 df, f64 mean, f64 sd) {
  return NonStandardTPDF((x - mean) / sd, df, 0) / sd;
}

}  // namespace xoos::stats
