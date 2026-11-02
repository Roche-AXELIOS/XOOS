#include "xoos/gc_correct/loess.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/util/math.h>

namespace xoos::gc_correct {

/**
 * @brief Get the range of the top N closest elements to x in arr. Assume arr is sorted
 * @param arr
 * @param n
 * @param x
 * @return the first and last indices in arr representing N closest elements in arr to x. [first, last] is inclusive,
 * closed
 */
std::tuple<Eigen::Index, Eigen::Index> GetNClosestRange(const Eigen::VectorXd& arr, const Eigen::Index n, const f64 x) {
  // since arr is sorted, arr - x will be sorted too, but abs(arr - x) will not be sorted. Will have to find location of
  // min, then do a search from there
  Eigen::Index first = 0;
  Eigen::Index last = arr.size() - 1;
  while (last - first >= n) {
    // if the distance`first` is greater than the distance at `last`, then we can move `first` up
    if (std::abs(arr[first] - x) > std::abs(arr[last] - x)) {
      ++first;
    } else {
      --last;
    }
  }
  return {first, last};
}

namespace {
// Exponent used by the tricubic weight function (1 - |d|^3)^3.
constexpr int kTricubicExponent = 3;

Eigen::VectorXd Tricubic(const Eigen::VectorXd& arr) {
  return (1.0 - arr.array().pow(kTricubicExponent)).pow(kTricubicExponent);
}
}  // namespace

Eigen::VectorXd TricubicDistWeights(const Eigen::VectorXd& v, const f64 x) {
  const Eigen::VectorXd distances = (v.array() - x).abs();
  const f64 max_distance = distances.maxCoeff();
  // When all points coincide (max distance is 0), there is no relative distance to weight by, so every point gets the
  // full weight of 1.0.
  if (math::IsCloseToZero(max_distance)) {
    return Eigen::VectorXd::Ones(v.size());
  }
  const Eigen::VectorXd normed_distances = distances / max_distance;
  return Tricubic(normed_distances);
}

namespace {
Eigen::MatrixXd MakePoly(const Eigen::VectorXd& x, const size_t degree) {
  Eigen::MatrixXd x_poly(x.size(), degree + 1);
  // fill each column (idx=c) with successive polynomials d of x, 0<=d<=2
  for (Eigen::Index i = 0; i < ToSigned(degree + 1); ++i) {
    x_poly.col(i) = x.array().pow(static_cast<f64>(i));
  }
  return x_poly;
}
}  // namespace

// x/y (known points) and new_x (targets) are conceptually distinct despite sharing a type; reordering this
// public API is out of scope for a Sonar-driven cleanup.
Eigen::VectorXd WeightedPolyFit(const Eigen::VectorXd& x,  // NOSONAR(cpp:S5419)
                                const Eigen::VectorXd& y,
                                const Eigen::VectorXd& new_x,
                                const Eigen::VectorXd& weights,
                                const size_t degree) {
  const Eigen::MatrixXd x_poly = MakePoly(x, degree);
  Eigen::VectorXd beta(degree + 1);
  if (weights.size() != 0) {
    // Weighted least squares via row-scaling: scale each row of X and y by sqrt(weight), then solve the
    // equivalent unweighted problem directly on the scaled design matrix. Solving via normal equations
    // (INV(T(X)WX)T(X)WY, as in e.g. https://www.stat.cmu.edu/~cshalizi/uADA/16/lectures/08.pdf Eq 7.4) squares
    // the condition number of X, which amplifies small input perturbations (from upstream coverage/GC data)
    // into much larger errors in beta when X is ill-conditioned — as it can be here, since X's rows come from a
    // narrow local window of the closest GC-bias points. QR avoids ever forming XtWX, so the amplification is
    // bounded by cond(X) instead of cond(X)^2.
    Eigen::MatrixXd x_poly_w = x_poly;
    Eigen::VectorXd y_w = y;
    for (Eigen::Index i = 0; i < weights.size(); ++i) {
      const f64 sqrt_weight = std::sqrt(weights[i]);
      x_poly_w.row(i) *= sqrt_weight;
      y_w[i] *= sqrt_weight;
    }
    beta = x_poly_w.colPivHouseholderQr().solve(y_w);
  } else {
    // Direct QR solve on X instead of forming/inverting XtX (see rationale above).
    beta = x_poly.colPivHouseholderQr().solve(y);
  }
  const Eigen::MatrixXd poly_x = MakePoly(new_x, degree);  // nx3
  return poly_x * beta;                                    // nx3 * 3x1 -> nx1
}

namespace {
// x_to_predict, sorted_xx, sorted_yy are conceptually distinct despite sharing a type; this is a private
// helper with a single call site, so restructuring for this alone isn't worthwhile.
f64 SmoothPoint(const Eigen::VectorXd& x_to_predict,  // NOSONAR(cpp:S5419)
                const Eigen::VectorXd& sorted_xx,
                const Eigen::VectorXd& sorted_yy,
                const Eigen::Index n,
                const s32 degree,
                const Eigen::Index i) {
  auto [start, end] = GetNClosestRange(sorted_xx, n, x_to_predict[i]);
  if (n > 1) {
    const Eigen::VectorXd xx_subvec = sorted_xx.segment(start, end - start + 1);
    const Eigen::VectorXd yy_subvec = sorted_yy.segment(start, end - start + 1);
    const Eigen::VectorXd weights = TricubicDistWeights(xx_subvec, x_to_predict[i]);
    Eigen::VectorXd x_pred(1);
    x_pred << x_to_predict[i];
    const f64 y_pred = WeightedPolyFit(xx_subvec, yy_subvec, x_pred, weights, ToUnsigned(degree))[0];
    return y_pred;
  }
  return sorted_yy[start];
}
}  // namespace

/**
 * @brief Perform LOESS (Locally Estimated Scatterplot Smoothing) regression to smooth noisy data using local polynomial
 * fitting.
 * @param x_to_predict - x values where predictions are desired
 * @param x            - known x values
 * @param y            - known y values
 * @param first_span   - span parameter representing the percent of data used to smooth each point
 * @param degree       - degree of the local polynomial to fit (e.g., 1 for linear, 2 for quadratic)
 * @return Eigen::VectorXd of smoothed y values corresponding to x_to_predict
 * @note x/y (known points) and x_to_predict (prediction targets) are conceptually distinct despite sharing a type;
 * reordering this public API is out of scope for a Sonar-driven cleanup.
 */
Eigen::VectorXd Loess(const Eigen::VectorXd& x_to_predict,  // NOSONAR(cpp:S5419)
                      const Eigen::VectorXd& x,
                      const Eigen::VectorXd& y,
                      const f64 first_span,
                      const s32 degree) {
  // Sorting indices indirectly (rather than sorting x/y in place) keeps x and y aligned without an extra data copy.
  Eigen::VectorXd ret;
  std::vector<Eigen::Index> sort_order(ToUnsigned(x.size()));
  // Seed with sequential indices [0, size) so x can be sorted indirectly.
  std::ranges::iota(sort_order, 0);
  std::ranges::sort(sort_order, [&x](const Eigen::Index a, const Eigen::Index b) { return x[a] < x[b]; });
  Eigen::VectorXd sorted_xx(x.size());
  Eigen::VectorXd sorted_yy(y.size());
  for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(sort_order.size()); ++i) {
    sorted_xx[i] = x[sort_order[ToUnsigned(i)]];
    sorted_yy[i] = y[sort_order[ToUnsigned(i)]];
  }
  const auto n = static_cast<Eigen::Index>(first_span * static_cast<f64>(x.size()));
  if (n == 0) {
    throw std::runtime_error("Too few data points per window. Try increasing the --first-span parameter.");
  }
  ret.resize(x_to_predict.size());
  for (Eigen::Index i = 0; i < x_to_predict.size(); ++i) {
    ret[i] = SmoothPoint(x_to_predict, sorted_xx, sorted_yy, n, degree, i);
  }
  return ret;
}
}  // namespace xoos::gc_correct
