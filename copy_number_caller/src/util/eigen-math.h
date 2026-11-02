#pragma once

// Eigen numerical helpers for the copy number caller: reductions, quantiles,
// sorting, and other free functions with no direct Eigen equivalent.
// eigen-math-tests.cpp pins every helper against reference output on shared
// inputs. This header is PERMANENT — it is the module's Eigen numeric toolbox.

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::cnc::eigen {

// Column vector of unsigned 64-bit indices, used for genomic start/end
// coordinates. Eigen has no built-in unsigned-index vector typedef, so we
// define one here.
using VectorXu64 = Eigen::Matrix<u64, Eigen::Dynamic, 1>;

// Half, used in several arithmetic shims below. Named to satisfy
// cppcoreguidelines-avoid-magic-numbers while keeping the value exact.
constexpr f64 kHalf = 0.5;

// Stride for the interleaved pairwise accumulation summation order (two
// accumulators advanced two elements at a time).
constexpr Eigen::Index kPairStride = 2;

// Minimum element count for a defined unbiased variance; returns 0 below.
constexpr Eigen::Index kMinVarianceCount = 2;

/// Deterministic summation (non-fast-math path): two interleaved accumulators
/// summed pairwise, then combined. Eigen's own reductions use a different
/// (pairwise/vectorized) order, and plain left-to-right summation also differs;
/// either introduces last-bit rounding differences that downstream RNG-based
/// segmentation amplifies into different segment boundaries. This fixed
/// summation order keeps segment boundaries stable.
inline f64 ArmaAccumulate(const Eigen::VectorXd& v) {
  const Eigen::Index n = v.size();
  f64 acc1 = 0.0;
  f64 acc2 = 0.0;
  Eigen::Index j = 1;
  for (; j < n; j += kPairStride) {
    acc1 += v[j - 1];
    acc2 += v[j];
  }
  if ((j - 1) < n) {
    acc1 += v[j - 1];
  }
  return acc1 + acc2;
}

/// Sample mean.
inline f64 Mean(const Eigen::VectorXd& v) {
  return v.size() == 0 ? 0.0 : ArmaAccumulate(v) / static_cast<f64>(v.size());
}

/// Sum of elements.
inline f64 Sum(const Eigen::VectorXd& v) {
  return ArmaAccumulate(v);
}

/// Unbiased sample variance (normalised by N-1): mean via ArmaAccumulate, then
/// interleaved accumulation of the squared deviations (acc2) and their sum
/// (acc3), combined as (acc2 - acc3^2/n)/(n-1). The fixed arithmetic order
/// keeps segment std-dev stable. Returns 0 for N<2.
inline f64 Var(const Eigen::VectorXd& v) {
  const Eigen::Index n = v.size();
  if (n < kMinVarianceCount) {
    return 0.0;
  }
  const f64 sample_mean = Mean(v);
  f64 acc2 = 0.0;
  f64 acc3 = 0.0;
  Eigen::Index j = 1;
  for (; j < n; j += kPairStride) {
    const f64 tmpi = sample_mean - v[j - 1];
    const f64 tmpj = sample_mean - v[j];
    acc2 += (tmpi * tmpi) + (tmpj * tmpj);
    acc3 += tmpi + tmpj;
  }
  if ((j - 1) < n) {
    const f64 tmpi = sample_mean - v[j - 1];
    acc2 += tmpi * tmpi;
    acc3 += tmpi;
  }
  return (acc2 - ((acc3 * acc3) / static_cast<f64>(n))) / static_cast<f64>(n - 1);
}

/// Unbiased sample standard deviation (normalised by N-1).
inline f64 Stddev(const Eigen::VectorXd& v) {
  return std::sqrt(Var(v));
}

/// Median: mean of the two middle values for even N.
inline f64 Median(const Eigen::VectorXd& v) {
  const Eigen::Index n = v.size();
  if (n == 0) {
    return 0.0;
  }
  vec<f64> s(v.begin(), v.end());
  std::ranges::sort(s);
  constexpr size_t kTwo = 2;
  const size_t mid = static_cast<size_t>(n) / kTwo;
  if (static_cast<size_t>(n) % kTwo == 1) {
    return s[mid];
  }
  // val1 = middle element (s[mid]), val2 = max of the lower half (s[mid-1]),
  // combined as A + (B-A)/2. std::lerp is deliberately NOT used here — it may
  // contract to an FMA and change the last bit, which eigen-math-tests.cpp pins.
  const f64 val1 = s[mid];
  const f64 val2 = s[mid - 1];
  return val1 + ((val2 - val1) * kHalf);  // NOSONAR std::lerp may contract to FMA and change last bit
}

/// Sample quantiles (Hyndman & Fan type 5, alpha=0.5). Input need not be sorted.
inline Eigen::VectorXd Quantile(const Eigen::VectorXd& v, const Eigen::VectorXd& probs) {
  Eigen::VectorXd out(probs.size());
  const Eigen::Index n = v.size();
  if (n == 0) {
    out.setConstant(std::numeric_limits<f64>::quiet_NaN());
    return out;
  }
  vec<f64> y(v.begin(), v.end());
  std::ranges::sort(y);
  const auto big_n = static_cast<f64>(n);
  // type 5
  constexpr f64 kAlpha = kHalf;
  const f64 p_min = (1.0 - kAlpha) / big_n;
  const f64 p_max = (big_n - kAlpha) / big_n;
  for (Eigen::Index i = 0; i < probs.size(); ++i) {
    const f64 p_i = probs(i);
    f64 out_val = 0.0;
    if (p_i < p_min) {
      out_val = (p_i < 0.0) ? -std::numeric_limits<f64>::infinity() : y.front();
    } else if (p_i > p_max) {
      out_val = (p_i > 1.0) ? std::numeric_limits<f64>::infinity() : y.back();
    } else {
      const auto k = static_cast<size_t>(std::floor((big_n * p_i) + kAlpha));
      const f64 p_k = (static_cast<f64>(k) - kAlpha) / big_n;
      const f64 w = (p_i - p_k) * big_n;
      const f64 y_lo = y[k - 1];
      // When p_i lands exactly on the upper interior bound, k == n and y[k]
      // would read one past the end; w is 0 there so the upper sample is unused.
      const f64 y_hi = (k < static_cast<size_t>(n)) ? y[k] : y_lo;
      out_val = ((1.0 - w) * y_lo) + (w * y_hi);
    }
    out(i) = out_val;
  }
  return out;
}

/// Cumulative sum.
inline Eigen::VectorXd Cumsum(const Eigen::VectorXd& v) {
  Eigen::VectorXd out(v.size());
  f64 acc = 0.0;
  for (Eigen::Index i = 0; i < v.size(); ++i) {
    acc += v(i);
    out(i) = acc;
  }
  return out;
}

/// Trapezoidal integral of y over x (returns scalar).
inline f64 Trapz(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
  f64 acc = 0.0;
  for (Eigen::Index i = 1; i < x.size(); ++i) {
    acc += kHalf * (x(i) - x(i - 1)) * (y(i) + y(i - 1));
  }
  return acc;
}

/// Ascending sort indices: the permutation that sorts v ascending. NaNs are
/// placed after all non-NaN values. Ties (including NaN-NaN) are resolved by
/// original index (lowest-index-first), making the output deterministic.
inline Eigen::Matrix<u64, Eigen::Dynamic, 1> SortIndex(const Eigen::VectorXd& v) {
  const auto n = static_cast<size_t>(v.size());
  vec<u64> idx(n);
  std::ranges::iota(idx, 0U);
  std::ranges::sort(idx, [&v](const u64 a, const u64 b) {
    const f64 left = v(static_cast<Eigen::Index>(a));
    const f64 right = v(static_cast<Eigen::Index>(b));
    const bool left_nan = std::isnan(left);
    const bool right_nan = std::isnan(right);
    if (left_nan != right_nan) {
      return !left_nan;
    }
    if (left < right) {
      return true;
    }
    if (right < left) {
      return false;
    }
    return a < b;
  });
  Eigen::Matrix<u64, Eigen::Dynamic, 1> out(v.size());
  for (size_t i = 0; i < n; ++i) {
    out(static_cast<Eigen::Index>(i)) = idx[i];
  }
  return out;
}

/// Index of the maximum element (first occurrence on ties). Eigen's
/// maxCoeff(&idx) keeps the first index whose value is strictly greater than
/// the running maximum, so ties resolve to the earliest index. Empty input is
/// undefined.
inline u64 IndexMax(const Eigen::VectorXd& v) {
  Eigen::Index idx = 0;
  v.maxCoeff(&idx);
  return static_cast<u64>(idx);
}

/// Linearly spaced vector of `num` points in [start, end].
inline Eigen::VectorXd Linspace(const f64 start, const f64 end, const u64 num) {
  return Eigen::VectorXd::LinSpaced(static_cast<Eigen::Index>(num), start, end);
}

/// (row, col) of the maximum matrix element. Scans in column-major order
/// (column outer, row inner) keeping the first position whose value is strictly
/// greater than the running maximum, so on ties the earliest column-major
/// position wins. Empty input is undefined.
inline std::pair<u64, u64> ArgMaxRowCol(const Eigen::MatrixXd& m) {
  Eigen::Index best_row = 0;
  Eigen::Index best_col = 0;
  f64 best = m(0, 0);
  for (Eigen::Index col = 0; col < m.cols(); ++col) {
    for (Eigen::Index row = 0; row < m.rows(); ++row) {
      if (m(row, col) > best) {
        best = m(row, col);
        best_row = row;
        best_col = col;
      }
    }
  }
  return {static_cast<u64>(best_row), static_cast<u64>(best_col)};
}

/// Replace every NaN entry of a matrix with `value`. Non-NaN entries (including
/// +/-inf) are left untouched.
inline void ReplaceNaN(Eigen::MatrixXd& m, const f64 value) {
  m = m.array().isNaN().select(value, m);
}

/// Indices of finite (non-NaN, non-Inf) elements.
inline Eigen::Matrix<u64, Eigen::Dynamic, 1> FindFinite(const Eigen::VectorXd& v) {
  vec<u64> found;
  for (Eigen::Index i = 0; i < v.size(); ++i) {
    if (std::isfinite(v(i))) {
      found.push_back(static_cast<u64>(i));
    }
  }
  Eigen::Matrix<u64, Eigen::Dynamic, 1> out(static_cast<Eigen::Index>(found.size()));
  for (size_t i = 0; i < found.size(); ++i) {
    out(static_cast<Eigen::Index>(i)) = found[i];
  }
  return out;
}

// Constants for normpdf / normcdf. sqrt(2*pi) has no std::numbers entry, so it
// is written as an explicit literal.
constexpr f64 kSqrt2 = std::numbers::sqrt2_v<f64>;
// NOLINTNEXTLINE(modernize-use-std-numbers) — no std entry for sqrt(2*pi)
constexpr f64 kSqrt2Pi =
    2.5066282746310005024157652848110452530069867406099383166299235763422936546078419749465958383780572661;

/// Standard normal PDF, elementwise (mean 0, sd 1):
/// exp(-0.5 * x*x) / (sigma * sqrt2pi).
inline Eigen::VectorXd Normpdf(const Eigen::VectorXd& v) {
  Eigen::VectorXd out(v.size());
  for (Eigen::Index i = 0; i < v.size(); ++i) {
    const f64 x = v(i);
    out(i) = std::exp(-kHalf * (x * x)) / (1.0 * kSqrt2Pi);
  }
  return out;
}

/// Standard normal CDF, elementwise (mean 0, sd 1):
/// 0.5 * erfc(x / (sigma * -sqrt2)).
inline Eigen::VectorXd Normcdf(const Eigen::VectorXd& v) {
  Eigen::VectorXd out(v.size());
  for (Eigen::Index i = 0; i < v.size(); ++i) {
    const f64 tmp = v(i) / (1.0 * (-kSqrt2));
    out(i) = kHalf * std::erfc(tmp);
  }
  return out;
}

/// Vertically concatenate two column vectors.
inline Eigen::VectorXd JoinCols(const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
  Eigen::VectorXd out(a.size() + b.size());
  out.head(a.size()) = a;
  out.tail(b.size()) = b;
  return out;
}

/// Vertically concatenate two u64 column vectors (genomic coordinate vectors).
inline VectorXu64 JoinCols(const VectorXu64& a, const VectorXu64& b) {
  VectorXu64 out(a.size() + b.size());
  out.head(a.size()) = a;
  out.tail(b.size()) = b;
  return out;
}

/// Gather elements at the given indices — builds a new vector
/// v[idxs[0]], v[idxs[1]], ... in index order.
inline Eigen::VectorXd Elem(const Eigen::VectorXd& v, const Eigen::Matrix<u64, Eigen::Dynamic, 1>& idxs) {
  Eigen::VectorXd out(idxs.size());
  for (Eigen::Index i = 0; i < idxs.size(); ++i) {
    out(i) = v(static_cast<Eigen::Index>(idxs(i)));
  }
  return out;
}

/// Gather overload for a std::vector of indices (library-neutral call sites).
inline Eigen::VectorXd Elem(const Eigen::VectorXd& v, const vec<u64>& idxs) {
  Eigen::VectorXd out(static_cast<Eigen::Index>(idxs.size()));
  for (size_t i = 0; i < idxs.size(); ++i) {
    out(static_cast<Eigen::Index>(i)) = v(static_cast<Eigen::Index>(idxs[i]));
  }
  return out;
}

/// Construct an Eigen column vector from a std::vector<f64> (element-order copy).
inline Eigen::VectorXd FromStdVec(const vec<f64>& v) {
  Eigen::VectorXd out(static_cast<Eigen::Index>(v.size()));
  for (size_t i = 0; i < v.size(); ++i) {
    out(static_cast<Eigen::Index>(i)) = v[i];
  }
  return out;
}

/// Construct an Eigen u64 column vector (VectorXu64) from a std::vector<u64>
/// (element-order copy), used for genomic start/end coordinates. Named
/// distinctly from FromStdVec to avoid brace-init overload ambiguity.
inline VectorXu64 FromStdVecU64(const vec<u64>& v) {
  VectorXu64 out(static_cast<Eigen::Index>(v.size()));
  for (size_t i = 0; i < v.size(); ++i) {
    out(static_cast<Eigen::Index>(i)) = v[i];
  }
  return out;
}

/// Indices where a boolean mask is true, in ascending order. The caller
/// pre-computes the mask.
inline Eigen::Matrix<u64, Eigen::Dynamic, 1> Find(const Eigen::Matrix<bool, Eigen::Dynamic, 1>& mask) {
  vec<u64> found;
  for (Eigen::Index i = 0; i < mask.size(); ++i) {
    if (mask(i)) {
      found.push_back(static_cast<u64>(i));
    }
  }
  Eigen::Matrix<u64, Eigen::Dynamic, 1> out(static_cast<Eigen::Index>(found.size()));
  for (size_t i = 0; i < found.size(); ++i) {
    out(static_cast<Eigen::Index>(i)) = found[i];
  }
  return out;
}

// ---------------------------------------------------------------------------
// Random number generation
//
// The segmentation permutation tests draw from Armadillo's global RNG via
// `arma::arma_rng::set_seed(i)` followed by `arma::randi(...)`. With
// ARMA_USE_CXX11_RNG (the build config in this repo — see
// armadillo_bits/arma_rng.hpp) that RNG is exactly a `std::mt19937_64`:
//
//   * set_seed(v)                -> producer.seed(v)
//   * randi(distr_param(a, b))   -> std::uniform_int_distribution<int>(a, b)(producer)
//
// The scalar and fill paths both go through
// `std::uniform_int_distribution<int>::operator()(std::mt19937_64&)`
// (armadillo_bits/arma_rng.hpp, arma_rng::randi<eT>::fill, CXX11 branch, and
// fn_randi.hpp which routes every randi(distr_param) call — scalar or vector —
// through that same fill). Because both the Armadillo build and this shim link
// the same libstdc++ in the same Docker image, replaying the identical
// std:: engine and distribution reproduces the draw sequence bit-for-bit, so
// segmentation breakpoints are unchanged. rng-tests.cpp pins this against a
// captured golden sequence. This header is PERMANENT; it stays after Armadillo is removed.
//
// Armadillo's seed_type is std::mt19937_64::result_type (u64); callers pass
// small loop counters, so the value is taken as-is.
class ArmaRng {
 public:
  /// Reproduces arma::arma_rng::set_seed(val): reseed the mt19937_64 producer.
  void SetSeed(const u64 val) {
    _engine.seed(val);
  }

  /// One draw of arma::randi(arma::distr_param(a, b)) — a uniform integer in
  /// the closed range [a, b]. Matches the scalar/fill CXX11 path exactly.
  /// CTAD deduces std::uniform_int_distribution<s32>, matching arma's <int>.
  s32 Randi(const s32 a, const s32 b) {
    std::uniform_int_distribution distr(a, b);
    return distr(_engine);
  }

  /// N draws of arma::randi<arma::uvec>(N, arma::distr_param(a, b)) as a u64
  /// index vector. Armadillo fills element by element from the same producer in
  /// order (arma_rng::randi<eT>::fill), so the sequence matches bit-for-bit.
  VectorXu64 RandiFill(const u64 n, const s32 a, const s32 b) {
    std::uniform_int_distribution distr(a, b);
    VectorXu64 out(static_cast<Eigen::Index>(n));
    for (Eigen::Index i = 0; i < out.size(); ++i) {
      out(i) = static_cast<u64>(distr(_engine));
    }
    return out;
  }

  /// Reproduces arma::shuffle(v) for a column vector exactly (op_shuffle_meat.hpp,
  /// vector / non-alias branch). Armadillo does NOT use Fisher-Yates here; it
  /// builds one sort key per element from the unbounded randi<int>() draw
  ///     val[i] = int(double(producer()) * scale),  scale = INT_MAX / mt19937_64::max()
  /// (arma_rng::randi<int>, CXX11 branch), assigns index[i] = i, sorts the
  /// (val, index) packets ascending by val only (arma_sort_index_helper_ascend,
  /// a non-stable std::sort), then gathers out[i] = v[index[i]]. Matching the
  /// key generation, the val-only comparator, and std::sort's ordering
  /// reproduces the permutation bit-for-bit. rng-tests.cpp pins this against
  /// a captured golden sequence.
  Eigen::VectorXd Shuffle(const Eigen::VectorXd& v) {
    const Eigen::Index n = v.size();
    Eigen::VectorXd out(n);
    if (n == 0) {
      return out;
    }
    constexpr f64 kScale = static_cast<f64>(std::numeric_limits<int>::max()) / static_cast<f64>(std::mt19937_64::max());

    struct Packet {
      s32 val;  // matches arma's int sort key exactly (s32 is a typedef for int)
      Eigen::Index index;
    };

    std::vector<Packet> packets(static_cast<size_t>(n));
    for (Eigen::Index i = 0; i < n; ++i) {
      packets[static_cast<size_t>(i)].val = static_cast<s32>(static_cast<f64>(_engine()) * kScale);
      packets[static_cast<size_t>(i)].index = i;
    }
    std::ranges::sort(packets, [](const Packet& a, const Packet& b) { return a.val < b.val; });
    for (Eigen::Index i = 0; i < n; ++i) {
      out(i) = v(packets[static_cast<size_t>(i)].index);
    }
    return out;
  }

 private:
  // Fixed default seed is intentional: Armadillo's arma_rng starts from a
  // deterministic state until set_seed() is called, and reproducing that exact
  // sequence is the whole point of this shim. Callers reseed via SetSeed().
  // NOLINTNEXTLINE(cert-msc51-cpp,cert-msc32-c) - deterministic seed reproduces arma_rng semantics by design
  std::mt19937_64 _engine{std::mt19937_64::default_seed};
};

}  // namespace xoos::cnc::eigen
