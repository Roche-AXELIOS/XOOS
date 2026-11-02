/**
 * reference:
 * Olshen et al, 2004. "Circular binary segmentation for the analysis of array‐based DNA copy number data"
 * https://academic.oup.com/biostatistics/article/5/4/557/275197
 */

#include "segmentation/cbs-for-somatic.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <ranges>
#include <stack>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/stats/external/asa243/asa243.hpp>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/util/math.h>

#include "segmentation/boundaries.h"
#include "util/eigen-math.h"

namespace xoos::cnc::somatic {

using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;
using std::chrono::nanoseconds;
using std::chrono::seconds;

const size_t kCbsSiegmundPIters = 100;
const size_t kCbsMinObvsN = 10;
const size_t kCbsHybridK = 25;
constexpr f64 kSiegmundLeadingCoeff = 0.25;
constexpr size_t kPrecomputedBoundaryPermutations = 10000;
constexpr f64 kPrecomputedBoundaryMaxP = 0.01;
const size_t kCbsMaxStackDepth = 1024;
// minimum segment size (as a multiple of kCbsHybridK) required to use the approximate p-value method
const size_t kCbsHybridMinSizeMultiple = 4;

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
namespace {
inline f64 Variance(const f64 m, const f64 mean, const f64 sq_sum) {
  return (sq_sum / m) - std::pow(mean, 2);
}
}  // namespace

/**
 * @brief Gets the sub-segment in obvs that has the maximum t-statistic when compared to the observations outside of the
 * subsegment
 * @param obvs observation std::vector
 * @param max_n maximum sub-segment size to test. Set to 0 if no upper limit is needed
 * @param min_n minimum sub-segment size to test.
 * @return Segment contain start and end indices within obvs of the segment w/ that maximizes the t-statistic
 */
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
Segment GetMaxTSegCumSum(const Eigen::VectorXd& obvs, const size_t max_n, const size_t min_n) {
  f64 max_t = 0;
  f64 max_t_mean = 0;
  f64 max_t_m = 0;
  const auto m = static_cast<u64>(obvs.size());
  u64 max_t_i = 0;
  u64 max_t_j = 0;
  Eigen::VectorXd obvs_psums(eigen::Cumsum(obvs));
  Eigen::VectorXd obvs_psqsums(eigen::Cumsum(obvs.array().square().matrix()));
  for (u64 i = 0; i < m - 1; ++i) {
    auto inner_loop_cond = [m, max_n](const u64 i, const u64 j) {
      if (max_n > 0) {
        // if max_n is specified, then stop when the segment exceeds this mch
        return j < m && j - i < max_n;
      }
      return j < m;
    };
    for (u64 j = i + min_n; inner_loop_cond(i, j); ++j) {  // NOSONAR(cpp:S5311) counter j lives in condition lambda
      const auto m_v1 = static_cast<f64>(j - i);
      const auto m_v2 = static_cast<f64>(i + m - j);
      // skip any time either segment is too small
      if (j - i < min_n || i + m - j < min_n) {
        continue;
      }
      const auto ei = static_cast<Eigen::Index>(i);
      const auto ej = static_cast<Eigen::Index>(j);
      const auto em = static_cast<Eigen::Index>(m);
      const f64 mean_v1 = (obvs_psums[ej - 1] - (i > 0 ? obvs_psums[ei - 1] : 0)) / m_v1;
      const f64 var_v1 = Variance(m_v1, mean_v1, obvs_psqsums[ej - 1] - (i > 0 ? obvs_psqsums[ei - 1] : 0));
      f64 mean_v2 = NAN;
      f64 var_v2 = NAN;
      if (i > 0 && j < m) {
        mean_v2 = (obvs_psums[ei - 1] + obvs_psums[em - 1] - obvs_psums[ej - 1]) / m_v2;
        var_v2 = Variance(m_v2, mean_v2, obvs_psqsums[ei - 1] + obvs_psqsums[em - 1] - obvs_psqsums[ej - 1]);
      } else if (i == 0 && j < m) {
        mean_v2 = (obvs_psums[em - 1] - obvs_psums[ej - 1]) / m_v2;
        var_v2 = Variance(m_v2, mean_v2, obvs_psqsums[em - 1] - obvs_psqsums[ej - 1]);
      } else if (i > 0 && j == m) {
        mean_v2 = (obvs_psums[ei - 1]) / m_v2;
        var_v2 = Variance(m_v2, mean_v2, obvs_psqsums[ei - 1]);
      } else {
        throw std::runtime_error("TStat: i==0 and j==m; cannot make segment");
      }
      const f64 t = TStat(m_v1, mean_v1, var_v1, m_v2, mean_v2, var_v2);
      if (std::abs(t) > max_t) {
        max_t = std::abs(t);
        max_t_i = i;
        max_t_j = j;
        max_t_mean = mean_v1;
        max_t_m = m_v1;
      }
    }
  }
  Segment ret;

  ret.start = max_t_i;
  ret.end = max_t_j;
  ret.t = max_t;
  ret.mean = max_t_mean;
  ret.n = max_t_m;
  return ret;
}

/**
 * @brief permute null distribution describing scenario where the maximal
 * t-value of a circular parition of the observation std::vector is not significantly
 * deviant from the maximal t-value of a circular partitioning of any random
 * permutation of the observations
 * acccept null if (max_p*n_permutations) permutations' max-t statistics exceed t
 * reject null if above condition fails, or if boundaries are provided and the boundary test fails at any permutation
 * @param obvs observation std::vector of doubles
 * @param n_permutations max number of permutations to test
 * @param max_p maximum p-value for rejecting null
 * @param t maximal t-statistic for circular partitioning of observation std::vector
 * @param max_obvs_size maximum number of observations allowed per segment
 * @param boundaries pre-calculated boundaries for early-stopping
 * @param min_obs_size minimum number of observations allowed per segment
 * @return true if null is rejected else false
 */
bool CircularPermuteNullAndTestSignificance(const Eigen::VectorXd& obvs,
                                            size_t n_permutations,
                                            const f64 max_p,
                                            const f64 t,
                                            const size_t max_obs_size,
                                            const vec<u64>& boundaries,
                                            const size_t min_obs_size,
                                            eigen::ArmaRng& rng) {
  const auto max_n_gt_t = static_cast<size_t>(max_p * static_cast<f64>(n_permutations));
  Logging::Debug("Estimating null over maximum {} iters using early stopping measure. hard stop at {} excds",
                 n_permutations,
                 max_n_gt_t + 1);
  size_t n_gt_t = 0;
  // k denotes start-index in boundaries std::vector
  size_t k = max_n_gt_t * (max_n_gt_t + 1) / 2;
  Logging::Debug("boundary (std::vector size {}) start index is: {}", boundaries.size(), k);
  for (size_t i = 0; i < n_permutations; ++i) {  // NOSONAR(cpp:S5313) loop counter is not modified in body
    rng.SetSeed(i);
    const Eigen::VectorXd obvs_perm = rng.Shuffle(obvs);
    const Segment perm_max_t_seg = GetMaxTSegCumSum(obvs_perm, max_obs_size, min_obs_size);
    if (t < perm_max_t_seg.t) {
      n_gt_t += 1;
    }
    // if the number of non-reject permutations exceeds <max_p*n_permutations>, then we can automatically say that the
    // segment isn't significant
    if (n_gt_t > max_n_gt_t) {
      Logging::Debug(
          "stopped permuting for segment of size {} at {} iters {} excds. no reject null", obvs.size(), i, n_gt_t);
      return false;
    }
    if (!boundaries.empty() && boundaries[k + n_gt_t] < i) {
      Logging::Debug("stopped permuting for segment of size {} at {} iters. reject null", obvs.size(), i);
      return true;
    }
  }
  Logging::Debug("permuted maximum amount of times for segment of size {} at {} iters", obvs.size(), n_permutations);
  return static_cast<f64>(n_gt_t) / static_cast<f64>(n_permutations) < max_p;
}

/**
 * @brief permute null distribution describing scenario where the t-value of a
 * binary parition of the observation std::vector is not significantly deviant from
 * the maximal t-value of a proportionally similar binary partitioning of any
 * random permutation of the observations
 * acccept null if (max_p*n_permutations) permutations' max-t statistics exceed t
 * reject null if above condition fails
 * @param obvs observation std::vector of doubles
 * @param n_permutations max number of permutations to test
 * @param max_p maximum p-value for rejecting null
 * @param seg1 first segment in the binary partition
 * @param seg2 second segment in the binary partition
 * @return true if null is rejected else false
 */
bool BinaryPermuteNullAndTestSignificance(const Eigen::VectorXd& obvs,
                                          const size_t n_permutations,
                                          const f64 max_p,
                                          const Segment& seg1,
                                          const Segment& seg2,
                                          eigen::ArmaRng& rng) {
  Eigen::VectorXd v1 =
      obvs.segment(static_cast<Eigen::Index>(seg1.start), static_cast<Eigen::Index>(seg1.end - seg1.start));
  Eigen::VectorXd v2 =
      obvs.segment(static_cast<Eigen::Index>(seg2.start), static_cast<Eigen::Index>(seg2.end - seg2.start));
  const Eigen::VectorXd to_cmp = eigen::JoinCols(v1, v2);
  const auto n_v1 = static_cast<u64>(v1.size());
  f64 nexcd = 0.0;
  const f64 max_excd = (static_cast<f64>(n_permutations) * max_p) + 1;
  const f64 t = std::abs(TStatFromVecs(v1, v2));
  for (size_t i = 0; i < n_permutations; ++i) {
    const Eigen::VectorXd perm = rng.Shuffle(to_cmp);
    v1 = perm.segment(0, static_cast<Eigen::Index>(n_v1));
    v2 = perm.segment(static_cast<Eigen::Index>(n_v1), perm.size() - static_cast<Eigen::Index>(n_v1));
    const f64 perm_t = std::abs(TStatFromVecs(v1, v2));
    if (perm_t > t) {
      nexcd += 1;
    }
    if (nexcd >= max_excd) {
      Logging::Debug("adjacent segs of size {} / mean {} and size {} / mean {} are not significantly different, p={}",
                     seg1.end - seg1.start,
                     eigen::Mean(obvs.segment(static_cast<Eigen::Index>(seg1.start),
                                              static_cast<Eigen::Index>(seg1.end - seg1.start))),
                     seg2.end - seg2.start,
                     eigen::Mean(obvs.segment(static_cast<Eigen::Index>(seg2.start),
                                              static_cast<Eigen::Index>(seg2.end - seg2.start))),
                     nexcd / static_cast<f64>(n_permutations));
      return false;
    }
  }
  return nexcd / static_cast<f64>(n_permutations) < max_p;
}

const f64 kCbsSiegmundPMinChange = 1e6;
const f64 kCbsSiegmundPMaxIter = 100'000;

/**
 * @brief described in Siegmund, "Approximate Tail Probabilities for the Maxima of Some Random Fields", 1988
 * @param x the parameter for Siegmund apprxomiation Nu function
 * @return the "nu" parameter for the Siegment p-value approximation equation
 *
 */
f64 SiegmundPNuFull(const f64 x) {
  auto fn = [&x](const f64 l) {
    const f64 power = std::pow(l, -1);
    const f64 half_x = -0.5 * x;
    const f64 sqrt_l = std::pow(l, 0.5);
    const Eigen::VectorXd vec = Eigen::VectorXd::Constant(1, half_x * sqrt_l);
    return power * eigen::Normcdf(vec)[0];
  };
  auto change_frac = [](const f64 p, const f64 n) {
    return math::IsCloseToZero(p) ? std::numeric_limits<f64>::max() : std::abs(p - n) / p;
  };
  // need to initialize ptotal and total to avoid div by 0
  f64 ptotal = fn(1);
  f64 total = ptotal + fn(2);
  f64 l = 3;
  while (change_frac(ptotal, total) > kCbsSiegmundPMinChange && l < kCbsSiegmundPMaxIter) {
    ptotal = total;
    total += fn(l);
    l += 1;
  }
  f64 ret = 2 * std::pow(x, -2);
  ret *= std::exp(-2 * total);
  return ret;
}

/**
 * @brief P-value approximation described in Siegmund, "Approximate Tail
 * Probabilities for the Maxima of Some Random Fields", 1986. This is used for
 * when you want a p-value from comparing the maxima of two samples. For
 * intuition, if you're comparing the means instead of the maxima then you can
 * get away with just using the normal CDF.
 * @param t test statistic for which to calculate significance
 * @param m size of sample
 * @param k scale parameter
 */
f64 SiegmundP(const f64 t, const f64 m, const f64 k) {
  const f64 pdf_t = eigen::Normpdf(Eigen::VectorXd::Constant(1, t))[0];
  f64 p = kSiegmundLeadingCoeff * std::pow(t, 3) * pdf_t;
  const f64 delta = k / m;
  const f64 incr = (1 - delta - 0.5) / static_cast<f64>(kCbsSiegmundPIters);
  Eigen::VectorXd fn_to_intg_x(kCbsSiegmundPIters);
  Eigen::VectorXd fn_to_intg_y(kCbsSiegmundPIters);
  for (size_t i = 0; i < kCbsSiegmundPIters; ++i) {
    const f64 k = 0.5 + (static_cast<f64>(i) * incr);
    const f64 k_1_min_k = k * (1 - k);
    const f64 nu_sq = std::pow(SiegmundPNuFull(t / std::pow(m * k_1_min_k, 0.5)), 2);
    const f64 y = nu_sq / std::pow(k_1_min_k, 2);
    fn_to_intg_x[static_cast<Eigen::Index>(i)] = k;
    fn_to_intg_y[static_cast<Eigen::Index>(i)] = y;
  }
  const f64 integral_area = eigen::Trapz(fn_to_intg_x, fn_to_intg_y);
  p = integral_area * p;
  return p;
}

/**
 * @brief the "hybrid" test for significance of a circular segmentation. If the
 * two segments are large enough, then use Siegmund's approximation to calculate
 * the portion of the p-value associated with >k sized segments. If the
 * partial-p is significant, use permutation testing to calculate the portion of
 * the p-value associated with <k sized segments and make a final decision on
 * significance. If the one of the segments is too small, default to permutation
 * testing for the whole p-value.
 * @param obvs observation std::vector
 * @param n_permutations max number of permutations to test
 * @param max_p maximum p-value for rejecting null
 * @param t maximal t-statistic for circular partitioning of observation std::vector
 * @param boundaries pre-calculated boundaries for early-stopping
 * @param min_obs_size minimum number of observations allowed per segment
 * @return true if the segmentation was significant according to the p-value calculation, false if not
 */
bool HybridTestForSignificance(const Eigen::VectorXd& obvs,
                               const size_t n_permutations,
                               const f64 max_p,
                               const f64 t,
                               const vec<u64>& boundaries,
                               const size_t min_obs_size,
                               eigen::ArmaRng& rng) {
  if (static_cast<u64>(obvs.size()) < kCbsHybridK * kCbsHybridMinSizeMultiple) {
    Logging::Debug("{} elem segment: using default P method", obvs.size());
    return CircularPermuteNullAndTestSignificance(obvs, n_permutations, max_p, t, 0, boundaries, min_obs_size, rng);
  }
  Logging::Debug("{} elem segment: using approx P method", obvs.size());
  // get null distribution for all splits minor_arc>=k && major_arc>=k
  f64 p = SiegmundP(t, static_cast<f64>(obvs.size()), static_cast<f64>(kCbsHybridK));
  if (p < max_p) {
    // get null distribution for all splits minor_arc<k to verify that p actually is significant
    Logging::Debug("{} elem segment: approx p is significant ({}), will also run P method", obvs.size(), p);
    return CircularPermuteNullAndTestSignificance(
        obvs, n_permutations, max_p - p, t, kCbsHybridK, boundaries, min_obs_size, rng);
  }
  return false;
}

/**
 * @brief Find the best scoring sub-segment (by t-statistic) and, if
 * statistically significant, return all segments created from this best-scoring
 * sub-segment's breakpoints. If not significant, return an empty std::vector
 * @param obvs observation std::vector
 * @param max_p maximum p-value for rejecting null
 * @param n_permutations max number of permutations to test
 * @param boundaries pre-calculated boundaries for early-stopping
 * @param min_obs_size minimum number of observations allowed per segment
 * @return a maximum of three segments w.r.t to obvs. These segments should span the entirety of obvs
 */
vec<Segment> GetSegments(const Eigen::VectorXd& obvs,
                         const f64 max_p,
                         const size_t n_permutations,
                         const vec<u64>& boundaries,
                         const size_t min_obs_size,
                         eigen::ArmaRng& rng) {
  // get the highest t-stat and corresponding segment for the real observations
  Logging::Debug("finding max-T sub-segment for {} obvs segment", obvs.size());
  const Segment max_t_seg = GetMaxTSegCumSum(obvs, 0, min_obs_size);
  Logging::Debug("found max-T sub-segment for {} obvs segment", obvs.size());
  if (math::IsCloseToZero(max_t_seg.t)) {
    Logging::Debug("segment with {} obvs is completely uniform", obvs.size());
    return {};
  }
  if (max_t_seg.end - max_t_seg.start < kCbsMinObvsN) {
    Logging::Debug("segment with {} obvs has very small max-t sub-segment", obvs.size());
    return {};
  }
  const bool is_significant =
      HybridTestForSignificance(obvs, n_permutations, max_p, max_t_seg.t, boundaries, min_obs_size, rng);
  Logging::Debug("found significance for max-T sub-segment for {} obvs segment", obvs.size());
  if (!is_significant) {
    return vec<Segment>{};
  }
  vec<Segment> out_segments;
  const auto make_segment = [&obvs](const size_t start, const size_t end) {
    const auto n = end - start;
    return Segment{.start = start,
                   .end = end,
                   .n = static_cast<f64>(n),
                   .mean = eigen::Mean(obvs.segment(static_cast<Eigen::Index>(start), static_cast<Eigen::Index>(n)))};
  };
  if (max_t_seg.start != 0) {
    out_segments.push_back(make_segment(0, max_t_seg.start));
  }
  out_segments.push_back(make_segment(max_t_seg.start, max_t_seg.end));
  if (max_t_seg.end != static_cast<u64>(obvs.size())) {
    out_segments.push_back(make_segment(max_t_seg.end, static_cast<size_t>(obvs.size())));
  }
  return out_segments;
}

/**
 * @brief Segment an observation std::vector using the Circular Binary Segmentation
 * Algorithm described in Olshen et al., 2004, "Circular binary segmentation
 * for the analysis of array- based DNA copy number data", 2004.
 * and implementing the speed-ups described in Venkatraman and Olshen, 2007,
 * "A faster circular binary segmentation algorithm for the analysis of array
 * CGH data", 2007.
 * @param obvs observation std::vector
 * @param max_p maximum p-value for rejecting null
 * @param n_permutations max number of permutations to test
 * @param min_obs_size minimum number of observations allowed per segment
 * @return all significant segments found in obvs
 */
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
vec<Segment> CircularBinarySegmentation(const Eigen::VectorXd& obvs,
                                        const f64 max_p,
                                        const size_t n_permutations,
                                        const size_t min_obs_size,
                                        const bool no_single_obvs_subsegments) {
  auto time1 = high_resolution_clock::now();
  if (obvs.size() == 0) {
    return {};
  }
  // One RNG object per CBS run: each taskflow worker runs a single CBS run, so state carried between permute calls
  // (e.g. BinaryPermute without reseed) must persist within the run but not across runs.
  eigen::ArmaRng rng;
  vec<Segment> ret;
  std::stack<Segment> seg_stack;
  const Segment first_seg = {.start = 0, .end = static_cast<u64>(obvs.size())};
  seg_stack.push(first_seg);
  while (!seg_stack.empty()) {
    Logging::Debug("stack has {} elements", seg_stack.size());
    const Segment& segment = seg_stack.top();
    seg_stack.pop();
    // don't process a segment if it's too small
    if (segment.end - segment.start < kCbsMinObvsN) {
      Logging::Debug("adding small segment with {} obvs", (segment.end - segment.start));
      ret.push_back(segment);
    } else {
      const Eigen::VectorXd sub_obvs = obvs.segment(static_cast<Eigen::Index>(segment.start),
                                                    static_cast<Eigen::Index>(segment.end - segment.start));
      // "recursive" step - partition segment into sub-segments
      vec<u64> boundaries;
      if (n_permutations == kPrecomputedBoundaryPermutations &&
          math::IsEqualWithTolerance(max_p, kPrecomputedBoundaryMaxP)) {
        boundaries = kBoundaryP10000A01Eta05;
      }
      Logging::Debug("Segmenting on {} obs", sub_obvs.size());
      vec<Segment> sub_segments = GetSegments(sub_obvs, max_p, n_permutations, boundaries, min_obs_size, rng);
      Logging::Debug("Segmented on {} obs", sub_obvs.size());
      // merge back single-obs sub-segments
      if (no_single_obvs_subsegments) {
        sub_segments = MergeSingleObvsEndSegment(sub_segments, min_obs_size);
      }
      if (sub_segments.size() <= 1) {
        // if this segment cannot be broken down, then return it as one segment
        Logging::Debug("adding un-breakable segment with {} obvs", (segment.end - segment.start));
        ret.push_back(segment);
      } else {
        // "verify" that the segments are actually different from each other by doing the binary segmentation t-test
        // procedure .
        // First push all segments onto a stack in reverse order, such that the first segment is on the top of the
        // stack, etc.
        std::stack<Segment> segment_merge_stack;
        for (const auto it : sub_segments | std::views::reverse) {
          segment_merge_stack.push(it);
        }
        sub_segments.clear();
        // while there are at least 2 elements in the stack, compare the top two
        // Segments in the stack. If they fail the binary segmentation test of
        // significance, then merge these two segments, and add the merged
        // segment onto the stack. If they are significantly different, then add
        // the top Segment to the array, and keep the second Segment in the
        // stack, so that it can be tested against the next Segment in the
        // following iteration
        // The loop necessarily stops when there is only 1 item in the stack.
        // This item is then added to the return std::vector
        while (segment_merge_stack.size() > 1) {
          const Segment seg_to_merge1 = segment_merge_stack.top();
          segment_merge_stack.pop();
          const Segment seg_to_merge2 = segment_merge_stack.top();
          if (!BinaryPermuteNullAndTestSignificance(
                  sub_obvs, n_permutations, max_p, seg_to_merge1, seg_to_merge2, rng)) {
            const Segment new_sub_segment(seg_to_merge1.start, seg_to_merge2.end);
            // remove seg_to_merge2 from top
            segment_merge_stack.pop();
            segment_merge_stack.push(new_sub_segment);
          } else {
            sub_segments.push_back(seg_to_merge1);
            // seg_to_merge2 is kept at the top of the stack
          }
        }
        sub_segments.push_back(segment_merge_stack.top());
        // case for when *all* the subsegments ended up being merged back together
        if (sub_segments.size() == 1) {
          Logging::Debug("adding un-breakable segment with {} obvs", (segment.end - segment.start));
          ret.push_back(segment);
        } else {
          // these sub-segments make no assumption about which "parent" its observations come from. Therefore, we have
          // to re-adjust the start/end fields of the sub-segments with respect to the parent std::vector's coordinate
          // system
          for (auto new_segment : sub_segments) {
            new_segment.start += segment.start;
            new_segment.end += segment.start;
            if (new_segment.end > static_cast<u64>(obvs.size())) {
              throw error::Error(
                  "cbs-for-somatic: sub-segment end exceeds observation count after coordinate adjustment");
            }
            Logging::Debug("pushing segment with {} obvs", (new_segment.end - new_segment.start));
            seg_stack.push(new_segment);
            if (seg_stack.size() > kCbsMaxStackDepth) {
              throw std::runtime_error("cbs: max depth exceeded!");
            }
          }
        }
      }
    }
  }
  std::ranges::sort(ret, [](const Segment& l, const Segment& r) { return l.start < r.start; });
  const auto time2 = high_resolution_clock::now();
  auto total_time = duration_cast<milliseconds>(time2 - time1).count();
  Logging::Debug("found {} segments over {} elements in {}ms", ret.size(), obvs.size(), total_time);
  return ret;
}

/**
 * @brief Given a set of 2-3 subsegments, if the first or the last sub-segment has fewer than min_obs_size observations,
 * then merge it with its adjacent sub-segment
 * @param sub_segments a std::vector of 2-3 sub-segments
 * @param min_obs_size  minimum number of observations allowed per segment
 * @return
 */
vec<Segment> MergeSingleObvsEndSegment(const vec<Segment>& sub_segments, size_t min_obs_size) {
  if (sub_segments.size() <= 1) {
    return sub_segments;
  }
  vec<Segment> res;
  if (sub_segments.size() < 4) {
    auto first_seg = sub_segments.begin();
    auto second_seg = sub_segments.begin() + 1;
    if (sub_segments.begin()->end - sub_segments.begin()->start < min_obs_size) {
      // first segment is small
      first_seg = sub_segments.begin();
      second_seg = sub_segments.begin() + 1;
    } else if ((sub_segments.end() - 1)->end - (sub_segments.end() - 1)->start < min_obs_size) {
      // last segment is small
      first_seg = sub_segments.end() - 2;
      second_seg = sub_segments.end() - 1;
    } else {
      // neither edge is "small", can return as-is
      return sub_segments;
    }
    const f64 merged_n = first_seg->n + second_seg->n;
    const Segment merged_segment = {
        .start = first_seg->start,
        .end = second_seg->end,
        .t = second_seg->t,
        .p = second_seg->p,
        .n = merged_n,
        .mean = ((first_seg->n * first_seg->mean) + (second_seg->n * second_seg->mean)) / merged_n,
    };
    if (sub_segments.size() == 2) {
      res.resize(1);
      res[0] = merged_segment;
      return res;
    }  // size == 3
    res.resize(2);
    if (sub_segments.begin()->end - sub_segments.begin()->start < min_obs_size) {
      // first segment is small
      res[0] = merged_segment;
      res[1] = sub_segments[2];
    } else if ((sub_segments.end() - 1)->end - (sub_segments.end() - 1)->start < min_obs_size) {
      // last segment is small
      res[0] = sub_segments[0];
      res[1] = merged_segment;
    }
    return res;
  }
  Logging::Error("Will not prune segment with more than three sub-segments");
  throw std::runtime_error("Will not prune segment with more than three sub-segments");
}
}  // namespace xoos::cnc::somatic
