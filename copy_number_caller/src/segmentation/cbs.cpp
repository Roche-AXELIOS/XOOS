/**
 * reference:
 * Olshen et al, 2004. "Circular binary segmentation for the analysis of array‐based DNA copy number data"
 * https://academic.oup.com/biostatistics/article/5/4/557/275197
 */

#include "segmentation/cbs.h"

#include <Eigen/Dense>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <ranges>
#include <stack>
#include <stdexcept>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/stats/external/asa243/asa243.hpp>
#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/util/math.h>

#include "segmentation/boundaries.h"
#include "segmentation/cbs-util.h"
#include "segmentation/genomic-segments.h"
#include "segmentation/max-t.h"
#include "segmentation/segmentation-options.h"
#include "util/eigen-math.h"

namespace xoos::cnc {

using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;
using std::chrono::nanoseconds;
using std::chrono::seconds;

const size_t kCbsSiegmundPIters = 100;
const size_t kCbsMinObvsForRecursion = 10;
const size_t kCbsHybridK = 25;
constexpr f64 kSiegmundLeadingCoeff = 0.25;
constexpr f64 kTrimLowerQuantile = 0.01;
constexpr f64 kTrimUpperQuantile = 0.99;
const size_t kCbsMinObvsForShortcut = 10;
const f64 kCbsMaxTForAutomaticNonRejectNull = 0.1;
const size_t kCbsMaxStackDepth = 1024;
const size_t kCbsDefaultNPermutations = 10'000;
const f64 kCbsDefaultMaxP = 0.01;
// degrees-of-freedom offset for the two-partition t-statistic denominator (n - 2)
const f64 kTStatDofOffset = 2;
// minimum segment size (as a multiple of kCbsHybridK) required to use the approximate p-value method
const size_t kCbsHybridMinSizeMultiple = 4;

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
 * @param min_obs_per_segment minimum number of observations allowed per segment
 * @return true if null is rejected else false
 */
bool CircularPermuteNullAndTestSignificance(const Eigen::VectorXd& obvs,
                                            size_t n_permutations,
                                            const f64 max_p,
                                            const Segment& max_t_seg,
                                            const size_t max_obs_size,
                                            const vec<u64>& boundaries,
                                            const size_t min_obs_per_segment,
                                            eigen::ArmaRng& rng,
                                            const CbsMaxTMethod max_t_method) {
  const auto max_n_gt_t = static_cast<size_t>(max_p * static_cast<f64>(n_permutations));
  Logging::Debug("Estimating null over maximum {} iters using early stopping measure. hard stop at {} excds",
                 n_permutations,
                 max_n_gt_t + 1);
  size_t n_gt_t = 0;
  // k denotes start-index in boundaries std::vector
  size_t k = max_n_gt_t * (max_n_gt_t + 1) / 2;
  Eigen::VectorXd shuffled_obvs(obvs);
  const auto obvs_size = static_cast<size_t>(obvs.size());
  Logging::Debug("boundary (std::vector size {}) start index is: {}", boundaries.size(), k);
  const Eigen::VectorXd centered_obvs = obvs.array() - eigen::Mean(obvs);
  const f64 tss = eigen::Sum(centered_obvs.array().square().matrix());
  for (size_t i = 0; i < n_permutations; ++i) {  // NOSONAR(cpp:S5313) loop counter is not modified in body
    rng.SetSeed(i);
    FisherYatesShuffleInPlace(shuffled_obvs, rng);
    Segment perm_max_t_seg = GetMaxTBreakpoints(shuffled_obvs, max_obs_size, min_obs_per_segment, max_t_method);
    const f64 new_tss = tss < perm_max_t_seg.t ? perm_max_t_seg.t + 1 : tss;
    perm_max_t_seg.t =
        perm_max_t_seg.t / ((new_tss - perm_max_t_seg.t) / (static_cast<f64>(obvs_size) - kTStatDofOffset));
    if (max_t_seg.t < perm_max_t_seg.t) {
      n_gt_t += 1;
    }
    // if the number of non-reject permutations exceeds <max_p*n_permutations>, then we can automatically say that the
    // segment isn't significant
    if (n_gt_t > max_n_gt_t) {
      Logging::Debug(
          "stopped permuting for segment of size {} at {} iters {} excds. no reject null", obvs_size, i, n_gt_t);
      return false;
    }
    if (!boundaries.empty() && boundaries[k + n_gt_t] < i) {
      Logging::Debug("stopped permuting for segment of size {} at {} iters. reject null", obvs_size, i);
      return true;
    }
  }
  Logging::Debug("permuted maximum amount of times for segment of size {} at {} iters", obvs_size, n_permutations);
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
                                          const f64 min_t_for_automatic_reject_null,
                                          eigen::ArmaRng& rng) {
  const Eigen::VectorXd v1 =
      obvs.segment(static_cast<Eigen::Index>(seg1.start), static_cast<Eigen::Index>(seg1.end - seg1.start));
  const Eigen::VectorXd v2 =
      obvs.segment(static_cast<Eigen::Index>(seg2.start), static_cast<Eigen::Index>(seg2.end - seg2.start));
  const f64 t = std::abs(TStatFromVecs(v1, v2));
  // shortcut: reject null if both segments are large enough and t is large
  if (t > min_t_for_automatic_reject_null && static_cast<size_t>(v1.size()) > kCbsMinObvsForShortcut &&
      static_cast<size_t>(v2.size()) > kCbsMinObvsForShortcut) {
    return true;
  }
  // otherwise if both segments are large enough, we can use the student's t distribution to assess if both segments are
  // signiciantly different from each other
  if (static_cast<size_t>(v1.size()) > kCbsHybridK && static_cast<size_t>(v2.size()) > kCbsHybridK) {
    const auto df = static_cast<f64>(v1.size() + v2.size() - 2);
    const f64 p = 2 * stats::tnc(-std::abs(t), df, 0);
    return p < max_p;
    // if one or both segments are very small, then we have to use permutation testing because it's not enough
    // observations to assume that the data follows a normal distribution
  }
  Eigen::VectorXd to_cmp = eigen::JoinCols(v1, v2);
  const f64 total_sum = eigen::Sum(to_cmp);
  const f64 total_sum_of_squares = eigen::Sum(to_cmp.array().square().matrix());
  auto minor_arc_len = static_cast<f64>(std::min(v1.size(), v2.size()));
  auto major_arc_len = static_cast<f64>(std::max(v1.size(), v2.size()));
  f64 nexcd = 0.0;
  const f64 max_excd = (static_cast<f64>(n_permutations) * max_p) + 1;
  for (size_t i = 0; i < n_permutations; ++i) {
    // since we know the minor arc is small, we don't need to permute the
    // entire array. We can just sample minor_arc_len points, calculate sums for those directly and infer the sums for
    // the major arc from the minor arc sum
    rng.SetSeed(i);
    const eigen::VectorXu64 minor_arc_sample_idxs =
        rng.RandiFill(static_cast<u64>(minor_arc_len), 0, static_cast<s32>(to_cmp.size() - 1));
    f64 minor_arc_sum = 0;
    f64 minor_arc_sum_of_squares = 0;
    for (Eigen::Index s = 0; s < minor_arc_sample_idxs.size(); ++s) {
      const auto idx = static_cast<Eigen::Index>(minor_arc_sample_idxs(s));
      minor_arc_sum += to_cmp(idx);
      minor_arc_sum_of_squares += std::pow(to_cmp(idx), 2);
    }
    const f64 major_arc_sum = total_sum - minor_arc_sum;
    const f64 major_arc_sum_of_squares = total_sum_of_squares - minor_arc_sum_of_squares;
    // calculate the t statistic
    const f64 minor_arc_mean = minor_arc_sum / minor_arc_len;
    const f64 minor_arc_var = Variance(minor_arc_len, minor_arc_mean, minor_arc_sum_of_squares);
    const f64 major_arc_mean = major_arc_sum / major_arc_len;
    const f64 major_arc_var = Variance(major_arc_len, major_arc_mean, major_arc_sum_of_squares);
    const f64 perm_t =
        std::abs(TStat(minor_arc_len, minor_arc_mean, minor_arc_var, major_arc_len, major_arc_mean, major_arc_var));
    if (perm_t > t) {
      nexcd += 1;
    }
    if (nexcd >= max_excd) {
      Logging::Debug("adjacent segs of size {} / mean {} and size {} / mean {} are not significantly different, p={}",
                     seg1.end - seg1.start,
                     eigen::Mean(v1),
                     seg2.end - seg2.start,
                     eigen::Mean(v2),
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
    Eigen::VectorXd arg(1);
    arg << half_x * sqrt_l;
    return power * eigen::Normcdf(arg)[0];
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
  Eigen::VectorXd t_vec(1);
  t_vec << t;
  const f64 pdf_t = eigen::Normpdf(t_vec)[0];
  f64 p = kSiegmundLeadingCoeff * std::pow(t, 3) * pdf_t;
  const f64 delta = k / m;
  const f64 incr = (1 - delta - 0.5) / static_cast<f64>(kCbsSiegmundPIters);
  Eigen::VectorXd fn_to_intg_x(static_cast<Eigen::Index>(kCbsSiegmundPIters));
  Eigen::VectorXd fn_to_intg_y(static_cast<Eigen::Index>(kCbsSiegmundPIters));
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
 * @param min_obs_per_segment minimum number of observations allowed per segment
 * @return true if the segmentation was significant according to the p-value calculation, false if not
 */
bool HybridTestForSignificance(const Eigen::VectorXd& obvs,
                               const size_t n_permutations,
                               const f64 max_p,
                               Segment max_t_seg,
                               const vec<u64>& boundaries,
                               const size_t min_obs_per_segment,
                               eigen::ArmaRng& rng,
                               const CbsMaxTMethod max_t_method,
                               f64 min_t_for_automatic_reject_null) {
  const auto obvs_size = static_cast<size_t>(obvs.size());
  // if the two arcs are large enough, then we can reject null if t is a large-ish number
  if (max_t_seg.end - max_t_seg.start > kCbsMinObvsForShortcut &&
      obvs_size - max_t_seg.end + max_t_seg.start > kCbsMinObvsForShortcut &&
      max_t_seg.t > min_t_for_automatic_reject_null) {
    Logging::Debug("{} elem segment because t ({}) > {}", obvs_size, max_t_seg.t, min_t_for_automatic_reject_null);
    return true;
    // we can also not-reject null automaticaly if t is realy small
  }
  if (max_t_seg.t < kCbsMaxTForAutomaticNonRejectNull) {
    Logging::Debug("not segmenting {} elem segment because t ({}) < {}",
                   obvs_size,
                   max_t_seg.t,
                   kCbsMaxTForAutomaticNonRejectNull);
    return false;
  }
  // if there aren't enough observations, we cannnot calculate approximate p. We have to do permutation testing on the
  // whole thing
  if (obvs_size < kCbsHybridK * kCbsHybridMinSizeMultiple) {
    Logging::Debug("{} elem segment: using default P method", obvs_size);
    return CircularPermuteNullAndTestSignificance(
        obvs, n_permutations, max_p, max_t_seg, 0, boundaries, min_obs_per_segment, rng, max_t_method);
  }
  // if there enough observations, we can approximate p for all the partitions where both arcs are big enough
  Logging::Debug("{} elem segment: using approx P method", obvs_size);
  // get null distribution for all splits minor_arc>=k && major_arc>=k
  f64 p = SiegmundP(max_t_seg.t, static_cast<f64>(obvs_size), static_cast<f64>(kCbsHybridK));
  // if this is significant, then we calculate the part of the p-value for all the paritions where one arc is small
  if (p < max_p) {
    Logging::Debug("{} elem segment: approx p is significant ({}), will also run P method", obvs_size, p);
    if (max_t_method == CbsMaxTMethod::kFast) {
      return CircularPermuteNullAndTestSignificance(obvs,
                                                    n_permutations,
                                                    max_p - p,
                                                    max_t_seg,
                                                    kCbsHybridK,
                                                    boundaries,
                                                    min_obs_per_segment,
                                                    rng,
                                                    CbsMaxTMethod::kFastMaxN);
    }
    return CircularPermuteNullAndTestSignificance(obvs,
                                                  n_permutations,
                                                  max_p - p,
                                                  max_t_seg,
                                                  kCbsHybridK,
                                                  boundaries,
                                                  min_obs_per_segment,
                                                  rng,
                                                  CbsMaxTMethod::kBruteForce);
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
 * @param min_obs_per_segment minimum number of observations allowed per segment
 * @return a maximum of three segments w.r.t to obvs. These segments should span the entirety of obvs
 */
vec<Segment> GetSegments(const Eigen::VectorXd& obvs,
                         f64 max_p,
                         size_t n_permutations,
                         const vec<u64>& boundaries,
                         size_t min_obs_per_segment,
                         eigen::ArmaRng& rng,
                         CbsMaxTMethod max_t_method,
                         f64 min_t_for_automatic_reject_null) {
  const auto obvs_size = static_cast<size_t>(obvs.size());
  // get the highest t-stat and corresponding segment for the real observations
  Logging::Debug("finding max-T sub-segment for {} obvs segment", obvs_size);
  const Eigen::VectorXd centered_obvs = obvs.array() - eigen::Mean(obvs);
  f64 tss = eigen::Sum(centered_obvs.array().square().matrix());
  Segment max_t_seg = GetMaxTBreakpoints(centered_obvs, 0, min_obs_per_segment, max_t_method);
  tss = tss < max_t_seg.t ? max_t_seg.t + 1 : tss;
  max_t_seg.t = max_t_seg.t / ((tss - max_t_seg.t) / (static_cast<f64>(obvs_size) - kTStatDofOffset));
  Logging::Debug("found max-T sub-segment for {} obvs segment", obvs_size);
  // conditions for immediate stop on recursion (no change, segmentation too small)
  if (math::IsCloseToZero(max_t_seg.t)) {
    Logging::Debug("segment with {} obvs is completely uniform", obvs_size);
    return {};
  }
  if (max_t_seg.end - max_t_seg.start < min_obs_per_segment) {
    Logging::Debug("segment with {} obvs has very small max-t sub-segment", obvs_size);
    return {};
  }
  if (!HybridTestForSignificance(obvs,
                                 n_permutations,
                                 max_p,
                                 max_t_seg,
                                 boundaries,
                                 min_obs_per_segment,
                                 rng,
                                 max_t_method,
                                 min_t_for_automatic_reject_null)) {
    return vec<Segment>{};
  }
  Logging::Debug("found significance for max-T sub-segment for {} obvs segment", obvs_size);
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
  if (max_t_seg.end != obvs_size) {
    out_segments.push_back(make_segment(max_t_seg.end, obvs_size));
  }
  return out_segments;
}

Eigen::VectorXd TruncateOutliers(const Eigen::VectorXd& obvs) {
  Eigen::VectorXd ret(obvs);
  // get the 0.1 and .99
  Eigen::VectorXd probs(2);
  probs(0) = kTrimLowerQuantile;
  probs(1) = kTrimUpperQuantile;
  const Eigen::VectorXd quantiles = eigen::Quantile(obvs, probs);
  const f64 interquantile_range = quantiles[1] - quantiles[0];
  const f64 lower_bound = quantiles[0] - (1.5 * interquantile_range);
  const f64 upper_bound = quantiles[1] + (1.5 * interquantile_range);
  for (Eigen::Index idx = 0; idx < ret.size(); ++idx) {
    // truncate if the observation is above or below the quantiles
    ret(idx) = std::max(ret(idx), lower_bound);
    ret(idx) = std::min(ret(idx), upper_bound);
  }
  return ret;
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
 * @param min_obs_per_segment minimum number of observations allowed per segment
 * @return all significant segments found in obvs
 */
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
vec<Segment> CircularBinarySegmentation(const Eigen::VectorXd& obvs, const SegmentationOptions& options) {
  auto n = static_cast<size_t>(obvs.size());
  // One RNG per CBS run: each SegmentByOneSegment call runs on a single taskflow
  // worker and performs exactly one CBS run. The germline permutation tests
  // reseed every iteration, so this object's carried state never crosses
  // iterations here, but threading a single instance keeps the model identical
  // to the somatic path.
  eigen::ArmaRng rng;
  // we need to truncate outliers to the same max value because one or two very extreme outliers can cause significance
  // testing to fail, causing severe under-segmentation
  Eigen::VectorXd use_these_obvs(obvs);
  if (options.truncate_outliers) {
    Logging::Info("truncating outliers");
    use_these_obvs = TruncateOutliers(obvs);
  }
  const auto time1 = high_resolution_clock::now();
  if (n == 0) {
    Logging::Info("no observations to segment!");
    return {};
  }
  vec<Segment> ret;
  std::stack<Segment> seg_stack;
  const Segment first_seg = {.start = 0, .end = n};
  seg_stack.push(first_seg);
  while (!seg_stack.empty()) {
    Logging::Debug("stack has {} elements", seg_stack.size());
    const Segment& segment = seg_stack.top();
    seg_stack.pop();
    // If segment is small enough that we can't perform permutation testing if we break it up, then don't recurse
    if (segment.end - segment.start < kCbsMinObvsForRecursion) {
      Logging::Debug("adding small segment with {} obvs", (segment.end - segment.start));
      ret.push_back(segment);
    } else {
      // NOTE: we can probably make this faster if we pass the indices into GetSegments instead of copying
      const Eigen::VectorXd sub_obvs = use_these_obvs.segment(static_cast<Eigen::Index>(segment.start),
                                                              static_cast<Eigen::Index>(segment.end - segment.start));
      // "recursive" step - partition segment into sub-segments
      vec<u64> boundaries;
      if (options.n_permutations == kCbsDefaultNPermutations &&
          math::IsEqualWithTolerance(options.max_p, kCbsDefaultMaxP)) {
        boundaries = kBoundaryP10000A01Eta05;
      }
      Logging::Debug("Segmenting on {} obs", sub_obvs.size());
      vec<Segment> sub_segments = GetSegments(sub_obvs,
                                              options.max_p,
                                              options.n_permutations,
                                              boundaries,
                                              options.min_obs_per_segment,
                                              rng,
                                              options.cbs_method,
                                              options.min_t_for_automatic_segmentation);
      Logging::Debug("Segmented on {} obs", sub_obvs.size());
      // merge back single-obs sub-segments
      if (options.no_single_obvs_subsegments) {
        sub_segments = MergeSingleObvsEndSegment(sub_segments, options.min_obs_per_segment);
      }
      if (sub_segments.size() <= 1) {
        // if this segment cannot be broken down, then return it as one segment
        Logging::Debug("adding un-breakable segment with {} obvs", (segment.end - segment.start));
        ret.push_back(segment);
      } else if (!options.disable_merging) {
        Logging::Debug("considering whether to merge back segments");
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
          if (!BinaryPermuteNullAndTestSignificance(sub_obvs,
                                                    options.n_permutations,
                                                    options.max_p,
                                                    seg_to_merge1,
                                                    seg_to_merge2,
                                                    options.min_t_for_automatic_segmentation,
                                                    rng)) {
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
            if (new_segment.end > n) {
              throw error::Error("cbs: sub-segment end exceeds observation count after coordinate adjustment");
            }
            Logging::Debug("pushing segment with {} obvs", (new_segment.end - new_segment.start));
            seg_stack.push(new_segment);
            if (seg_stack.size() > kCbsMaxStackDepth) {
              throw std::runtime_error("cbs: max depth exceeded!");
            }
          }
        }
      } else {
        for (auto new_segment : sub_segments) {
          new_segment.start += segment.start;
          new_segment.end += segment.start;
          if (new_segment.end > n) {
            throw error::Error("cbs: sub-segment end exceeds observation count after coordinate adjustment");
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
  std::ranges::sort(ret, [](const Segment& l, const Segment& r) { return l.start < r.start; });
  const auto time2 = high_resolution_clock::now();
  auto total_time = duration_cast<milliseconds>(time2 - time1).count();
  Logging::Debug("found {} segments over {} elements in {}ms", ret.size(), n, total_time);
  return ret;
}

/**
 * @brief Given a set of 2-3 subsegments, if the first or the last sub-segment has fewer than min_obs_per_segment
 * observations, then merge it with its adjacent sub-segment
 * @param sub_segments a std::vector of 2-3 sub-segments
 * @param min_obs_per_segment  minimum number of observations allowed per segment
 * @return
 */
vec<Segment> MergeSingleObvsEndSegment(const vec<Segment>& sub_segments, size_t min_obs_per_segment) {
  if (sub_segments.size() <= 1) {
    return sub_segments;
  }
  vec<Segment> res;
  if (sub_segments.size() < 4) {
    auto first_seg = sub_segments.begin();
    auto second_seg = sub_segments.begin() + 1;
    if (sub_segments.begin()->end - sub_segments.begin()->start < min_obs_per_segment) {
      // first segment is small
      first_seg = sub_segments.begin();
      second_seg = sub_segments.begin() + 1;
    } else if ((sub_segments.end() - 1)->end - (sub_segments.end() - 1)->start < min_obs_per_segment) {
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
    if (sub_segments.begin()->end - sub_segments.begin()->start < min_obs_per_segment) {
      // first segment is small
      res[0] = merged_segment;
      res[1] = sub_segments[2];
    } else if ((sub_segments.end() - 1)->end - (sub_segments.end() - 1)->start < min_obs_per_segment) {
      // last segment is small
      res[0] = sub_segments[0];
      res[1] = merged_segment;
    }
    return res;
  }
  Logging::Error("Will not prune segment with more than three sub-segments");
  throw std::runtime_error("Will not prune segment with more than three sub-segments");
}
}  // namespace xoos::cnc
