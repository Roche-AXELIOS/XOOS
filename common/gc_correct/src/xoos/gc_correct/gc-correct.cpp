#include "xoos/gc_correct/gc-correct.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <xoos/log/logging.h>
#include <xoos/types/int.h>

#include "xoos/gc_correct/loess.h"

namespace xoos::gc_correct {

// NOTE: This impl file uses Eigen internally. The header file avoids including Eigen in order to mitigate the cognitive
// complexity of using Eigen in user code. Eigen allows vec to be wrapped and have operations performed on the
// memory stored by vec directly.

/**
 * @brief determine which indices are able to be corrected for gc-correction
 * for GC correction
 * @param average_coverage    - average coverage vector (each index corresponds to a region)
 * @param gc_bias            - gc bias vector (each index corresponds to a region)
 * @return
 */
namespace {
vec<std::ptrdiff_t> GetValidRegions(const vec<f64>& average_coverage, const vec<f64>& gc_bias) {
  vec<std::ptrdiff_t> is_valid_indices;
  for (std::ptrdiff_t i = 0; i < ToSigned(average_coverage.size()); ++i) {
    if (average_coverage[ToUnsigned(i)] > 0.0 && gc_bias[ToUnsigned(i)] >= 0.0) {
      is_valid_indices.push_back(i);
    }
  }
  return is_valid_indices;
}

/**
 * @brief return the indices from which to train the loess model for gc correction. These are values that are in the
 * 0-99th quantile of the average coverage, and in the 0.001-0.999 quantile of the gc bias
 * @param average_coverage    - average coverage vector (each index corresponds to a region)
 * @param gc_bias            - gc bias vector (each index corresponds to a region)
 * @param is_valid_indices - list of indexes of "valid" regions. Use GetValidRegions to obtain this list
 * @return list of indexes of "ideal" regions
 */
vec<std::ptrdiff_t> GetIdealRegions(const vec<f64>& average_coverage,
                                    const vec<f64>& gc_bias,
                                    const vec<std::ptrdiff_t>& is_valid_indices) {
  vec<std::ptrdiff_t> is_ideal_indices;
  // calculate quantiles for avg coverage and gc bias
  vec<f64> valid_average_coverage(is_valid_indices.size());
  vec<f64> valid_gc_bias(is_valid_indices.size());
  for (std::ptrdiff_t i = 0; i < ToSigned(is_valid_indices.size()); ++i) {
    valid_average_coverage[i] = average_coverage[ToUnsigned(is_valid_indices.at(ToUnsigned(i)))];
    valid_gc_bias[i] = gc_bias[ToUnsigned(is_valid_indices.at(ToUnsigned(i)))];
  }
  // ranges for filtering out coverage and gc content outliers
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic) - contiguous buffer end iterator
  vec<f64> cov_vec(valid_average_coverage.data(), valid_average_coverage.data() + valid_average_coverage.size());
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic) - contiguous buffer end iterator
  vec<f64> gc_vec(valid_gc_bias.data(), valid_gc_bias.data() + valid_gc_bias.size());
  std::ranges::sort(cov_vec);
  std::ranges::sort(gc_vec);
  // Out-of-bounds access is possible here for edge cases with very few valid regions; not currently guarded against.
  // std::nth_element or another partial-sort method could improve performance over a full sort here.
  const f64 cov_min = cov_vec[0];
  const f64 cov_max = cov_vec[static_cast<size_t>(0.99 * static_cast<f64>(cov_vec.size()))];
  const f64 gc_min = gc_vec[static_cast<size_t>(0.001 * static_cast<f64>(gc_vec.size()))];
  const f64 gc_max = gc_vec[static_cast<size_t>((1 - 0.001) * static_cast<f64>(gc_vec.size()))];
  // gather "ideal" regions - i.e. valid regions w/o outliers
  for (const std::ptrdiff_t i : is_valid_indices) {
    if (average_coverage[ToUnsigned(i)] >= cov_min && average_coverage[ToUnsigned(i)] <= cov_max &&
        gc_bias[ToUnsigned(i)] >= gc_min && gc_bias[ToUnsigned(i)] <= gc_max) {
      is_ideal_indices.push_back(i);
    }
  }
  return is_ideal_indices;
}

GCCorrectLoessInput PrepareDataForLoess(const vec<f64>& average_coverage, const vec<f64>& gc_bias) {
  GCCorrectLoessInput gc_correct_loess_input;
  gc_correct_loess_input.is_valid_idxs = GetValidRegions(average_coverage, gc_bias);
  if (gc_correct_loess_input.is_valid_idxs.empty()) {
    Logging::Info("no valid regions found for GC correction");
    return {};
  }
  gc_correct_loess_input.is_ideal_idxs =
      GetIdealRegions(average_coverage, gc_bias, gc_correct_loess_input.is_valid_idxs);
  gc_correct_loess_input.ideal_average_coverage.resize(gc_correct_loess_input.is_ideal_idxs.size());
  gc_correct_loess_input.ideal_gc_bias.resize(gc_correct_loess_input.is_ideal_idxs.size());
  for (std::ptrdiff_t i = 0; i < ToSigned(gc_correct_loess_input.is_ideal_idxs.size()); ++i) {
    gc_correct_loess_input.ideal_average_coverage[i] =
        average_coverage[ToUnsigned(gc_correct_loess_input.is_ideal_idxs.at(ToUnsigned(i)))];
    gc_correct_loess_input.ideal_gc_bias[i] =
        gc_bias[ToUnsigned(gc_correct_loess_input.is_ideal_idxs.at(ToUnsigned(i)))];
  }
  gc_correct_loess_input.valid_average_coverage.resize(gc_correct_loess_input.is_valid_idxs.size());
  gc_correct_loess_input.valid_gc_bias.resize(gc_correct_loess_input.is_valid_idxs.size());
  for (std::ptrdiff_t i = 0; i < ToSigned(gc_correct_loess_input.is_valid_idxs.size()); ++i) {
    gc_correct_loess_input.valid_average_coverage[i] =
        average_coverage[ToUnsigned(gc_correct_loess_input.is_valid_idxs.at(ToUnsigned(i)))];
    gc_correct_loess_input.valid_gc_bias[i] =
        gc_bias[ToUnsigned(gc_correct_loess_input.is_valid_idxs.at(ToUnsigned(i)))];
  }
  return gc_correct_loess_input;
}

/**
 * @brief using total coverages and region sizes, calculate the mean per-base coverage for each region
 * @return Eigen::VectorXd of average coverages
 */
vec<f64> CalculateAverageCoverage(const vec<f64>& coverage, const vec<f64>& region_sizes) {
  vec<f64> avg_cov = coverage;
  // first, map all the vecs to Eigen arrays in order to take advantage of vectorized operations
  Eigen::Map<const Eigen::ArrayXd> region_sizes_v(region_sizes.data(), ToSigned(region_sizes.size()));
  Eigen::Map<const Eigen::ArrayXd> coverage_v(coverage.data(), ToSigned(coverage.size()));
  Eigen::Map<Eigen::ArrayXd> avg_cov_v(avg_cov.data(), ToSigned(avg_cov.size()));
  // then it's a simple element-wise division
  avg_cov_v = coverage_v / region_sizes_v;
  return avg_cov;
}

GCCorrectInput PrepareCoverageRecordsForGCCorrect(const vec<CoverageRecord>& coverage) {
  GCCorrectInput ret;
  ret.counts.resize(coverage.size());
  for (std::ptrdiff_t i = 0; i < ToSigned(coverage.size()); ++i) {
    ret.counts[i] = std::round(coverage.at(ToUnsigned(i)).read_count);
  }
  ret.total_coverage.resize(coverage.size());
  for (std::ptrdiff_t i = 0; i < ToSigned(coverage.size()); ++i) {
    ret.total_coverage[i] = std::round(coverage.at(ToUnsigned(i)).total_coverage);
  }
  ret.region_sizes.resize(coverage.size());
  for (std::ptrdiff_t i = 0; i < ToSigned(coverage.size()); ++i) {
    const auto& region = coverage.at(ToUnsigned(i)).region;
    ret.region_sizes[i] = static_cast<f64>(region.end - region.start);
  }
  ret.average_coverage = CalculateAverageCoverage(ret.total_coverage, ret.region_sizes);
  return ret;
}

}  // namespace

GCCorrectModel TrainGCCorrectModel(const vec<CoverageRecord>& coverage, const vec<f64>& gc_bias, const f64 first_span) {
  GCCorrectInput input = PrepareCoverageRecordsForGCCorrect(coverage);
  // Passing input.counts directly instead of calculating average_coverage might work equally well, but this is
  // unverified.
  const GCCorrectLoessInput loess_input = PrepareDataForLoess(input.average_coverage, gc_bias);
  // Guard against empty LOESS input: without any ideal or valid regions the downstream Loess() call would fail with a
  // misleading "too few data points per window" error that points at first_span rather than the real cause.
  if (loess_input.ideal_average_coverage.empty() || loess_input.is_valid_idxs.empty()) {
    Logging::Error("GC correction has no valid or ideal regions to train on (ideal regions: {}, valid regions: {})",
                   loess_input.ideal_average_coverage.size(),
                   loess_input.is_valid_idxs.size());
    throw std::runtime_error("no valid or ideal regions available for GC correction");
  }
  Eigen::VectorXd gc_to_est = Eigen::VectorXd::LinSpaced(
      static_cast<int>((kGCCorrectLoessGCMax - kGCCorrectLoessGCMin) / kGCCorrectLoessGCStep) + 1,
      kGCCorrectLoessGCMin,
      kGCCorrectLoessGCMax);
  // round 1 : find a representative coverage value for each GC bin (gc_to_est)
  auto ideal_gc_bias_v =
      Eigen::Map<const Eigen::VectorXd>(loess_input.ideal_gc_bias.data(), ToSigned(loess_input.ideal_gc_bias.size()));
  auto ideal_average_coverage_v = Eigen::Map<const Eigen::VectorXd>(
      loess_input.ideal_average_coverage.data(), ToSigned(loess_input.ideal_average_coverage.size()));
  Eigen::VectorXd smoothed_coverage_v =
      Loess(gc_to_est, ideal_gc_bias_v, ideal_average_coverage_v, first_span, kGCCorrectLoessDegree);
  // round 2: use the bins from round 1 to smooth all the input coverages
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic) - contiguous buffer end iterator
  return {.x = vec<f64>(gc_to_est.data(), gc_to_est.data() + gc_to_est.size()),
          // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic) - contiguous buffer end iterator
          .y = vec<f64>(smoothed_coverage_v.data(), smoothed_coverage_v.data() + smoothed_coverage_v.size()),
          .mean = ideal_average_coverage_v.mean(),
          .input = input,
          .loess_input = loess_input};
}

namespace {
GCCorrectCoverage ApplyGCCorrectionModel(const GCCorrectModel& model,
                                         const vec<std::ptrdiff_t>& valid_idxs,
                                         const vec<f64>& gc_bias,
                                         const GCCorrectInput& input) {
  Eigen::VectorXd x_to_predict(valid_idxs.size());
  for (std::ptrdiff_t i = 0; i < ToSigned(valid_idxs.size()); ++i) {
    x_to_predict[i] = gc_bias[ToUnsigned(valid_idxs.at(ToUnsigned(i)))];
  }
  auto x_v = Eigen::Map<const Eigen::VectorXd>(model.x.data(), ToSigned(model.x.size()));
  auto y_v = Eigen::Map<const Eigen::VectorXd>(model.y.data(), ToSigned(model.y.size()));
  auto round2_smoothed_cov = Loess(x_to_predict, x_v, y_v, kGCCorrectSecondLoessSpan, kGCCorrectLoessDegree);
  Eigen::VectorXd scale_factor = round2_smoothed_cov / model.mean;
  s64 j = 0;
  vec<CoverageRecord> coverage_records(input.average_coverage.size());
  for (const s64 i : valid_idxs) {
    const auto idx = ToUnsigned(i);
    if (scale_factor[j] > 0.0) {
      coverage_records[idx].read_count = input.counts[idx] / scale_factor[j];
      coverage_records[idx].total_coverage = input.total_coverage[idx] / scale_factor[j];
      coverage_records[idx].average_coverage = coverage_records[idx].total_coverage / input.region_sizes[idx];
    } else {
      coverage_records[idx].read_count = 0;
      coverage_records[idx].average_coverage = 0;
      coverage_records[idx].total_coverage = 0;
    }
    j += 1;
  }
  return {.corrected_coverage = coverage_records,
          .gc_bins = {},
          .smoothed_coverage_at_gc_bins = {},
          // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic) - contiguous buffer end iterator
          .scale_factor = vec<f64>(scale_factor.data(), scale_factor.data() + scale_factor.size()),
          .valid_indices = valid_idxs};
}

// Restore genomic coordinates onto GC-corrected records. ApplyGCCorrectionModel default-constructs the output records
// and populates only the numeric fields, so the region of each record must be copied back from the input coverage.
void BackfillRegions(vec<CoverageRecord>& corrected, const vec<CoverageRecord>& source) {
  if (corrected.size() != source.size()) {
    Logging::Error(
        "size of corrected coverage does not match size of input coverage: {} vs {}", corrected.size(), source.size());
    throw std::runtime_error("mismatch between size of corrected coverage and size of input coverage");
  }
  for (size_t i = 0; i < corrected.size(); ++i) {
    corrected[i].region = source[i].region;
  }
}
}  // namespace

/**
 * @brief correct counts and total coverage values for GC bias using two rounds of LOESS smoothing.
 * @param coverage      - vector of CoverageRecord structs
 * @param gc_bias      - vector of GC bias values (each index corresponds to a region). Must be same size as coverage
 * @param first_span   - span parameter for first round of LOESS smoothing. Cannot be 0. Represents the percent of data
 * needed used to smooth each point. For small datasets, use 1.0
 * @return a GCCorrectRet struct, which contains a vector of corrected CoverageRecord structs and the scale factors used
 */
GCCorrectCoverage GCCorrect(const vec<CoverageRecord>& coverage, const vec<f64>& gc_bias, const f64 first_span) {
  auto model = TrainGCCorrectModel(coverage, gc_bias, first_span);
  auto ret = ApplyGCCorrectionModel(model, model.loess_input.is_valid_idxs, gc_bias, model.input);
  BackfillRegions(ret.corrected_coverage, coverage);
  return ret;
}

GCCorrectCoverage GCCorrectFromPrecomputedModel(const GCCorrectModel& model,
                                                const vec<CoverageRecord>& coverage,
                                                const vec<f64>& gc_bias) {
  GCCorrectInput input = PrepareCoverageRecordsForGCCorrect(coverage);
  auto loess_input = PrepareDataForLoess(input.average_coverage, gc_bias);
  auto ret = ApplyGCCorrectionModel(model, loess_input.is_valid_idxs, gc_bias, input);
  BackfillRegions(ret.corrected_coverage, coverage);
  return ret;
}

}  // namespace xoos::gc_correct
