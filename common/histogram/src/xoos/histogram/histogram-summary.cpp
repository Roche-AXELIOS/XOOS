#include "histogram-summary.h"

#include <optional>
#include <string>

#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::histogram {

const std::string kNA = "NA";
const vec<u64> kDefaultPercentiles{10, 25, 50, 75, 90};

const HistogramSummaryConfig kDefaultCoverageHistogramSummaryConfig = {.mean = "mean_coverage",
                                                                       .min = "min_coverage",
                                                                       .max = "max_coverage",
                                                                       .median = "median_coverage",
                                                                       .zero = "zero_bases",
                                                                       .non_zero = "non_zero_bases"};

const HistogramSummaryConfig kDefaultReadLengthHistogramSummaryConfig = {.mean = "mean_read_length",
                                                                         .min = "min_read_length",
                                                                         .max = "max_read_length",
                                                                         .median = "median_read_length",
                                                                         .output_zero = false,
                                                                         .output_non_zero = false,
                                                                         .output_percentile_ratios = false};

std::optional<f64> GetRatioPercentile(const PercentileMap& percentiles, const u64 numerator, const u64 denominator) {
  // check if the numerator and denominator are valid percentiles
  if (!percentiles.contains(numerator) || !percentiles.contains(denominator)) {
    return std::nullopt;
  }
  const std::optional<u64>& numerator_opt = percentiles.at(numerator);
  const std::optional<u64>& denominator_opt = percentiles.at(denominator);
  if (!numerator_opt.has_value() || !denominator_opt.has_value()) {
    return std::nullopt;
  }
  // if the denominator is zero, return std::nullopt to avoid division by zero
  if (denominator_opt.value() == 0) {
    return std::nullopt;
  }
  return static_cast<f64>(numerator_opt.value()) / static_cast<f64>(denominator_opt.value());
}

bool OmitsFirstBin(const HistogramBinOutput mode) {
  return mode == HistogramBinOutput::kOmitFirstBin ||
         mode == HistogramBinOutput::kOmitFirstBinAndMaxLastBinWithOutlier ||
         mode == HistogramBinOutput::kOmitFirstBinAndMaxLastBin;
}

bool FoldsOutliersIntoLastBin(const HistogramBinOutput mode) {
  return mode == HistogramBinOutput::kMaxLastBin || mode == HistogramBinOutput::kOmitFirstBinAndMaxLastBin;
}

}  // namespace xoos::histogram
