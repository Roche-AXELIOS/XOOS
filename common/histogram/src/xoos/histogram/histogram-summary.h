#pragma once

#include <algorithm>
#include <limits>
#include <map>
#include <optional>
#include <ranges>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>

#include <csv.hpp>

#include <xoos/histogram/histogram.h>
#include <xoos/io/metadata-util.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::histogram {

extern const std::string kNA;
extern const vec<u64> kDefaultPercentiles;
constexpr u8 kMaxPercentile = 100;
constexpr u8 kMedianPercentile = 50;
// lower and upper percentiles used to compute percentile ratios (e.g. ratio_90_to_10percentile)
constexpr u8 kLowerRatioPercentile = 10;
constexpr u8 kUpperRatioPercentile = 90;

// HistogramBinOutput defines how histogram bins should be outputted when writing to a file
enum class HistogramBinOutput : u8 {
  // Default behavior: output all bins individually
  kDefault,
  // Omit the first bin
  kOmitFirstBin,
  // Add "+" to the last bin (when the last bin contains all values greater than or equal to the last bin)
  kMaxLastBin,
  // Omit first bin AND add "+" to the last bin
  kOmitFirstBinAndMaxLastBin,
  // Omit first bin AND the last bin is the outlier bin (bin size with "+" suffix)
  kOmitFirstBinAndMaxLastBinWithOutlier
};

// key is percentile, value is the count
using PercentileMap = std::map<u64, std::optional<u64>>;

// key is cutoff, value is the percentage of counts for bins that are greater than or equal to the given cutoff
using PercentageCutoffMap = std::map<u64, std::optional<f64>>;

// Note: remove this functionality when done with refactoring alignment metrics
struct HistogramSummaryConfig {
  std::string mean = "mean";
  std::string min = "min";
  std::string max = "max";
  std::string median = "median";
  std::string zero = "zero";
  std::string non_zero = "non_zero";
  bool output_mean{true};
  bool output_min{true};
  bool output_max{true};
  bool output_median{true};
  bool output_zero{true};
  bool output_non_zero{true};
  bool output_percentiles{true};
  bool output_percentile_ratios{true};
  bool output_cutoff_percentages{true};
};

// Provides a summary with a variable number of percentiles and cutoff percentages
// std::optional is used to represent invalid values as std::nullopt (e.g. empty histogram stats or divide by zero
// result)
template <typename T>
  requires Numeric<T>
struct HistogramSummary {
  // count is the sum of all bin values (of type T) in the histogram
  T count;
  std::optional<u64> median;
  std::optional<f64> mean;
  std::optional<f64> stddev;
  std::optional<u64> min;
  std::optional<u64> max;

  // Note: remove this functionality when done with refactoring alignment metrics
  u64 count_zero = 0;
  u64 count_non_zero = 0;

  std::map<u64, std::optional<u64>> percentiles;
  std::map<u64, std::optional<f64>> cutoff_percentage;

  // Note: remove this functionality when done with refactoring alignment metrics
  HistogramSummaryConfig config{};
};

// A vector representing histogram summaries for distributions that share the same binning
// Each tuple contains the name of the histogram and the summary
template <typename T>
  requires Numeric<T>
using HistogramSummaries = vec<std::tuple<std::string, HistogramSummary<T>>>;

/**
 * Compute the value of all data points (i.e. the sum of all bin values multiplied by their respective bin index) in the
 * histogram.
 */
template <typename T>
  requires Numeric<T>
T ComputeValueSum(const Histogram<T>& histogram, const bool exclude_zero_bins = false) {
  T total_value_sum = 0;
  for (u64 bin_index = 0; std::cmp_less(bin_index, histogram.counts.size()); ++bin_index) {
    const T bin_count = histogram.At(bin_index);
    if (exclude_zero_bins && bin_index == 0) {
      continue;
    }
    total_value_sum += bin_count * bin_index;
  }
  // count outliers
  for (const auto& [key, value] : histogram.outliers) {
    total_value_sum += value * key;
  }
  return total_value_sum;
}

/**
 * Compute the total number of data points (i.e. the sum of all bin values) in the histogram.
 */
template <typename T>
  requires Numeric<T>
T ComputeCount(const Histogram<T>& histogram, const bool exclude_zero_bins = false) {
  T total_count = 0;
  for (u64 bin_index = 0; std::cmp_less(bin_index, histogram.counts.size()); ++bin_index) {
    const T bin_count = histogram.At(bin_index);
    if (exclude_zero_bins && bin_index == 0) {
      continue;
    }
    total_count += bin_count;
  }
  // count outliers
  for (const auto& [key, value] : histogram.outliers) {
    total_count += value;
  }
  return total_count;
}

/**
 *  Compute the stddev bin of the histogram.
 *  We calculate it in the same loop as the (unrounded) mean to avoid multiple passes over the histogram.
 */
template <typename T>
  requires Numeric<T>
f64 ComputeStddev(const Histogram<T>& histogram, const u64 start_index = 0) {
  if (histogram.IsEmpty(start_index)) {
    throw std::invalid_argument("Cannot compute stddev of an empty histogram");
  }
  // we have to recalculate the mean since we need it for stddev calculation and we don't want it rounded
  T weighted_sum = 0;
  T total_count = 0;
  for (u64 bin_index = start_index; std::cmp_less(bin_index, histogram.counts.size()); ++bin_index) {
    const T bin_count = histogram.At(bin_index);
    total_count += bin_count;
    weighted_sum += 1LL * bin_index * bin_count;
  }
  // count outliers
  for (const auto& [bin_index, bin_count] : histogram.outliers) {
    total_count += bin_count;
    weighted_sum += 1LL * bin_index * bin_count;
  }

  // calculate variance
  const f64 mean = static_cast<f64>(weighted_sum) / static_cast<f64>(total_count);
  f64 variance_sum = 0.0;
  for (u64 bin_index = start_index; std::cmp_less(bin_index, histogram.counts.size()); ++bin_index) {
    const T bin_count = histogram.At(bin_index);
    const f64 deviation = static_cast<f64>(bin_index) - mean;
    variance_sum += bin_count * deviation * deviation;
  }
  // count outliers
  for (const auto& [bin_index, bin_count] : histogram.outliers) {
    const f64 deviation = static_cast<f64>(bin_index) - mean;
    variance_sum += bin_count * deviation * deviation;
  }
  const f64 variance = variance_sum / static_cast<f64>(total_count);
  return std::sqrt(variance);
}

/**
 * Compute the mean bin of the histogram.
 *
 * NOTE: The histogram is assumed to be non-empty (at least one bin has a non-zero count)
 */
template <typename T>
  requires Numeric<T>
f64 ComputeMean(const Histogram<T>& histogram, u64 start_index) {
  if (histogram.IsEmpty(start_index)) {
    throw std::invalid_argument("Cannot compute mean of an empty histogram");
  }
  T weighted_sum = 0;
  T total_count = 0;
  for (u64 bin_index = start_index; std::cmp_less(bin_index, histogram.counts.size()); ++bin_index) {
    const T bin_count = histogram.At(bin_index);
    total_count += bin_count;
    weighted_sum += 1LL * bin_index * bin_count;
  }
  // count outliers
  for (const auto& [bin_index, bin_count] : histogram.outliers) {
    total_count += bin_count;
    weighted_sum += 1LL * bin_index * bin_count;
  }
  return static_cast<f64>(weighted_sum) / static_cast<f64>(total_count);
}

template <typename T>
  requires Numeric<T>
f64 ComputeMean(const Histogram<T>& histogram) {
  return ComputeMean(histogram, 0);
}

/**
 * Compute the min bin of the histogram. The min bin is the first bin with non-zero count.
 * If there is no bin with non-zero count, it returns the first outlier bin.
 * If there is no outlier bin, it returns the last bin index (which is the size of the histogram counts).
 * NOTE: The histogram is assumed to be non-empty (at least one bin has a non-zero count)
 */
template <typename T>
  requires Numeric<T>
u64 ComputeMin(const Histogram<T>& histogram, u64 start_index) {  // NOSONAR(S5536) template instantiated in other TUs
  if (histogram.IsEmpty(start_index)) {
    throw std::invalid_argument("Cannot compute min of an empty histogram");
  }

  // search counts only when start_index is within bounds
  if (start_index < histogram.counts.size()) {
    auto it = std::ranges::find_if(histogram.counts.begin() + ToSigned(start_index),
                                   histogram.counts.end(),
                                   [](const T bin_count) { return bin_count != 0; });
    auto min = std::ranges::distance(histogram.counts.begin(), it);
    if (std::cmp_less(min, histogram.counts.size())) {
      return static_cast<u64>(min);
    }
  }

  // check outliers at or after start_index
  auto outlier_it = histogram.outliers.lower_bound(start_index);
  if (outlier_it != histogram.outliers.end()) {
    return outlier_it->first;
  }
  return static_cast<u64>(histogram.counts.size());
}

/// @overload ComputeMin with start_index defaulting to 0.
template <typename T>
  requires Numeric<T>
u64 ComputeMin(const Histogram<T>& histogram) {  // NOSONAR(S5536) template instantiated in other TUs
  return ComputeMin(histogram, 0);
}

/**
 * Compute the max bin of the histogram. The max bin is the last bin with non-zero count.
 * NOTE: The histogram is assumed to be non-empty (at least one bin has a non-zero count) and max_counts is the last
 * histogram.
 */
template <typename T>
  requires Numeric<T>
u64 ComputeMax(const Histogram<T>& histogram) {
  if (histogram.IsEmpty()) {
    throw std::invalid_argument("Cannot compute max of an empty histogram");
  }
  // find the last bin with non-zero count by reversing the histogram
  auto it_reversed = std::ranges::find_if(histogram.counts | std::ranges::views::reverse,
                                          [](const T bin_count) { return bin_count != 0; });
  auto max = std::ranges::distance(it_reversed, histogram.counts.rend() - 1);
  // if the histogram has outliers, we need to check if the last outlier exceeds the counts range
  if (!histogram.outliers.empty() && histogram.outliers.rbegin()->first >= histogram.counts.size()) {
    return histogram.outliers.rbegin()->first;
  }
  return static_cast<u64>(max);
}

/**
 * Compute the percentile of the histogram.
 *
 * NOTE: The histogram is assumed to be non-empty (at least one bin has a non-zero count). The percentile must be
 * between 0 and 100.
 *
 * @param histogram a vector representing a histogram where the index is the bin and the value at a given index is the
 * count
 * @param percentile the percentile to calculate
 */
template <typename T>
  requires Numeric<T>
u64 ComputePercentile(const Histogram<T>& histogram, const u8 percentile, u64 start_index = 0) {
  if (histogram.IsEmpty(start_index)) {
    throw std::invalid_argument("Cannot compute percentile of an empty histogram");
  }
  if (percentile > kMaxPercentile) {
    throw std::out_of_range(fmt::format("Percentile must be between 0 and {}", kMaxPercentile));
  }
  // create new histogram with the starting index
  // this is useful when we want to compute percentiles for histograms that start at a specific index
  // e.g. when we want to compute percentiles for histograms that start at 2, we need to create a new histogram
  // with the starting index of 2
  histogram::Histogram<T> updated_histogram_with_new_start_index;
  updated_histogram_with_new_start_index = histogram::Histogram<T>(histogram.counts.size() - start_index, 0);
  for (u64 bin_index = start_index; std::cmp_less(bin_index, histogram.counts.size()); ++bin_index) {
    updated_histogram_with_new_start_index.At(bin_index - start_index) = histogram.At(bin_index);
  }
  // add outliers to the new histogram
  for (const auto& [bin_index, bin_count] :
       histogram.outliers) {  // NOSONAR(S5311) structured bindings are idiomatic C++17
    updated_histogram_with_new_start_index.outliers[bin_index - start_index] = bin_count;
  }

  u64 percentile_bin = 0;
  if (percentile == 0) {
    percentile_bin = ComputeMin(updated_histogram_with_new_start_index);
  } else if (percentile == kMaxPercentile) {
    percentile_bin = ComputeMax(updated_histogram_with_new_start_index);
  } else {
    T total_count = ComputeCount(updated_histogram_with_new_start_index);
    T cumulative_count = 0;
    for (u64 bin_index = 0; std::cmp_less(bin_index, updated_histogram_with_new_start_index.counts.size());
         ++bin_index) {
      cumulative_count += updated_histogram_with_new_start_index.At(bin_index);
      if (cumulative_count >= static_cast<f64>(total_count) * percentile / static_cast<f64>(kMaxPercentile)) {
        // add the start index to the percentile bin to get the correct bin index
        return bin_index + start_index;
      }
    }
    // handle outliers
    for (const auto& [bin_index, bin_count] : updated_histogram_with_new_start_index.outliers) {
      cumulative_count += bin_count;
      if (cumulative_count >= static_cast<f64>(total_count) * percentile / static_cast<f64>(kMaxPercentile)) {
        // add the start index to the percentile bin to get the correct bin index
        return bin_index + start_index;
      }
    }
  }
  // add the start index to the percentile bin to get the correct bin index
  return percentile_bin + start_index;
}

/**
 * Compute percentage of counts for bins that are greater than or equal to the given cutoff bin.
 *
 * NOTE: The histogram is assumed to be non-empty (at least one bin has a non-zero count)
 *
 * @param histogram a vector representing a histogram where the index is the bin and the value at a given index is the
 * count
 * @param cutoff cutoff value
 */
template <typename T>
  requires Numeric<T>
f64 ComputeCutoffPercentage(const Histogram<T>& histogram, const u64 cutoff) {
  if (histogram.IsEmpty()) {
    throw std::invalid_argument("Cannot compute cutoff percentage of an empty histogram");
  }
  const u64 total_count = ComputeCount(histogram);
  u64 total_count_above_cutoff = 0;
  // if bin_index is greater than or equal to cutoff, we add the count to the total count above cutoff
  for (u64 bin_index = 0; std::cmp_less(bin_index, histogram.counts.size()); ++bin_index) {
    const u64 bin_count = histogram.counts.at(bin_index);
    if (bin_index >= cutoff) {
      total_count_above_cutoff += bin_count;
    }
  }
  for (const auto& [bin_index, bin_count] : histogram.outliers) {
    if (bin_index >= cutoff) {
      total_count_above_cutoff += bin_count;
    }
  }
  // calculate percentage — cast before multiplying to avoid u64 overflow
  const auto percentage_above_cutoff =
      static_cast<f64>(total_count_above_cutoff) * static_cast<f64>(kMaxPercentile) / static_cast<f64>(total_count);
  return percentage_above_cutoff;
}

/**
 * Calculate summary statistics of the histogram.
 *
 * It calculates count, median, mean, min, max, percentiles and cutoff percentage. If the histogram is empty (all bins
 * are zero), then the respective metrics are represented by std::nullopt.
 *
 * @param histogram a vector representing a histogram where the index is the bin and the value at a given index is the
 * count.
 * @param cutoffs cutoff values to calculate the percentage
 * @param percentiles percentiles to calculate
 * @param exclude_zero_bins if true, exclude zero bins when calculating percentiles (i.e. start from bin index 1)
 */
template <typename T>
  requires Numeric<T>
HistogramSummary<T> CalculateHistogramSummary(const Histogram<T>& histogram,
                                              const vec<u64>& cutoffs = {},
                                              const vec<u64>& percentiles = kDefaultPercentiles,
                                              bool exclude_zero_bins = false) {
  T count = ComputeCount(histogram, exclude_zero_bins);
  // If the histogram is empty (either because original histogram is empty or because we are excluding zero bins and the
  // histogram only has zero bins), return empty summary with std::nullopt values
  if (count == 0) {
    PercentileMap empty_percentiles;
    for (const auto& percentile : percentiles) {
      empty_percentiles[percentile] = std::nullopt;
    }
    PercentageCutoffMap empty_cutoffs;
    for (const auto& cutoff : cutoffs) {
      empty_cutoffs[cutoff] = std::nullopt;
    }
    // count is still returned as the total count (including zero bins)
    return HistogramSummary<T>{.count = count,
                               .median = std::nullopt,
                               .mean = std::nullopt,
                               .stddev = std::nullopt,
                               .min = std::nullopt,
                               .max = std::nullopt,
                               .percentiles = std::move(empty_percentiles),
                               .cutoff_percentage = std::move(empty_cutoffs)};
  }

  const u64 median = exclude_zero_bins ? ComputePercentile(histogram, kMedianPercentile, 1)
                                       : ComputePercentile(histogram, kMedianPercentile);

  const f64 mean = ComputeMean(histogram, exclude_zero_bins ? 1 : 0);
  const f64 stddev = ComputeStddev(histogram, exclude_zero_bins ? 1 : 0);
  const u64 min = ComputeMin(histogram);
  const u64 max = ComputeMax(histogram);

  PercentileMap percentiles_map;
  for (const auto& percentile : percentiles) {
    percentiles_map[percentile] =
        exclude_zero_bins ? ComputePercentile(histogram, percentile, 1) : ComputePercentile(histogram, percentile);
  }

  PercentageCutoffMap cutoff_percentage_map;
  for (const auto& cutoff : cutoffs) {
    cutoff_percentage_map[cutoff] = ComputeCutoffPercentage(histogram, cutoff);
  }

  return HistogramSummary<T>{.count = count,
                             .median = median,
                             .mean = mean,
                             .stddev = stddev,
                             .min = min,
                             .max = max,
                             .percentiles = std::move(percentiles_map),
                             .cutoff_percentage = std::move(cutoff_percentage_map)};
}

/**
 * Calculates the ratio of two specified percentiles from a PercentileMap
 *
 * NOTE: Handles invalid operations (i.e. division by zero, or bad access) by returning std::nullopt.
 * */
std::optional<f64> GetRatioPercentile(const PercentileMap& percentiles, u64 numerator, u64 denominator);

// helper function to get string with 2 decimal places
template <typename T>
  requires Numeric<T>
std::string ValueToString(const T& val, const u8 precision = 0) {
  return fmt::format("{:.{}f}", static_cast<f64>(val), precision);
}

// helper function to convert optional values to appropriate string
// must use stringstream to get the correct precision
template <typename T>
  requires Numeric<T>
std::string ValueToStringOrNA(const std::optional<T>& val, const u8 precision = 0) {
  if (val.has_value()) {
    return ValueToString(val.value(), precision);
  }
  return kNA;
}

/**
 * @brief Append a formatted numeric value to a TSV row.
 *
 * Converts @p value to a string and appends it to @p row. If @p precisions is non-empty, the value is formatted with
 * the precision at @p precision_index (which is then incremented). If @p precisions is empty, the value is formatted
 * with default precision (0 decimal places).
 *
 * This is a helper used by WriteHistogramBinRows, WriteOutlierRows, and WriteFoldedOutlierRow to keep formatting
 * logic in one place.
 *
 * @param row the row to append the formatted value to
 * @param value the numeric value to format and append
 * @param precisions per-histogram decimal precision values; empty means use default precision
 * @param precision_index current index into @p precisions, incremented after use
 */
template <typename T>
  requires Numeric<T>
void AppendValue(vec<std::string>& row,  // NOSONAR(S5536)
                 const T value,
                 const vec<u8>& precisions,
                 size_t& precision_index) {
  if (precisions.empty()) {
    row.push_back(ValueToString<T>(value));
  } else {
    const auto precision = precisions.at(precision_index);
    ++precision_index;
    row.push_back(ValueToString<T>(value, precision));
  }
}

/**
 * @brief Check whether the output mode skips bin 0.
 *
 * Returns true for kOmitFirstBin, kOmitFirstBinAndMaxLastBinWithOutlier, and kOmitFirstBinAndMaxLastBin. These modes
 * are used when bin 0 represents an uninteresting or dominant category (e.g., zero-coverage regions) that would
 * obscure the rest of the distribution.
 *
 * @param mode the histogram bin output mode to check
 * @return true if the mode omits the first bin
 */
bool OmitsFirstBin(HistogramBinOutput mode);

/**
 * @brief Check whether the output mode folds outlier bins into the last regular bin.
 *
 * Returns true for kMaxLastBin and kOmitFirstBinAndMaxLastBin. In these modes, the last bin label is suffixed with "+"
 * and its count includes all outlier values at or beyond that bin index. This produces a compact output where the tail
 * of the distribution is collapsed into a single row.
 *
 * @param mode the histogram bin output mode to check
 * @return true if the mode folds outliers into the last bin
 */
bool FoldsOutliersIntoLastBin(HistogramBinOutput mode);

/**
 * @brief Write the main histogram bin rows to a TSV writer.
 *
 * Iterates over each bin index from 0 to @p num_bins and writes one TSV row per bin. The first column is the bin label
 * and subsequent columns are the histogram values (one per histogram in @p histograms).
 *
 * Bin output behavior depends on @p histogram_bin_output:
 * - If OmitsFirstBin() is true, bin 0 is skipped.
 * - If FoldsOutliersIntoLastBin() is true, the last bin label is suffixed with "+" and its value includes all outlier
 *   counts at or beyond that bin index.
 *
 * @param writer a csv::TSVWriter to write rows to
 * @param histograms the histograms to write (all must share the same bin count)
 * @param num_bins the number of bins (from the first histogram's counts vector)
 * @param histogram_bin_output controls bin skipping and outlier folding behavior
 * @param histogram_precisions per-histogram decimal precision values; empty means default precision
 */
template <typename T>
  requires Numeric<T>
void WriteHistogramBinRows(auto& writer,  // NOSONAR(S5536)
                           const Histograms<T>& histograms,
                           const u64 num_bins,
                           const HistogramBinOutput histogram_bin_output,
                           const vec<u8>& histogram_precisions) {
  for (u64 bin_index = 0; bin_index < num_bins; ++bin_index) {
    if (OmitsFirstBin(histogram_bin_output) && bin_index == 0) {
      continue;
    }
    auto bin = std::to_string(bin_index);
    if (bin_index == num_bins - 1 && FoldsOutliersIntoLastBin(histogram_bin_output)) {
      bin += "+";
    }
    vec row{bin};
    size_t precision_index = 0;
    for (const auto& [_, histogram] : histograms) {
      auto histogram_value = histogram.At(bin_index);
      if (bin_index == num_bins - 1 && FoldsOutliersIntoLastBin(histogram_bin_output)) {
        for (const auto& [outlier_bin_index, outlier_value] : histogram.outliers) {
          if (outlier_bin_index >= bin_index) {
            histogram_value += outlier_value;
          }
        }
      }
      AppendValue(row, histogram_value, histogram_precisions, precision_index);
    }
    writer << row;
  }
}

/**
 * @brief Write explicit outlier bin rows to a TSV writer.
 *
 * Used when the output mode is kDefault or kOmitFirstBin — modes where outliers are not folded into the last regular
 * bin. Writes one row per outlier bin index, from @p num_bins up to the maximum outlier bin index across all
 * histograms. Bins with no outlier entry for a given histogram are written as zero.
 *
 * @param writer a csv::TSVWriter to write rows to
 * @param histograms the histograms whose outlier maps are iterated
 * @param num_bins the number of regular bins (outlier rows start at this index)
 * @param histogram_precisions per-histogram decimal precision values; empty means default precision
 */
template <typename T>
  requires Numeric<T>
void WriteOutlierRows(auto& writer,  // NOSONAR(S5536)
                      const Histograms<T>& histograms,
                      const u64 num_bins,
                      const vec<u8>& histogram_precisions) {
  u64 max_outlier_bin_index = 0;
  for (const auto& [_, histogram] : histograms) {
    for (const auto& [outlier_bin_index, _] : histogram.outliers) {
      max_outlier_bin_index = std::max(max_outlier_bin_index, outlier_bin_index);
    }
  }
  for (u64 bin_index = num_bins; bin_index < max_outlier_bin_index + 1; ++bin_index) {
    vec<std::string> row{};
    row.push_back(ValueToString<u64>(bin_index));
    // Each histogram contributes one column; missing outlier entries are written as zero.
    size_t precision_index = 0;
    for (const auto& [_, histogram] : histograms) {
      T bin_value = T{0};
      if (histogram.outliers.contains(bin_index)) {
        bin_value = histogram.outliers.at(bin_index);
      }
      AppendValue(row, bin_value, histogram_precisions, precision_index);
    }
    writer << row;
  }
}

/**
 * @brief Write a single summary row that folds all outliers into a "N+" bin.
 *
 * Used when the output mode is kOmitFirstBinAndMaxLastBinWithOutlier. Writes one row labeled "{num_bins}+" whose value
 * is the sum of all outlier counts across all outlier bin indices for each histogram. This differs from
 * FoldsOutliersIntoLastBin modes where outliers are added to the last regular bin — here they get their own separate
 * row.
 *
 * @param writer a csv::TSVWriter to write the row to
 * @param histograms the histograms whose outlier maps are summed
 * @param num_bins the number of regular bins (used to label the row as "{num_bins}+")
 * @param histogram_precisions per-histogram decimal precision values; empty means default precision
 */
template <typename T>
  requires Numeric<T>
void WriteFoldedOutlierRow(auto& writer,  // NOSONAR(S5536)
                           const Histograms<T>& histograms,
                           const u64 num_bins,
                           const vec<u8>& histogram_precisions) {
  const auto bin = fmt::format("{}+", num_bins);
  vec row{bin};
  // Each histogram contributes one column with the sum of all its outlier values.
  size_t precision_index = 0;
  for (const auto& [_, histogram] : histograms) {
    auto histogram_value = T{0};
    for (const auto& [outlier_bin_index, outlier_value] : histogram.outliers) {
      histogram_value += outlier_value;
    }
    AppendValue(row, histogram_value, histogram_precisions, precision_index);
  }
  writer << row;
}

/**
 * @brief Write histogram data to a TSV file with structured TsvMetadata.
 *
 * Writes the histogram(s) to a file in tab-separated format. Each row is a bin and each column is a value from the
 * respective histogram (aside from the first column, which is for the bins themselves). This provides the user with a
 * single file for multiple distributions (if applicable). Precisions of each column can be specified with
 * histogram_precisions.
 *
 * NOTE: All histograms must share the same x-axis (bins). Additional formatting can be specified with
 * histogram_bin_output (i.e. omit the first bin, add a "+" to the last bin).
 *
 * Writes the metadata comment header (a `##RocheCommandLine=<...>` line and, when sample information is
 * present, a `##Sample=<...>` line) via io::WriteTsvMetadata.
 *
 * Example output:
 * @code
 * ##RocheCommandLine=<ID=tool_name,Version="1.0.0",CommandLine="tool_name --input foo.bam">
 * ##Sample=<Name=SAMPLE_A,LibraryId=SBX-D.SAMPLE_A,Platform=SBX,PlatformModel=AXELIOS_1>
 * cluster_size  total_clusters
 * 0  0
 * 1  5
 * 2  3
 * 3+  1
 * @endcode
 *
 * @param histograms a vector of (name, Histogram) tuples; all histograms must share the same bin count
 * @param output_histogram the path to the output TSV file
 * @param binning_name the label for the first column (e.g., "depth", "cluster_size")
 * @param histogram_bin_output controls bin skipping and outlier handling (see HistogramBinOutput)
 * @param histogram_precisions per-histogram decimal precision values; empty means default precision. Must be empty or
 *        have exactly one entry per histogram.
 * @param tsv_metadata command line info and optional sample information written as the comment header
 *
 * @throws std::invalid_argument if histogram_precisions is non-empty and its size does not match the number of
 *         histograms
 */
template <typename T>
  requires Numeric<T>
void WriteHistogramsToTsv(const Histograms<T>& histograms,  // NOSONAR(S5536)
                          const fs::path& output_histogram,
                          const std::string& binning_name,
                          const HistogramBinOutput histogram_bin_output,
                          const vec<u8>& histogram_precisions,
                          const io::TsvMetadata& tsv_metadata) {
  if (histograms.empty()) {
    return;
  }
  if (!histogram_precisions.empty() && histogram_precisions.size() != histograms.size()) {
    throw std::invalid_argument("Histogram precisions must be empty or match the number of histograms");
  }
  auto out = std::ofstream{output_histogram};
  auto writer = csv::make_tsv_writer(out).set_auto_flush(false);
  io::WriteTsvMetadata(out, tsv_metadata);
  const auto num_bins = std::get<1>(histograms.at(0)).counts.size();
  for (const auto& [name, histogram] : histograms) {
    if (histogram.counts.size() != num_bins) {
      throw std::invalid_argument(
          fmt::format("All histograms must have the same bin count. '{}' has {} bins, expected {}",
                      name,
                      histogram.counts.size(),
                      num_bins));
    }
  }

  vec header{binning_name};
  for (const auto& [histogram_name, histogram] : histograms) {
    header.push_back(histogram_name);
  }
  writer << header;

  WriteHistogramBinRows(writer, histograms, num_bins, histogram_bin_output, histogram_precisions);

  if (histogram_bin_output == HistogramBinOutput::kDefault ||
      histogram_bin_output == HistogramBinOutput::kOmitFirstBin) {
    WriteOutlierRows(writer, histograms, num_bins, histogram_precisions);
  }

  if (histogram_bin_output == HistogramBinOutput::kOmitFirstBinAndMaxLastBinWithOutlier) {
    WriteFoldedOutlierRow(writer, histograms, num_bins, histogram_precisions);
  }
}

constexpr u32 kMeanPrecision = 2;
constexpr u32 kStddevPrecision = 2;

/**
 * @brief Write summary histogram data to a TSV file with structured TsvMetadata.
 *
 * Writes the summary of the histogram(s) to a file in tab-separated format. Each row is a metric and each column is a
 * histogram (aside from the first column, which is for the metric names). This provides the user with a single, summary
 * report for multiple distributions (if applicable).
 *
 * NOTE: All histograms must share the same x-axis (bins). If the histogram is empty (all bins are zero), then some
 * statistics will be represented by the constant kNA (i.e. not applicable). Integer statistics (count, median,
 * min, max, percentiles) are written without decimal places. Mean and stddev are written with kMeanPrecision /
 * kStddevPrecision decimal places (currently 2). Cutoff percentages are written with 2 decimal places.
 *
 * Writes the metadata comment header (a `##RocheCommandLine=<...>` line and, when sample information is
 * present, a `##Sample=<...>` line) via io::WriteTsvMetadata.
 *
 * Example output:
 * @code
 * ##RocheCommandLine=<ID=tool_name,Version="1.0.0",CommandLine="tool_name --input foo.bam">
 * ##Sample=<Name=SAMPLE_A,LibraryId=SBX-D.SAMPLE_A,Platform=SBX,PlatformModel=AXELIOS_1>
 * metric_name  depth_distribution
 * count  6
 * median 1
 * mean   1.17
 * min    0
 * max    3
 * @endcode
 *
 * @param summary_histograms a vector of (name, HistogramSummary) tuples containing pre-computed summary statistics
 * @param output_histogram_summary the path to the output TSV file
 * @param output_percentiles a vector of percentiles to include in the summary (e.g., {25, 50, 75})
 * @param output_cutoffs a vector of cutoffs to include in the summary. If specified, the summary will contain the
 *        percentage of counts for bins that are greater than or equal to the given cutoff.
 * @param output_total_count if true, the total count of all bins will be included in the summary
 * @param tsv_metadata command line info and optional sample information written as the comment header
 */
template <typename T>
  requires Numeric<T>
void WriteSummaryHistogramsToTsv(const HistogramSummaries<T>& summary_histograms,  // NOSONAR(S5536)
                                 const fs::path& output_histogram_summary,
                                 const std::vector<u64>& output_percentiles,
                                 const std::vector<u64>& output_cutoffs,
                                 const bool output_total_count,
                                 const io::TsvMetadata& tsv_metadata) {
  constexpr u8 kPercentagePrecision = 2;
  auto out = std::ofstream{output_histogram_summary};
  auto writer = csv::make_tsv_writer(out).set_auto_flush(false);
  io::WriteTsvMetadata(out, tsv_metadata);

  // Combine all summary histograms by metric
  vec<std::string> header{"metric_name"};
  vec<std::string> count{"count"};
  vec<std::string> median{"median"};
  vec<std::string> mean{"mean"};
  vec<std::string> stddev{"stddev"};
  vec<std::string> min{"min"};
  vec<std::string> max{"max"};
  for (const auto& [summary_histogram_name, summary_histogram_stats] : summary_histograms) {
    header.push_back(summary_histogram_name);
    count.push_back(ValueToString<T>(summary_histogram_stats.count));
    median.push_back(ValueToStringOrNA<u64>(summary_histogram_stats.median));
    mean.push_back(ValueToStringOrNA<f64>(summary_histogram_stats.mean, kMeanPrecision));
    stddev.push_back(ValueToStringOrNA<f64>(summary_histogram_stats.stddev, kStddevPrecision));
    min.push_back(ValueToStringOrNA<u64>(summary_histogram_stats.min));
    max.push_back(ValueToStringOrNA<u64>(summary_histogram_stats.max));
  }

  // Write the header
  writer << header;
  if (output_total_count) {
    writer << count;
  }
  writer << median;
  writer << mean;
  writer << stddev;
  writer << min;
  writer << max;

  // Custom percentiles are added if present
  if (!output_percentiles.empty()) {
    for (const auto& output_percentile : output_percentiles) {
      vec percentile{fmt::format("percentile_{}", output_percentile)};
      for (const auto& [_, summary_histogram_stats] : summary_histograms) {
        if (summary_histogram_stats.percentiles.contains(output_percentile)) {
          percentile.push_back(ValueToStringOrNA<u64>(summary_histogram_stats.percentiles.at(output_percentile)));
        } else {
          percentile.push_back(kNA);
        }
      }
      writer << percentile;
    }
  }

  // Custom percentage cutoffs are added if present
  if (!output_cutoffs.empty()) {
    for (const auto& output_cutoff : output_cutoffs) {
      vec cutoff{fmt::format("percentage_{}", output_cutoff)};
      for (const auto& [_, summary_histogram_stats] : summary_histograms) {
        if (summary_histogram_stats.cutoff_percentage.contains(output_cutoff)) {
          cutoff.push_back(ValueToStringOrNA<f64>(summary_histogram_stats.cutoff_percentage.at(output_cutoff),
                                                  kPercentagePrecision));
        } else {
          cutoff.push_back(kNA);
        }
      }
      writer << cutoff;
    }
  }
}

// Note: remove this functionality when done with refactoring alignment metrics
extern const HistogramSummaryConfig kDefaultCoverageHistogramSummaryConfig;

extern const HistogramSummaryConfig kDefaultReadLengthHistogramSummaryConfig;

template <typename T>
  requires Numeric<T>
// Note: cognitive complexity 32 > 25; legacy function slated for removal (see "Note: remove" above).
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
HistogramSummary<T> CalculateLegacyHistogramSummary(
    const Histogram<T>& histogram,
    const vec<u64>& cutoffs,
    const vec<u64>& percentiles,
    const HistogramSummaryConfig& config = kDefaultCoverageHistogramSummaryConfig) {
  u64 total = 0;
  u64 non_zero_count = 0;
  u64 zero_count = 0;
  u64 min_value = std::numeric_limits<u64>::max();
  u64 max_value = 0;

  Histogram<T> updated_histogram = Histogram<T>(histogram.counts.size(), 0);
  for (u32 value = 0; std::cmp_less(value, updated_histogram.counts.size()); ++value) {
    updated_histogram.At(value) = histogram.At(value);
  }
  // merge outliers
  for (const auto& [key, value] : histogram.outliers) {
    updated_histogram.outliers[key] += value;
  }

  for (u32 value = 0; std::cmp_less(value, updated_histogram.counts.size()); ++value) {
    const u64 count = updated_histogram.counts.at(value);
    total += 1LL * value * count;

    if (value > 0) {
      non_zero_count += count;
    } else {
      zero_count += count;
    }
    // min_value is the first bin with non-zero count in sorted histogram array
    if (min_value == std::numeric_limits<u64>::max() && count != 0) {
      min_value = value;
    }
    // max_value is the last bin with non-zero count in sorted histogram array, unless it is greater than max_counts
    if (value != 0 && count != 0) {
      max_value = value;
    }
  }

  // calculate min, max, total, zero_count, non_zero_count for outliers
  for (const auto& [key, value] : updated_histogram.outliers) {
    total += 1LL * key * value;
    // value is always greater than 0 for outliers, so we can safely add to non_zero_count
    non_zero_count += value;
    if (min_value == std::numeric_limits<u64>::max()) {
      min_value = key;
    }
    // the last outlier is the max
    max_value = key;
  }

  // for computing percentile, ignore the 0 bin
  Histogram<T> updated_histogram_without_zero = Histogram<T>(updated_histogram.counts.size() - 1, 0);
  for (u32 value = 1; std::cmp_less(value, updated_histogram.counts.size()); ++value) {
    updated_histogram_without_zero.At(value - 1) = histogram.At(value);
  }
  updated_histogram_without_zero.outliers = updated_histogram.outliers;

  // determine if the histogram is empty (all bins are zero) to set PercentileMap and PercentageCutoffMap values to 0
  bool is_empty_histogram = updated_histogram_without_zero.IsEmpty();

  PercentileMap percentiles_map;
  if (config.output_percentiles) {
    // these percentiles values are needed to calculate other metrics like ratio_90_to_10 and ratio_50_to_10
    std::set<u32> required_percentile = {kLowerRatioPercentile, kMedianPercentile, kUpperRatioPercentile};
    // merge cli percentile with require percentile
    for (const auto& percentile : percentiles) {
      required_percentile.insert(percentile);
    }

    for (const auto& percentile : required_percentile) {
      // need to add 1 to percentile to get the correct bin index (since the first bin is no longer 0)
      percentiles_map[percentile] =
          is_empty_histogram ? 0 : ComputePercentile(updated_histogram_without_zero, percentile) + 1;
    }
  }

  PercentageCutoffMap cutoff_percentage_map;
  if (!cutoffs.empty() && config.output_cutoff_percentages) {
    for (const auto& cutoff : cutoffs) {
      // need to subtract 1 from cutoff to get the correct bin index (since the first bin is no longer 0)
      cutoff_percentage_map[cutoff] =
          is_empty_histogram ? 0 : ComputeCutoffPercentage(updated_histogram_without_zero, cutoff - 1);
    }
  }

  const f64 mean_value = non_zero_count == 0 ? 0.0 : static_cast<f64>(total) / static_cast<f64>(non_zero_count);
  return HistogramSummary<T>{.count = static_cast<T>(zero_count + non_zero_count),
                             .mean = mean_value,
                             .min = min_value,
                             .max = max_value,
                             .count_zero = zero_count,
                             .count_non_zero = non_zero_count,
                             .percentiles = std::move(percentiles_map),
                             .cutoff_percentage = std::move(cutoff_percentage_map),
                             .config = config};
}

// Note: remove when removing coverage-histogram
template <typename T>
  requires Numeric<T>
void WriteLegacySummaryHistogram(const HistogramSummary<T>& summary_histogram,
                                 const fs::path& output_histogram_summary,
                                 const std::vector<u64>& output_percentiles) {
  auto out = std::ofstream{output_histogram_summary};
  auto writer = csv::make_tsv_writer(out).set_auto_flush(false);
  if (summary_histogram.percentiles.count(kMedianPercentile) &&
      summary_histogram.percentiles.at(kMedianPercentile).has_value()) {
    writer << vec<std::string>{summary_histogram.config.median,
                               std::to_string(summary_histogram.percentiles.at(kMedianPercentile).value())};
  }
  if (summary_histogram.mean.has_value()) {
    writer << vec<std::string>{summary_histogram.config.mean, std::to_string(summary_histogram.mean.value())};
  }
  if (summary_histogram.min.has_value()) {
    writer << vec<std::string>{summary_histogram.config.min, std::to_string(summary_histogram.min.value())};
  }
  if (summary_histogram.max.has_value()) {
    writer << vec<std::string>{summary_histogram.config.max, std::to_string(summary_histogram.max.value())};
  }

  // We write out default ratios if the config is set to true
  if (summary_histogram.config.output_percentile_ratios) {
    auto ratio_90_to_10percentile =
        GetRatioPercentile(summary_histogram.percentiles, kUpperRatioPercentile, kLowerRatioPercentile);
    writer << vec<std::string>{
        "ratio_90_to_10percentile",
        ratio_90_to_10percentile.has_value() ? std::to_string(ratio_90_to_10percentile.value()) : "0"};
    auto ratio_median_to_10percentile =
        GetRatioPercentile(summary_histogram.percentiles, kMedianPercentile, kLowerRatioPercentile);
    writer << vec<std::string>{
        "ratio_median_to_10percentile",
        ratio_median_to_10percentile.has_value() ? std::to_string(ratio_median_to_10percentile.value()) : "0"};
  }

  // Custom percentiles are added if present
  if (!output_percentiles.empty()) {
    for (const auto& key : output_percentiles) {
      if (summary_histogram.percentiles.count(key) && summary_histogram.percentiles.at(key).has_value()) {
        std::string updated_key = fmt::format("percentile_{}", key);
        writer << vec<std::string>{updated_key, std::to_string(summary_histogram.percentiles.at(key).value())};
      }
    }
  }

  if (summary_histogram.config.output_non_zero) {
    writer << vec<std::string>{summary_histogram.config.non_zero, std::to_string(summary_histogram.count_non_zero)};
  }
  if (summary_histogram.config.output_zero) {
    writer << vec<std::string>{summary_histogram.config.zero, std::to_string(summary_histogram.count_zero)};
  }
  if (summary_histogram.config.output_cutoff_percentages) {
    for (const auto& [key, value] : summary_histogram.cutoff_percentage) {
      if (value.has_value()) {
        std::string updated_key = fmt::format("percentage_{}X", key);
        writer << vec<std::string>{updated_key, std::to_string(value.value())};
      }
    }
  }
}

}  // namespace xoos::histogram
