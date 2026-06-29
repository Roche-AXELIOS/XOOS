#pragma once

#include <algorithm>
#include <filesystem>
#include <limits>
#include <map>
#include <optional>
#include <ranges>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <fmt/format.h>

#include <csv.hpp>

#include <xoos/io/metadata-util.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::histogram {

namespace fs = std::filesystem;

// templating with the concept Numeric to ensure that T is a numeric type (int, float, double, etc.)
template <typename T>
concept Numeric = std::is_arithmetic_v<T>;

const std::string kNA = "NA";
const vec<u64> kDefaultPercentiles{10, 25, 50, 75, 90};

// HistogramBinOutput defines how histogram bins should be outputted when writing to a file
enum HistogramBinOutput {
  kDefault,       // Default behavior: output all bins individually
  kOmitFirstBin,  // Omit the first bin
  kMaxLastBin,  // Add "+" to the last bin (when the last bin contains all values greater than or equal to the last bin)
  kOmitFirstBinAndMaxLastBin,            // Omit first bin AND add "+" to the last bin
  kOmitFirstBinAndMaxLastBinWithOutlier  // Omit first bin AND the last bin is the outliers which is size wth "+" to the
                                         // last bin
};

// key is percentile, value is the count
using PercentileMap = std::map<u64, std::optional<u64>>;

// key is cutoff, value is the percentage of counts for bins that are greater than or equal to the given cutoff
using PercentageCutoffMap = std::map<u64, std::optional<double>>;

// A vector representing a histogram where the index is the bin and the value at a given index is the count
// We use T (int, float, double, etc.) to represent the count (or frequency) of a bin
// Bins are always discrete values (0, 1, 2, ...), so we use u64 to represent the bin index
template <typename T>
  requires Numeric<T>
struct Histogram {
  // counts that includes the 0-vector up to the 1 - number of max bins declared
  vec<T> counts;
  // we assume that an outlier value can never be 0 and that the outlier key is always greater than the last bin index
  // outliers contain all values that are greater or equal to the number of bins, including the maximum bin size
  std::map<u64, T> outliers{};

  // Constructor to initialize the histogram with a given number of bins, inferred as the maximum bin size recorded.
  // The counts are the number of bins, which includes the 0-vector, but not the maximum bins, which are stored in the
  // outliers
  explicit Histogram(const size_t num_bins = 0, T initial_value = 0) : counts(num_bins, initial_value) {
  }

  // Accessor to get the bin count at a specific index
  size_t Size() const {
    // return the last index of the outliers map
    if (outliers.empty()) {
      return counts.size();
    }
    return std::max<size_t>(counts.size(), outliers.rbegin()->first);
  }

  // Accessor to get the bin count at a specific index
  T& At(u64 index) {
    if (index < counts.size()) {
      return counts.at(index);
    }
    if (outliers.contains(index)) {
      return outliers.at(index);
    }
    throw std::out_of_range("Histogram index bin is out of bounds");
  }

  const T& At(u64 index) const {
    if (index < counts.size()) {
      return counts.at(index);
    }
    if (outliers.contains(index)) {
      return outliers.at(index);
    }
    throw std::out_of_range("Histogram index is out of bounds");
  }

  Histogram operator+(const Histogram& other) const {
    if (counts.size() != other.counts.size()) {
      throw std::invalid_argument("Histograms must be the same size for + operator");
    }
    Histogram result{other.counts.size()};
    for (size_t i = 0; i < other.counts.size(); ++i) {
      result.counts[i] = counts[i] + other.counts[i];
    }
    for (const auto& [key, value] : other.outliers) {
      result.outliers[key] = outliers.at(key) + value;
    }
    return result;
  }

  Histogram& operator+=(const Histogram& other) {
    if (counts.size() != other.counts.size()) {
      throw std::invalid_argument("Histograms must be the same size for += operator");
    }
    for (size_t i = 0; i < counts.size(); ++i) {
      counts[i] += other.counts[i];
    }
    for (const auto& [key, value] : other.outliers) {
      outliers[key] += value;
    }
    return *this;
  }

  [[nodiscard]] bool IsEmpty(u64 start_index = 0) const {
    // check if all values are 0
    return std::ranges::all_of(counts.begin() + ToSigned(start_index),
                               counts.end(),
                               [](const T& bin_count) { return bin_count == 0; }) &&
           outliers.empty();
  }

  // Resets values to default, which will be the same size.
  // We have to assume the default value of 0 as we have no way to know the initial value.
  void Reset() {
    outliers.clear();
    // set all counts to 0
    std::fill(counts.begin(), counts.end(), T{0});
  }

  /**
   * Add a count to a histogram at the specified bin.
   * If it exceeds the max number of bins, it adds the count to the outliers
   */
  void AddCountToHistogram(u64 bin, T count) {
    if (bin < counts.size()) {
      counts.at(bin) += count;
    } else {
      outliers[bin] += count;
    }
  }
};

// TODO: remove this functionality when done with refactoring alignment metrics
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
  std::optional<u64> mean;
  std::optional<f64> stddev;
  std::optional<u64> min;
  std::optional<u64> max;

  // TODO: remove this functionality when done with refactoring alignment metrics
  u64 count_zero = 0;
  u64 count_non_zero = 0;

  std::map<u64, std::optional<u64>> percentiles;
  std::map<u64, std::optional<double>> cutoff_percentage;

  // TODO: remove this functionality when done with refactoring alignment metrics
  HistogramSummaryConfig config{};
};

// A vector representing histograms for distributions that share the same binning
template <typename T>
  requires Numeric<T>
using Histograms = vec<std::tuple<std::string, Histogram<T>>>;

// We sum the histograms with identical names.
// All histograms are assumed to have the same default size.
template <typename T>
  requires Numeric<T>
Histograms<T>& operator+=(Histograms<T>& lhs, const Histograms<T>& rhs) {
  if (lhs.size() != rhs.size()) {
    throw std::invalid_argument("Histograms must be the same size for += operator");
  }
  for (auto& [lhs_name, lhs_histogram] : lhs) {
    for (const auto& [rhs_name, rhs_histogram] : rhs) {
      if (lhs_name == rhs_name) {
        lhs_histogram += rhs_histogram;
      }
    }
  }
  return lhs;
}

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
    variance_sum += bin_count * std::pow(static_cast<f64>(bin_index) - mean, 2);
  }
  // count outliers
  for (const auto& [bin_index, bin_count] : histogram.outliers) {
    variance_sum += bin_count * std::pow(static_cast<f64>(bin_index) - mean, 2);
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
u64 ComputeMean(const Histogram<T>& histogram, u64 start_index = 0) {
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
  // Added std::round to avoid integer truncation
  return static_cast<u64>(std::round(static_cast<double>(weighted_sum) / static_cast<double>(total_count)));
}

/**
 * Compute the min bin of the histogram. The min bin is the first bin with non-zero count.
 * If there is no bin with non-zero count, it returns the first outlier bin.
 * If there is no outlier bin, it returns the last bin index (which is the size of the histogram counts).
 * NOTE: The histogram is assumed to be non-empty (at least one bin has a non-zero count)
 */
template <typename T>
  requires Numeric<T>
u64 ComputeMin(const Histogram<T>& histogram) {
  if (histogram.IsEmpty()) {
    throw std::invalid_argument("Cannot compute min of an empty histogram");
  }
  // find the first bin with non-zero count
  auto it = std::ranges::find_if(histogram.counts, [](const u64 bin_count) { return bin_count != 0; });
  auto min = std::ranges::distance(histogram.counts.begin(), it);

  if (min < histogram.counts.size()) {
    // if the min is within the histogram counts, we return it
    return static_cast<u64>(min);
  }

  // NOTE: unlikely, but if the histogram has outliers and the min is histogram.size(), we need to check if the first
  // outlier is an outlier
  if (!histogram.outliers.empty()) {
    return histogram.outliers.begin()->first;
  }
  return static_cast<u64>(histogram.counts.size());
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
                                          [](const u64 bin_count) { return bin_count != 0; });
  auto max = std::ranges::distance(it_reversed, histogram.counts.rend() - 1);
  // if the histogram has outliers, we need to check if the last bin is an outlier
  if (!histogram.outliers.empty() && histogram.outliers.rbegin()->first >= histogram.outliers.size()) {
    // we return the last outlier bin minus 1 since the outlier bin is not included in the histogram
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
    throw std::invalid_argument("Cannot compute percentily of an empty histogram");
  }
  if (percentile < 0 || percentile > 100) {
    throw std::out_of_range("Percentile must be between 0 and 100");
  }
  // create new histogram with the starting index
  // this is useful when we want to compute percentiles for histograms that start at a specific index
  // e.g. when we want to compute percentiles for histograms that start at 2, we need to create a new histogram
  // with the starting index of 2
  histogram::Histogram<u64> updated_histogram_with_new_start_index;
  updated_histogram_with_new_start_index = histogram::Histogram<u64>(histogram.counts.size() - start_index, 0);
  for (u32 bin_index = start_index; std::cmp_less(bin_index, histogram.counts.size()); ++bin_index) {
    updated_histogram_with_new_start_index.At(bin_index - start_index) = histogram.At(bin_index);
  }
  // add outliers to the new histogram
  for (const auto& [bin_index, bin_count] : histogram.outliers) {
    updated_histogram_with_new_start_index.outliers[bin_index - start_index] = bin_count;
  }

  u64 percentile_bin = 0;
  if (percentile == 0) {
    percentile_bin = ComputeMin(updated_histogram_with_new_start_index);
  } else if (percentile == 100) {
    percentile_bin = ComputeMax(updated_histogram_with_new_start_index);
  } else {
    T total_count = ComputeCount(updated_histogram_with_new_start_index);
    T cumulative_count = 0;
    for (u64 bin_index = 0; std::cmp_less(bin_index, updated_histogram_with_new_start_index.counts.size());
         ++bin_index) {
      cumulative_count += updated_histogram_with_new_start_index.At(bin_index);
      // we divide by 100.0 to get the ratio instead of multiplying by 100 in order to prevent potential overflow
      if (cumulative_count >= total_count * percentile / 100.0) {
        return bin_index + start_index;  // add the start index to the percentile bin to get the correct bin index
      }
    }
    // handle outliers
    for (const auto& [bin_index, bin_count] : updated_histogram_with_new_start_index.outliers) {
      cumulative_count += bin_count;
      if (cumulative_count >= total_count * percentile / 100.0) {
        return bin_index + start_index;  // add the start index to the percentile bin to get the correct bin index
      }
    }
  }
  return percentile_bin + start_index;  // add the start index to the percentile bin to get the correct bin index
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
double ComputeCutoffPercentage(const Histogram<T>& histogram, const u64 cutoff) {
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
  // calculate percentage
  const auto percentage_above_cutoff =
      static_cast<double>(total_count_above_cutoff * 100) / static_cast<double>(total_count);
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

  const u64 median = exclude_zero_bins ? ComputePercentile(histogram, 50, 1) : ComputePercentile(histogram, 50);

  const u64 mean = ComputeMean(histogram, exclude_zero_bins ? 1 : 0);  // mean is computed from the original histogram
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
std::optional<double> GetRatioPercentile(const PercentileMap& percentiles, u64 numerator, u64 denominator);

// helper function to get string with 2 decimal places
template <typename T>
  requires Numeric<T>
std::string ValueToString(const T& val, const u8 precision = 0) {
  std::stringstream ss;
  ss << std::fixed << std::setprecision(precision) << static_cast<double>(val);
  return ss.str();
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
 *
 * Writes the histogram(s) to a file in tab-separated format. Each row is a bin and each column is a value from the
 * respective histogram (aside from the first column, which is for the bins themselves). This provides the user with a
 * single file for multiple distributions (if applicable). Precisions of each column can be specified with
 * histogram_precisions.
 *
 * NOTE: All histograms must share the same x-axis (bins). Additional formatting can be specified with
 * histogram_bin_output (i.e. omit the first bin, add a "+" to the last bin).
 *
 * Step 1: Open the output file and compute the number of bins (rows).
 * Step 2: Write the header.
 * Step 3: For each bin, combine all histogram data into a suitable format (vec<std::string>) and write to output.
 *
 * Example:
 * cluster_size  total_clusters
 * 0  0
 * 1  5
 * 2  3
 * 3+  1
 *
 */
template <typename T>
  requires Numeric<T>
void WriteHistograms(const Histograms<T>& histograms,
                     const fs::path& output_histogram,
                     const std::string& binning_name,
                     const HistogramBinOutput histogram_bin_output = kDefault,
                     const vec<u8>& histogram_precisions = {},
                     const io::Comments& comments = {}) {
  // if the histogram does not contain any bins, there is nothing to write
  if (histograms.empty()) {
    return;
  }
  auto out = std::ofstream{output_histogram};
  auto writer = csv::make_tsv_writer_buffered(out);
  io::WriteTsvComments(writer, comments);
  auto num_bins = std::get<1>(histograms.at(0)).counts.size();

  // Write the header
  vec header{binning_name};
  for (const auto& [histogram_name, histogram] : histograms) {
    header.push_back(histogram_name);
  }
  writer << header;

  // Write the histogram data
  for (u64 bin_index = 0; std::cmp_less(bin_index, num_bins); ++bin_index) {
    if ((histogram_bin_output == kOmitFirstBin || histogram_bin_output == kOmitFirstBinAndMaxLastBinWithOutlier ||
         histogram_bin_output == kOmitFirstBinAndMaxLastBin) &&
        bin_index == 0) {
      continue;
    }
    auto bin = std::to_string(bin_index);
    // if the last bin contains counts from greater bins, add a "+" to the bin name
    if (bin_index == num_bins - 1 &&
        (histogram_bin_output == kMaxLastBin || histogram_bin_output == kOmitFirstBinAndMaxLastBin)) {
      bin += "+";
    }
    vec row{bin};
    size_t precision_index = 0;
    for (const auto& [_, histogram] : histograms) {
      // here we handle the last row differently if the histogram_bin_output is set to kMaxLastBin or
      // kOmitFirstBinAndMaxLastBin
      auto histogram_value = histogram.At(bin_index);
      if (bin_index == num_bins - 1 &&
          (histogram_bin_output == kMaxLastBin || histogram_bin_output == kOmitFirstBinAndMaxLastBin)) {
        // accumulate all values of the outliers
        for (const auto& [outlier_bin_index, outlier_value] : histogram.outliers) {
          if (outlier_bin_index >= bin_index) {
            histogram_value += outlier_value;
          }
        }
      }
      if (histogram_precisions.empty()) {
        row.push_back(ValueToString<T>(histogram_value));
      } else {
        auto precision = histogram_precisions[precision_index++];
        row.push_back(ValueToString<T>(histogram_value, precision));
      }
    }
    writer << row;
  }
  // if histogram_bin_output is set to default or omit first bin, we need to handle outliers in each histogram
  // separately
  if (histogram_bin_output == kDefault || histogram_bin_output == kOmitFirstBin) {
    // since outliers are not going to necessarily sync up between histograms, we need to get an ordered set of both
    std::set<u64> outlier_bins;
    u64 max_outlier_bin_index = 0;
    for (const auto& [_, histogram] : histograms) {
      for (const auto& [outlier_bin_index, _] : histogram.outliers) {
        outlier_bins.insert(outlier_bin_index);
        max_outlier_bin_index = std::max(max_outlier_bin_index, outlier_bin_index);
      }
    }

    // now we start and the bin_index from the last bin in the histograms.count
    for (u64 bin_index = num_bins; bin_index < max_outlier_bin_index + 1; ++bin_index) {
      // for each bin, we need to check if it is an outlier bin in each histogram
      vec<std::string> row{};
      row.push_back(ValueToString<u64>(bin_index));
      size_t precision_index = 0;
      for (const auto& [_, histogram] : histograms) {
        // if the histogram contains the bin_index, we can use it directly, otherwise we set it to 0
        // find the bin_index in the histogram
        T bin_value = T{0};
        if (histogram.outliers.contains(bin_index)) {
          bin_value = histogram.outliers.at(bin_index);
        }

        if (histogram_precisions.empty()) {
          row.push_back(ValueToString<T>(bin_value));
        } else {
          auto precision = histogram_precisions[precision_index++];
          row.push_back(ValueToString<T>(bin_value, precision));
        }
      }
      writer << row;
    }
  }

  // if we want to treat the very last bin for only outliers, we need to handle it separately
  if (histogram_bin_output == kOmitFirstBinAndMaxLastBinWithOutlier) {
    auto bin = std::to_string(num_bins) + "+";
    vec row{bin};
    size_t precision_index = 0;
    for (const auto& [_, histogram] : histograms) {
      auto histogram_value = T{0};
      // accumulate all values of the outliers
      for (const auto& [outlier_bin_index, outlier_value] : histogram.outliers) {
        histogram_value += outlier_value;
      }
      if (histogram_precisions.empty()) {
        row.push_back(ValueToString<T>(histogram_value));
      } else {
        auto precision = histogram_precisions[precision_index++];
        row.push_back(ValueToString<T>(histogram_value, precision));
      }
    }
    writer << row;
  }
}

constexpr u32 kStddevPrecision = 2;

/**
 *
 * Writes the summary of the histogram(s) to a file in tab-separated format. Each row is a metric and each column is a
 * histogram (aside from the first column, which is for the metric names). This provides the user with a single, summary
 * report for multiple distributions (if applicable).
 *
 * NOTE: All histograms must share the same x-axis (bins). If the histogram is empty (all bins are zero), then some
 * statistics will be represented by the constant kNA (i.e. not applicable). Statistics are represented by 2 decimal
 * places.
 *
 * Step 1: Open the output file.
 * Step 2: Combine all summary histogram data for the header, count, median, mean, min, and max into a suitable format
 * (vec<std::string>) for writing to output.
 * Step 3: Write the header and summary stats to the output file.
 * Step 4: If output percentiles and/or output cutoffs are provided, add each one as a new row to the output file for
 * all histograms.
 *
 * Example:
 * metric_name  depth_distribution
 * count  6
 * median 1
 * mean 1
 * min  0
 * max  3
 *
 * @param summary_histograms a vector of histogram summaries where each tuple contains the name of the histogram and the
 * summary statistics
 * @param output_histogram_summary the path to the output file where the summary will be written
 * @param output_percentiles a vector of percentiles to output in the summary (default is `kDefaultPercentiles`)
 * @param output_cutoffs a vector of cutoffs to output in the summary (default is empty, i.e. no cutoffs). If specified,
 * the summary will contain the percentage of counts for bins that are greater than or equal to the given cutoff.
 * @param output_total_count if true, the total count of all bins will be included in the summary (default is true)
 * @param comments optional comments to include in the output file (default is empty, i.e. no comments)
 *
 */
template <typename T>
  requires Numeric<T>
void WriteSummaryHistograms(const HistogramSummaries<T>& summary_histograms,
                            const fs::path& output_histogram_summary,
                            const std::vector<u64>& output_percentiles = kDefaultPercentiles,
                            const std::vector<u64>& output_cutoffs = {},
                            bool output_total_count = true,
                            const io::Comments& comments = {}) {
  auto out = std::ofstream{output_histogram_summary};
  auto writer = csv::make_tsv_writer_buffered(out);
  io::WriteTsvComments(writer, comments);

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
    median.push_back(ValueToStringOrNA<f64>(summary_histogram_stats.median));
    mean.push_back(ValueToStringOrNA<f64>(summary_histogram_stats.mean));
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
      vec percentile{"percentile_" + std::to_string(output_percentile)};
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
      vec cutoff{"percentage_" + std::to_string(output_cutoff)};
      for (const auto& [_, summary_histogram_stats] : summary_histograms) {
        // Set the output double precision to 2 since these are percentages
        if (summary_histogram_stats.cutoff_percentage.contains(output_cutoff)) {
          cutoff.push_back(ValueToStringOrNA<double>(summary_histogram_stats.cutoff_percentage.at(output_cutoff), 2));
        } else {
          cutoff.push_back(kNA);
        }
      }
      writer << cutoff;
    }
  }
}

// TODO: remove this functionality when done with refactoring alignment metrics
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

template <typename T>
  requires Numeric<T>
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
    max_value = key;  // the last outlier is the max
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
    std::set<u32> required_percentile = {10, 50, 90};
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

  // converting u64 (total) to u32 - which should be theoretical max value of histogram
  const u64 mean_value = non_zero_count == 0 ? 0 : ToSigned(total / non_zero_count);
  return HistogramSummary<T>{.mean = mean_value,
                             .min = min_value,
                             .max = max_value,
                             .count_zero = zero_count,
                             .count_non_zero = non_zero_count,
                             .percentiles = std::move(percentiles_map),
                             .cutoff_percentage = std::move(cutoff_percentage_map),
                             .config = config};
}

// TODO: remove when removing coverage-histogram
template <typename T>
  requires Numeric<T>
void WriteLegacySummaryHistogram(const HistogramSummary<T>& summary_histogram,
                                 const fs::path& output_histogram_summary,
                                 const std::vector<u64>& output_percentiles) {
  auto out = std::ofstream{output_histogram_summary};
  auto writer = csv::make_tsv_writer_buffered(out);
  writer << std::make_tuple(summary_histogram.config.median, summary_histogram.percentiles.at(50).value());
  writer << std::make_tuple(summary_histogram.config.mean, summary_histogram.mean.value());
  writer << std::make_tuple(summary_histogram.config.min, summary_histogram.min.value());
  writer << std::make_tuple(summary_histogram.config.max, summary_histogram.max.value());

  // We write out default ratios if the config is set to true
  if (summary_histogram.config.output_percentile_ratios) {
    auto ratio_90_to_10percentile = GetRatioPercentile(summary_histogram.percentiles, 90, 10);
    if (!ratio_90_to_10percentile.has_value()) {
      writer << std::make_tuple("ratio_90_to_10percentile", 0);
    } else {
      writer << std::make_tuple("ratio_90_to_10percentile", ratio_90_to_10percentile.value());
    }
    auto ratio_median_to_10percentile = GetRatioPercentile(summary_histogram.percentiles, 50, 10);
    if (!ratio_median_to_10percentile.has_value()) {
      writer << std::make_tuple("ratio_median_to_10percentile", 0);
    } else {
      writer << std::make_tuple("ratio_median_to_10percentile", ratio_median_to_10percentile.value());
    }
  }

  // Custom percentiles are added if present
  if (!output_percentiles.empty()) {
    for (const auto& key : output_percentiles) {
      auto value = summary_histogram.percentiles.at(key).value();
      std::string updated_key = fmt::format("percentile_{}", key);
      writer << std::make_tuple(updated_key, value);
    }
  }

  if (summary_histogram.config.output_non_zero) {
    writer << std::make_tuple(summary_histogram.config.non_zero, summary_histogram.count_non_zero);
  }
  if (summary_histogram.config.output_zero) {
    writer << std::make_tuple(summary_histogram.config.zero, summary_histogram.count_zero);
  }
  if (summary_histogram.config.output_cutoff_percentages) {
    for (const auto& [key, value] : summary_histogram.cutoff_percentage) {
      std::string updated_key = fmt::format("percentage_{}X", key);
      writer << std::make_tuple(updated_key, value.value());
    }
  }
}

}  // namespace xoos::histogram
