#pragma once

#include <algorithm>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>

#include <fmt/format.h>

#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::histogram {

// templating with the concept Numeric to ensure that T is a numeric type (int, float, f64, etc.)
template <typename T>
concept Numeric = std::is_arithmetic_v<T>;

// A vector representing a histogram where the index is the bin and the value at a given index is the count
// We use T (int, float, f64, etc.) to represent the count (or frequency) of a bin
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
  explicit Histogram(const size_t num_bins = 0,  // NOSONAR(S5536,S1712) template instantiated in other TUs
                     T initial_value = 0)        // NOSONAR(S1712)
      : counts(num_bins, initial_value) {
  }

  // Accessor to get the bin count at a specific index
  size_t Size() const {  // NOSONAR(S5536) template instantiated in other TUs
    // return the total number of bins including outliers
    if (outliers.empty()) {
      return counts.size();
    }
    return std::max<size_t>(counts.size(), outliers.rbegin()->first + 1);
  }

  // Accessor to get the bin count at a specific index
  T& At(u64 index) {  // NOSONAR(S5536) template instantiated in other TUs
    if (index < counts.size()) {
      return counts.at(index);
    }
    if (outliers.contains(index)) {
      return outliers.at(index);
    }
    throw std::out_of_range(fmt::format("Histogram bin index is out of bounds: {}", index));
  }

  const T& At(u64 index) const {  // NOSONAR(S5536) template instantiated in other TUs
    if (index < counts.size()) {
      return counts.at(index);
    }
    if (outliers.contains(index)) {
      return outliers.at(index);
    }
    throw std::out_of_range(fmt::format("Histogram bin index is out of bounds: {}", index));
  }

  // Member operator kept for symmetry with the class template; instantiated in other TUs
  Histogram operator+(const Histogram& other) const {  // NOSONAR(S5536,S2807,M23_365)
    if (counts.size() != other.counts.size()) {
      throw std::invalid_argument("Histograms must be the same size for + operator");
    }
    Histogram result{other.counts.size()};
    for (size_t i = 0; i < other.counts.size(); ++i) {
      result.counts[i] = counts[i] + other.counts[i];
    }
    for (const auto& [key, value] : outliers) {
      result.outliers[key] += value;
    }
    for (const auto& [key, value] : other.outliers) {
      result.outliers[key] += value;
    }
    return result;
  }

  Histogram& operator+=(const Histogram& other) {  // NOSONAR(S5536) template instantiated in other TUs
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

  bool IsEmpty(const u64 start_index = 0) const {  // NOSONAR(S5536,S1712) template instantiated in other TUs
    const auto begin = start_index < counts.size() ? counts.begin() + ToSigned(start_index) : counts.end();
    return std::ranges::all_of(begin, counts.end(), [](const T& bin_count) { return bin_count == 0; }) &&
           outliers.empty();
  }

  // Resets values to default, which will be the same size.
  // We have to assume the default value of 0 as we have no way to know the initial value.
  void Reset() {  // NOSONAR(S5536) template instantiated in other TUs
    outliers.clear();
    // set all counts to 0
    std::ranges::fill(counts, T{0});
  }

  /**
   * Add a count to a histogram at the specified bin.
   * If it exceeds the max number of bins, it adds the count to the outliers
   */
  void AddCountToHistogram(u64 bin, T count) {  // NOSONAR(S5536) template instantiated in other TUs
    if (bin < counts.size()) {
      counts.at(bin) += count;
    } else {
      outliers[bin] += count;
    }
  }
};

// A vector representing histograms for distributions that share the same binning
template <typename T>
  requires Numeric<T>
using Histograms = vec<std::tuple<std::string, Histogram<T>>>;

// We sum the histograms with identical names.
// All histograms are assumed to have the same default size.
template <typename T>
  requires Numeric<T>
Histograms<T>& operator+=(Histograms<T>& lhs,  // NOSONAR(S5536) template instantiated in other TUs
                          const Histograms<T>& rhs) {
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

}  // namespace xoos::histogram
