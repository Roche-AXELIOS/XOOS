#pragma once

#include <algorithm>
#include <iterator>
#include <optional>
#include <ranges>
#include <string>

namespace xoos::util::container {

template <typename C, typename K = C::key_type, typename V = C::value_type::second_type>
std::optional<V> Find(const C& elements, const K& k) {
  auto it = elements.find(k);
  if (it == std::end(elements)) {
    return std::nullopt;
  }
  return it->second;
}

template <typename C, typename T>
bool Contains(const C& elements, const T& value) {
  return std::find(std::cbegin(elements), std::cend(elements), value) != std::cend(elements);
}

/**
 * Find the second largest element in a range.
 *
 * This function finds the iterator to the second largest unique value.
 * For example, in the range [1, 3, 5, 5, 2], the second largest element is 3.
 *
 * @param first The iterator to the beginning of the range.
 * @param last The iterator to the exclusive end of the range.
 *
 * @return An iterator to the second largest element in the range.
 * If the range is empty or has only one unique element, returns last.
 */
template <std::forward_iterator ForwardIt>
ForwardIt FindSecondLargestUnique(ForwardIt first, ForwardIt last) {
  if (first == last || std::next(first) == last) {
    return last;  // Return last if the range is empty or has only one element
  }
  auto largest = first;
  auto second_largest = last;
  // Start searching from the second element
  auto curr = std::next(first);
  for (; curr != last; ++curr) {
    if (*curr > *largest) {
      second_largest = largest;
      largest = curr;
    } else if (*curr < *largest) {
      // Either we haven't found a second_largest yet, or curr is larger than the current second_largest
      if (second_largest == last || (*curr > *second_largest)) {
        second_largest = curr;
      }
    }
    // If *curr == *largest, we do nothing to ensure uniqueness
  }
  return second_largest;
}

/**
 * @brief permute a std::vector in-place according a user-defined order
 * @tparam T type of elements in the std::vector
 * @tparam U (integral) type of the elements in argsort
 * @param vec std::vector to permute
 * @param argsort new permutation (order of indices)
 * @return reference to the permuted std::vector
 */
template <std::ranges::random_access_range T, std::ranges::random_access_range U>
  requires std::integral<typename U::value_type>
T& Permute(T& vec, const U& argsort) {
  if (argsort.size() != vec.size()) {
    throw std::runtime_error("Permute: vec and argsort must be the same size!");
  }
  T nvec(vec.size());
  for (size_t i = 0; i < vec.size(); ++i) {
    nvec[i] = vec[argsort[i]];
  }
  vec = std::move(nvec);
  return vec;
}

}  // namespace xoos::util::container
