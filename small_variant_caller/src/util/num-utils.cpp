#include "num-utils.h"

#include <algorithm>

#include <xoos/error/error.h>

namespace xoos::svc {

/**
 * @brief Find the median of a vector of integers.
 * @param vals Vector of 32-bit unsigned integers
 * @return The median value
 */
u32 Median(vec<u32>& vals) {
  if (!vals.empty()) {
    const size_t num_vals = vals.size();
    const int half_size = static_cast<int>(num_vals) / 2;
    // Sort the vector just enough to find the middle value
    std::nth_element(vals.begin(), vals.begin() + half_size, vals.end());
    auto val = vals[half_size];
    if ((num_vals % 2) != 0u) {
      // Odd number of elements, return the middle value
      return val;
    } else {
      // Even number of elements, find the value before the middle value
      // Elements before `half_size` are smaller or equal, but they are not sorted
      auto max_it = std::max_element(vals.begin(), vals.begin() + half_size);
      return (*max_it + val) / 2;
    }
  }
  throw error::Error("Cannot calculate median from an empty vector");
}

}  // namespace xoos::svc
