#include "xoos/util/math.h"

#include <cmath>
#include <stdexcept>

#include <xoos/types/float.h>

namespace xoos::math {

static const f64 kMinErrorRate = 0.0;
static const f64 kMaxErrorRate = 1.0;

/**
 * Divides the numerator by the denominator and rounds the result to the nearest integer.
 * If the division results in a fractional value, it rounds up if the fractional part is
 * greater than or equal to 0.5, otherwise it rounds down.
 *
 * Note: If the denominator is zero, this function will cause a runtime error due to
 * division by zero, which is undefined behavior in C++.
 *
 * @param numerator The dividend in the division operation.
 * @param denominator The divisor in the division operation. Must not be zero.
 * @return The result of the division rounded to the nearest integer.
 */
u32 DivideAndRound(const u32 numerator, const u32 denominator) {
  if (denominator == 0) {
    throw std::runtime_error("Denominator cannot be zero in DivideAndRound.");
  }
  // Use static_cast to ensure proper rounding
  return static_cast<u32>(std::round(static_cast<f64>(numerator) / denominator));
}

/**
 *   We calculate the Phred score from the error rate using the formula:
 *   Phred score = -10 * log10(error_rate)
 *
 *   However, we normalize the maximum score to 93, which is the maximum Phred score
 *   that can be represented in the current implementation, though this may change.
 *
 *   For an error near 0, we limit the Phred score to a maximum phred score used for the given chemistry.
 */
f64 ErrorRateToPhred(const f64 error_rate, const f64 max_phred_score) {
  if (error_rate == kMinErrorRate) {
    // If error rate is 0 (the minimum), Phred score is the maximum phred score.
    return max_phred_score;
  }
  if (error_rate < kMinErrorRate || error_rate > kMaxErrorRate) {
    throw std::invalid_argument("Error rate must be in the range [0, 1].");
  }
  return std::fmin(max_phred_score, -10.0 * std::log10(error_rate));
}

}  // namespace xoos::math
