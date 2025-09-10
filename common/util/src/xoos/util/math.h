#pragma once

#include <concepts>

#include <xoos/types/float.h>
#include <xoos/types/int.h>

namespace xoos::math {

// This is the current maximum phred score associated with the current duplex chemistry for a concordant base
// though this is subject to change in the future
constexpr f64 kMaxPhredScore = 93.0;

/**
 * Subtract unsigned int `y` from unsigned int `x` and saturate the result to
 * the minimum value of the type `T`. That is, if the result underflows,
 * it returns the minimum value of `T`.
 */
template <std::unsigned_integral T>
T SatSub(const T x, const T y) {
  return x > y ? x - y : T{0};
}

/**
 * Divides the numerator by the denominator and rounds the result to the nearest integer.
 */
u32 DivideAndRound(u32 numerator, u32 denominator);

/**
 * Calculates Phred score from the error rate.
 * The Phred score is calculated as -10 * log10(error_rate) given a max Phred score.
 */
f64 ErrorRateToPhred(f64 error_rate, f64 max_phred_score = kMaxPhredScore);

}  // namespace xoos::math
