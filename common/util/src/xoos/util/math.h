#pragma once

#include <cmath>
#include <concepts>
#include <limits>

#include <xoos/types/float.h>
#include <xoos/types/int.h>

#include "xoos/types/vec.h"

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

/**
 * @brief Partially sort a non-empty vector of integers and find the median.
 * @details If the vector has an odd number of elements, return the middle element.
 * If the vector has an even number of elements, return the floor of the average of the two middle elements.
 * @param vals Vector of 32-bit unsigned integers
 * @return The median value
 */
u32 Median(vec<u32>& vals);

/**
 * @brief Determine if a floating point number is close to zero, within a tolerance defined as 100 times the machine
 * epsilon for the type.
 * @tparam T Floating point type
 * @param f Floating point number to check
 * @return True if the number is close to zero, false otherwise
 */
template <std::floating_point T>
bool IsCloseToZero(T f) {
  auto epsilon = std::numeric_limits<T>::epsilon() * 100;
  return std::fabs(f) < epsilon;
}

/**
 * @brief Determine if two floating point numbers are equal within a specified tolerance.
 * @tparam T Floating point type
 * @param f1 First floating point number
 * @param f2 Second floating point number
 * @param tolerance Tolerance for comparison (default is 1e-6)
 * @return True if the numbers are equal within the specified tolerance, false otherwise
 */
template <std::floating_point T>
bool IsEqualWithTolerance(T f1, T f2, T tolerance = T(1e-6)) {
  return std::abs(f1 - f2) < tolerance;
}

}  // namespace xoos::math
