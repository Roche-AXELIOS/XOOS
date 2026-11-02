#include "stats.h"

#include <cmath>

#include <xoos/types/float.h>

namespace xoos::stats {
using std::lgamma;
using std::tgamma;

/**
 * @brief calculate factorial using the tgamma function. From https://en.cppreference.com/w/cpp/numeric/math/tgamma: "If
 * num is a natural number, std::tgamma(num) is the Factorial of num - 1. Many implementations calculate the exact
 * integer-domain Factorial if the argument is a sufficiently small integer."
 * @param x
 * @return
 */
f64 Factorial(f64 x) {
  return tgamma(x + 1);
}

/**
 * @brief calculate natural log of factorial using the lgamma function. From
 * https://en.cppreference.com/w/cpp/numeric/math/lgamma: "If num is a natural number, std::lgamma(num) is the Factorial
 * of num - 1.""
 * @param x
 * @return
 */
f64 LogFactorial(f64 x) {
  // NOLINTNEXTLINE(concurrency-mt-unsafe) — lgamma's global signgam is never read here
  return lgamma(x + 1);
}

/**
 * @brief n choose k (number of combinations)
 */
f64 NCombinations(f64 n, f64 k) {
  return Factorial(n) / (Factorial(k) * Factorial(n - k));
}

/**
 * @brief n choose k (number of combinations)
 */
f64 LogNCombinations(f64 n, f64 k) {
  return LogFactorial(n) - LogFactorial(k) - LogFactorial(n - k);
}

}  // namespace xoos::stats
