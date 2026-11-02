#pragma once
#include <Eigen/Dense>

#include <xoos/types/float.h>

#include "util/eigen-math.h"

namespace xoos::cnc {
f64 TStat(f64 n_v1, f64 mean_v1, f64 var_v1, f64 n_v2, f64 mean_v2, f64 var_v2);
f64 TStatFromVecs(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2);
f64 Variance(f64 m, f64 mean, f64 sq_sum);
Eigen::VectorXd& FisherYatesShuffleInPlace(Eigen::VectorXd& arr, eigen::ArmaRng& rng);
f64 TStatSquaredMean(const Eigen::VectorXd& psums, size_t i, size_t j);
f64 TStatSquaredMean(f64 total_sum, f64 total_length, f64 minor_sum, f64 minor_length);
f64 TStatSquaredMean(f64 total_sum, size_t total_length, f64 minor_sum, size_t minor_length);
}  // namespace xoos::cnc
