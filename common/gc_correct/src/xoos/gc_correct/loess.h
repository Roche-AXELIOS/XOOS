#pragma once

#include <Eigen/Dense>
#include <cstddef>
#include <tuple>

#include <xoos/types/float.h>
#include <xoos/types/int.h>

namespace xoos::gc_correct {
// Performs LOESS (Locally Estimated Scatterplot Smoothing) regression to smooth noisy data using local polynomial
// fitting.
Eigen::VectorXd Loess(const Eigen::VectorXd& x_to_predict,
                      const Eigen::VectorXd& x,
                      const Eigen::VectorXd& y,
                      f64 first_span,
                      s32 degree);
// Fits a weighted polynomial to data points and evaluates it at new x locations for local smoothing approximation.
Eigen::VectorXd WeightedPolyFit(const Eigen::VectorXd& x,
                                const Eigen::VectorXd& y,
                                const Eigen::VectorXd& new_x,
                                const Eigen::VectorXd& weights,
                                size_t degree);

// Computes tricubic distance-based weights to emphasize data points closer to a target location in local regression.
Eigen::VectorXd TricubicDistWeights(const Eigen::VectorXd& v, f64 x);

// Identifies the indices of the n nearest data points to a target value for efficient local window selection.
std::tuple<Eigen::Index, Eigen::Index> GetNClosestRange(const Eigen::VectorXd& arr, Eigen::Index n, f64 x);
}  // namespace xoos::gc_correct
