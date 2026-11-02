#pragma once
#include <Eigen/Dense>

#include <xoos/types/float.h>
#include <xoos/types/vec.h>

namespace xoos::cnc {

struct KernelDensityEstimateRet {
  Eigen::VectorXd x;
  Eigen::VectorXd y;
};

KernelDensityEstimateRet GaussianKernelDensityEstimate(const Eigen::VectorXd& x, const Eigen::VectorXd& weights);
KernelDensityEstimateRet GaussianKernelDensityEstimate(const Eigen::VectorXd& x,
                                                       const Eigen::VectorXd& weights,
                                                       f64 bandwidth_multiplier);
Eigen::VectorXd GetBAFFromAD(const Eigen::VectorXd& ref_ads, const Eigen::VectorXd& alt_ads);
vec<size_t> FindLocalMaxima(const Eigen::VectorXd& x);
f64 GetPeakOfBAFDistributionThenMirror(const Eigen::VectorXd& ref_ads, const Eigen::VectorXd& alt_ads);
f64 WeightedMBAFSD(const Eigen::VectorXd& ref_ads, const Eigen::VectorXd& alt_ads, f64 mean_baf);
f64 GetBandwidthSilvermanRuleOfThumb(const Eigen::VectorXd& x);
}  // namespace xoos::cnc
