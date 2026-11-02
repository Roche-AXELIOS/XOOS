#pragma once

#include <Eigen/Dense>

#include <xoos/types/int.h>
#include <xoos/types/str-container.h>

#include "genomic-segments.h"

namespace xoos::cnc {
enum class CbsMaxTMethod : u8 {
  kBruteForce,
  kFast,
  kFastMaxN
};
extern const StrMap<CbsMaxTMethod> kStringToCbsMethod;
Segment GetMaxTBreakpoints(const Eigen::VectorXd& obvs, size_t max_n, size_t min_n, CbsMaxTMethod method);
Segment GetMaxTBreakpointsBruteForce(const Eigen::VectorXd& obvs, size_t max_n, size_t min_n);
Segment GetMaxTBreakpointsFast(const Eigen::VectorXd& obvs, size_t max_n, size_t min_n);
Segment GetMaxTBreakpointsFastMaxN(const Eigen::VectorXd& obvs, size_t max_n, size_t min_n);
}  // namespace xoos::cnc
