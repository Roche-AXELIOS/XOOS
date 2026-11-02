#include "coverage/coverage-check.h"

#include <xoos/types/float.h>
#include <xoos/types/int.h>

#include "coverage/coverage.h"
#include "util/eigen-math.h"

namespace xoos::cnc {

CoverageCheckResult CheckMedianCoverage(const CoverageRecords& coverage, const f64 threshold) {
  CoverageCheckResult result{};
  if (coverage.Empty()) {
    result.median_coverage = 0.0;
    result.pct_windows_below_threshold = 100.0;
    result.is_sufficient = false;
    return result;
  }
  result.median_coverage = eigen::Median(coverage.count);
  const auto n_below = static_cast<u64>((coverage.count.array() < threshold).count());
  result.pct_windows_below_threshold = 100.0 * static_cast<f64>(n_below) / static_cast<f64>(coverage.count.size());
  result.is_sufficient = result.median_coverage >= threshold;
  return result;
}

}  // namespace xoos::cnc
