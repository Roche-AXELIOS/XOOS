#pragma once
#include <optional>

#include <xoos/types/float.h>

#include "coverage/coverage.h"

namespace xoos::cnc {

constexpr f64 kMinMedianCoverageThreshold = 10.0;

struct CoverageCheckResult {
  f64 median_coverage;
  f64 pct_windows_below_threshold;
  bool is_sufficient;
};

/// Returns true when a coverage check was performed and the result was insufficient.
inline bool IsLowCoverage(const std::optional<CoverageCheckResult>& result) {
  return result.has_value() && !result->is_sufficient;
}

/**
 * @brief Check whether a sample has sufficient coverage for copy number calling.
 *
 * Computes the median of the count vector across all windows and the percentage
 * of windows with coverage below the threshold.
 *
 * @param coverage Coverage records from BAM
 * @param threshold Minimum median coverage required (default 10x)
 * @return CoverageCheckResult with median, percentage below threshold, and pass/fail
 */
CoverageCheckResult CheckMedianCoverage(const CoverageRecords& coverage, f64 threshold = kMinMedianCoverageThreshold);

}  // namespace xoos::cnc
