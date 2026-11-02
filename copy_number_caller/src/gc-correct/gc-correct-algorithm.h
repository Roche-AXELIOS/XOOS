#pragma once

#include <Eigen/Dense>

#include <xoos/types/float.h>

#include "coverage/coverage.h"
#include "intervals/augmented-intervals.h"

namespace xoos::cnc {

/**
 * @struct GCCorrectResult
 * @brief GC-corrected coverage, one entry per input region, in the same order as the input.
 */
struct GCCorrectResult {
  Eigen::VectorXd count;
  Eigen::VectorXd total_coverage;
  Eigen::VectorXd average_coverage;
};

/**
 * @brief Correct coverage for GC bias using xoos::gc_correct.
 *
 * The LOESS model is trained only on autosomal regions with mappability >= 1 (matching the legacy
 * behavior of this caller), then applied to every region in `coverage`. Regions that
 * xoos::gc_correct considers invalid for correction (non-positive average coverage or negative
 * GC bias) keep their original, uncorrected values instead of the zero-fill xoos::gc_correct
 * produces for them, matching the legacy pass-through behavior for those regions.
 *
 * @param coverage             - raw coverage, one entry per region
 * @param augmented_intervals  - per-region GC bias and mappability annotations, aligned with `coverage`
 * @param first_span           - span parameter (0, 1] for the first round of LOESS smoothing. See
 * xoos::gc_correct::TrainGCCorrectModel for details.
 * @return GC-corrected coverage, one entry per input region
 */
GCCorrectResult GCCorrect(const CoverageRecords& coverage,
                          const AugmentedIntervals& augmented_intervals,
                          f64 first_span);

}  // namespace xoos::cnc
