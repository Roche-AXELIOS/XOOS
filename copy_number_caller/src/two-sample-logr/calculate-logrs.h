#pragma once
#include "coverage/coverage.h"
#include "observations/observations.h"

namespace xoos::cnc {
const size_t kNormalMinCov = 30;
// We do not put a minimum coverage threshold on the tumour because we would not be able to pick up homozygous deletions
// in tumour cell-lines. An example of this would be a PTEN deletion (chr10:87940542-87952584) in COLO829T.
const size_t kTumorMinCov = 0;

void FilterAndNormalizeTumorNormal(CoverageRecords& tumor,
                                   CoverageRecords& normal,
                                   size_t normal_min_cov = kNormalMinCov,
                                   size_t tumor_min_cov = kTumorMinCov);
Observations CalculateLogrs(const CoverageRecords& tumor, const CoverageRecords& normal);
}  // namespace xoos::cnc
