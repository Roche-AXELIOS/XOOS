#pragma once

#include "coverage/coverage.h"
#include "observations/observations.h"

namespace xoos::cnc {

Observations ProcessTumorNormal(const CoverageRecords& tumor_cov,
                                const CoverageRecords& normal_cov,
                                size_t normal_min_coverage,
                                size_t tumor_min_coverage);
}

// namespace xoos::cnc
