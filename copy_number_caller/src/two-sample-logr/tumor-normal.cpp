#include "two-sample-logr/tumor-normal.h"

#include <cstddef>

#include "coverage/coverage.h"
#include "observations/observations.h"
#include "two-sample-logr/calculate-logrs.h"

namespace xoos::cnc {
Observations ProcessTumorNormal(const CoverageRecords& tumor_cov,
                                const CoverageRecords& normal_cov,
                                size_t normal_min_coverage,
                                size_t tumor_min_coverage) {
  CoverageRecords tumor_cov_cpy = tumor_cov;    // copy
  CoverageRecords normal_cov_cpy = normal_cov;  // copy
  FilterAndNormalizeTumorNormal(tumor_cov_cpy, normal_cov_cpy, normal_min_coverage, tumor_min_coverage);
  Observations ret = CalculateLogrs(tumor_cov_cpy, normal_cov_cpy);
  return ret;
}
}  // namespace xoos::cnc
