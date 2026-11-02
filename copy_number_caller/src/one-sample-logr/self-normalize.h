#pragma once

#include "coverage/coverage.h"
#include "observations/observations.h"

namespace xoos::cnc {
Observations SelfNormalizeCounts(const CoverageRecords& coverage);
}
