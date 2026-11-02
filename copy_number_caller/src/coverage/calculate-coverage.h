#pragma once

#include <xoos/io/alignment-reader.h>
#include <xoos/types/int.h>

#include "cli/copy-number-caller-options.h"
#include "coverage/coverage.h"

namespace xoos::cnc {

// this will handle all the IO logic and load data structures to be passed into CalculateCoverageRegion.
void CalculateCoverageMain(const CopyNumberCallerOptions& options);
void CalculateCoverageRegion(const io::AlignmentReader& reader,
                             CoverageRecords& res,
                             u32 exclude_flags,
                             bool ignore_dn,
                             s32 start_idx,
                             s32 step);
}  // namespace xoos::cnc
