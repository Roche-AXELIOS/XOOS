#pragma once
#include <xoos/types/fs.h>

#include "cli/copy-number-caller-options.h"
#include "intervals/augmented-intervals.h"

namespace xoos::cnc {

AugmentedIntervals LoadOrGenerateWholeGenomeIntervals(const CopyNumberCallerOptions& options,
                                                      const fs::path& augmented_intervals_out);

}  // namespace xoos::cnc
