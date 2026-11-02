#pragma once
#include <xoos/types/int.h>

#include "cli/copy-number-caller-options.h"
#include "intervals/augmented-intervals.h"

namespace xoos::cnc {
void AugmentIntervalsMain(const CopyNumberCallerOptions& options);
AugmentedIntervals GenerateAndAugmentOnAndOffTargetIntervals(const CopyNumberCallerOptions& options);
AugmentedIntervals GenerateAndAugmentWholeGenomeIntervals(const CopyNumberCallerOptions& options);
}  // namespace xoos::cnc
