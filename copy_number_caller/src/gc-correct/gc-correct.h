#pragma once
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "cli/copy-number-caller-options.h"

namespace xoos::cnc {

void GCCorrectMain(const CopyNumberCallerOptions& options);

}  // namespace xoos::cnc
