#pragma once

#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "io/alignment.h"

namespace xoos::read_collapser {

void DownsampleReadsInCluster(vec<AlignmentPtr>& reads, u32 max_reads);

}  // namespace xoos::read_collapser
