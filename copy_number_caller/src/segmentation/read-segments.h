#pragma once

#include <xoos/types/vec.h>

#include "segmentation/genomic-segments.h"
#include "segmentation/segment-type.h"
#include "segmentation/segments-header.h"

namespace xoos::cnc {

vec<GenomicSegment> ReadSegments(const fs::path& fname, SegmentType segment_type);
SegmentsHeader ReadHeaderFromSegments(const fs::path& fname);
}  // namespace xoos::cnc
