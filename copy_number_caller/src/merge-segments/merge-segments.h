#pragma once

#include <xoos/types/float.h>
#include <xoos/types/vec.h>

#include "segmentation/genomic-segments.h"

namespace xoos::cnc {
GenomicSegment MergeSegments(const vec<GenomicSegment>& segments_to_merge);
vec<GenomicSegment> MergeLowMapqSegments(const vec<GenomicSegment>& segments, f64 min_mapq_cutoff);
vec<GenomicSegment> MergeAdjacentEqualCopyNumberSegments(const vec<GenomicSegment>& segments);
}  // namespace xoos::cnc
