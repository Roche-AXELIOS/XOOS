#pragma once

#include <xoos/sex_predict/sex.h>
#include <xoos/types/float.h>
#include <xoos/types/vec.h>

#include "likelihood/likelihood-options.h"
#include "segmentation/genomic-segments.h"

namespace xoos::cnc {
void FlagCallsByMeanMapq(vec<GenomicSegment>& segments, f64 mapq_cutoff_for_calls);
void FlagCallsByExpectedTotalCopyNumber(vec<GenomicSegment>& segments,
                                        const sex_predict::Sex& sex,
                                        const NonIntegerTcnFilterOptions& options);
void FlagCallsByLength(vec<GenomicSegment>& segments, size_t cnv_length_flag_min_size);
void FlagChrYCallsIfFemale(vec<GenomicSegment>& segments);
}  // namespace xoos::cnc
