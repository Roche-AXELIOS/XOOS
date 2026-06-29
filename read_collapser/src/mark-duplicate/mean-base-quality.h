#pragma once

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/float.h>
#include <xoos/types/vec.h>

#include "io/alignment.h"

namespace xoos::read_collapser {

/**
 * Calculate the mean base quality of the bases in the alignment.
 */
f64 MeanBaseQ(const AlignmentPtr& alignment);

/**
 * Find the alignment in @p alignments with the maximum mean base quality. If multiple alignments have the same mean
 * base quality, the first alignment with the maximum mean base quality is returned. If @p alignments is empty, nullptr.
 *
 * This method is used to determine which alignment should not be marked as a duplicate during duplicate marking.
 */
AlignmentPtr FindAlignmentWithMaxMeanBaseQ(const vec<AlignmentPtr>& alignments);

}  // namespace xoos::read_collapser
