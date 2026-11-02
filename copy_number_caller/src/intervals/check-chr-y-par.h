#pragma once

#include <xoos/types/fs.h>
#include <xoos/types/vec.h>

#include "intervals/augmented-intervals.h"
#include "segmentation/genomic-segments.h"

namespace xoos::cnc {
/**
 * @brief Checks whether the chrY PAR region was masked during alignment by looking for any alignments in the chrY PAR
 * region defined in the seed segments file. If any alignments are found in the chrY PAR region, a warning message is
 * logged since this could significantly affect copy number predictions in allosomes. If the seed segments do not
 * specify any segments in the chrY PAR region, a warning is logged since we cannot check whether the chrY PAR region
 * was masked during alignment.
 * @param bam_file BAM file to check for alignments in the chrY PAR region
 * @param seed_segments Seed segments to check for segments in the chrY PAR region
 * @return true if alignments are found in the chrY PAR region, false otherwise
 */
bool CheckChrYPARUnmasked(const fs::path& bam_file, const vec<GenomicSegment>& seed_segments);

/**
 * @brief Removes intervals that overlap with the chrY PAR regions defined in the seed segments file if
 * CheckChrYPARUnmasked returns true. If CheckChrYPARUnmasked returns false, the original intervals are returned without
 * any modifications.
 * @param bam_file BAM file to check for alignments in the chrY PAR region
 * @param augmented_intervals AugmentedIntervals to filter if chrY PAR region is unmasked
 * @param seed_segments Seed segments to check for segments in the chrY PAR region
 * @return filtered AugmentedIntervals with intervals that overlap with the chrY PAR regions removed if
 * CheckChrYPARUnmasked returns true, original AugmentedIntervals otherwise
 */
AugmentedIntervals RemovePARIfChrYPARUnmasked(const fs::path& bam_file,
                                              const AugmentedIntervals& augmented_intervals,
                                              const vec<GenomicSegment>& seed_segments);
}  // namespace xoos::cnc
