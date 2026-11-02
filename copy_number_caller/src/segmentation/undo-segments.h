#pragma once
#include <Eigen/Dense>

#include <xoos/types/float.h>
#include <xoos/types/vec.h>

#include "observations/observations.h"
#include "segmentation/genomic-segments.h"

namespace xoos::cnc {
f64 GetTrimmedStandardDeviation(const Eigen::VectorXd& arr);
/**
 * @brief Undo genomic segments by merging adjacent segments if their median observation difference is less than a
 * factor (max_trimmed_sd_factor) of the trimmed standard deviation (trimmed_sd)
 * @param segments input segments
 * @param obvs observations vector
 * @param trimmed_sd trimmed standard deviation of the observations (use GetTrimmedStandardDeviation to calculate)
 * @param max_trimmed_sd_factor maximum factor of trimmed standard deviation to consider for merging segments
 * @return merged segments
 */
/// Merge two Observations by concatenating all fields in order: [obvs1, obvs2].
Observations MergeObservations(const Observations& obvs1, const Observations& obvs2);

vec<GenomicSegment> UndoSegments(const vec<GenomicSegment>& segments,
                                 const Observations& obvs,
                                 f64 trimmed_sd,
                                 f64 max_trimmed_sd_factor);
vec<GenomicSegment> IterativeUndoGenomicSegments(const vec<GenomicSegment>& seed_segments,
                                                 const vec<GenomicSegment>& segments,
                                                 const Observations& obvs,
                                                 f64 start_max_trimmed_sd_factor,
                                                 f64 increment_max_trimmed_sd_factor,
                                                 size_t max_num_segments);

/**
 * @brief Undo segments in a single, non-iterative merge pass at a fixed max_trimmed_sd_factor. Unlike
 * IterativeUndoGenomicSegments, the factor is not escalated and no target segment count is applied. Note that except
 * for the contig, start, and end fields of each GenomicSegment member, all other fields are invalidated.
 * @param seed_segments seed segments
 * @param segments input segments
 * @param obvs observations vector
 * @param max_trimmed_sd_factor max trimmed sd factor used for the single merge pass
 * @return merged segments
 */
vec<GenomicSegment> NonIterativeUndoGenomicSegments(const vec<GenomicSegment>& seed_segments,
                                                    const vec<GenomicSegment>& segments,
                                                    const Observations& obvs,
                                                    f64 max_trimmed_sd_factor);
}  // namespace xoos::cnc
