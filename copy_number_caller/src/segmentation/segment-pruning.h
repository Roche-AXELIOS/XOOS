#pragma once

#include <optional>

#include <xoos/types/float.h>
#include <xoos/types/vec.h>

#include "observations/observations.h"
#include "segmentation/genomic-segments.h"

namespace xoos::cnc {

vec<GenomicSegment> BreakLogrSegmentsByDHSegments(const vec<GenomicSegment>& logr_segments,
                                                  const Observations& logrs,
                                                  const vec<GenomicSegment>& dh_segments,
                                                  const Observations& dhs);
vec<GenomicSegment> PruneSegmentsByClustering(const vec<GenomicSegment>& logr_segments,
                                              bool joint_clustering,
                                              std::optional<f64> h,
                                              f64 baf_weight);
vec<GenomicSegment> PruneSegmentsByClusteringPerSeedSegment(const vec<GenomicSegment>& seed_segments,
                                                            const vec<GenomicSegment>& segments_to_prune,
                                                            bool joint_clustering,
                                                            std::optional<f64> h,
                                                            f64 baf_weight);

}  // namespace xoos::cnc
