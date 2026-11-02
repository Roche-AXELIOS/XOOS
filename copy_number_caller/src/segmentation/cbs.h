#pragma once

#include <Eigen/Dense>

#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "genomic-segments.h"
#include "segmentation/max-t.h"
#include "segmentation/segmentation-options.h"
#include "util/eigen-math.h"

namespace xoos::cnc {

bool CircularPermuteNullAndTestSignificance(const Eigen::VectorXd& obvs,
                                            size_t n_permutations,
                                            f64 max_p,
                                            const Segment& max_t_seg,
                                            size_t max_obs_size,
                                            const vec<u64>& boundaries,
                                            size_t min_obs_per_segment,
                                            eigen::ArmaRng& rng,
                                            CbsMaxTMethod max_t_method = CbsMaxTMethod::kBruteForce);
bool BinaryPermuteNullAndTestSignificance(const Eigen::VectorXd& obvs,
                                          size_t n_permutations,
                                          f64 max_p,
                                          const Segment& seg1,
                                          const Segment& seg2,
                                          f64 min_t_for_automatic_reject_null,
                                          eigen::ArmaRng& rng);
f64 SiegmundPNuFull(f64 x);
f64 SiegmundP(f64 t, f64 m, f64 k);
bool HybridTestForSignificance(const Eigen::VectorXd& obvs,
                               size_t n_permutations,
                               f64 max_p,
                               Segment max_t_seg,
                               const vec<u64>& boundaries,
                               size_t min_obs_per_segment,
                               eigen::ArmaRng& rng,
                               CbsMaxTMethod max_t_method = CbsMaxTMethod::kBruteForce,
                               f64 min_t_for_automatic_reject_null = kSegmentationDefaultMinTForAutomaticSegmentation);
vec<Segment> GetSegments(const Eigen::VectorXd& obvs,
                         f64 max_p,
                         size_t n_permutations,
                         const vec<u64>& boundaries,
                         size_t min_obs_per_segment,
                         eigen::ArmaRng& rng,
                         CbsMaxTMethod max_t_method = CbsMaxTMethod::kBruteForce,
                         f64 min_t_for_automatic_reject_null = kSegmentationDefaultMinTForAutomaticSegmentation);
vec<Segment> CircularBinarySegmentation(const Eigen::VectorXd& obvs, const SegmentationOptions& options);
vec<Segment> MergeSingleObvsEndSegment(const vec<Segment>& sub_segments, size_t min_obs_per_segment);
Eigen::VectorXd TruncateOutliers(const Eigen::VectorXd& obvs);

}  // namespace xoos::cnc
