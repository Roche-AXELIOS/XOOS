#pragma once

#include <Eigen/Dense>

#include <xoos/enum/enum-util.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "segmentation/genomic-segments.h"
#include "util/eigen-math.h"

namespace xoos::cnc::somatic {

f64 TStat(f64 n_v1, f64 mean_v1, f64 var_v1, f64 n_v2, f64 mean_v2, f64 var_v2);
f64 TStatFromVecs(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2);
Segment GetMaxTSegCumSum(const Eigen::VectorXd& obvs, size_t max_n, size_t min_n);
bool CircularPermuteNullAndTestSignificance(const Eigen::VectorXd& obvs,
                                            size_t n_permutations,
                                            f64 max_p,
                                            f64 t,
                                            size_t max_obs_size,
                                            const vec<u64>& boundaries,
                                            size_t min_obs_size,
                                            eigen::ArmaRng& rng);
bool BinaryPermuteNullAndTestSignificance(const Eigen::VectorXd& obvs,
                                          size_t n_permutations,
                                          f64 max_p,
                                          const Segment& seg1,
                                          const Segment& seg2,
                                          eigen::ArmaRng& rng);
f64 SiegmundPNuFull(f64 x);
f64 SiegmundP(f64 t, f64 m, f64 k);
bool HybridTestForSignificance(const Eigen::VectorXd& obvs,
                               size_t n_permutations,
                               f64 max_p,
                               f64 t,
                               const vec<u64>& boundaries,
                               size_t min_obs_size,
                               eigen::ArmaRng& rng);
vec<Segment> GetSegments(const Eigen::VectorXd& obvs,
                         f64 max_p,
                         size_t n_permutations,
                         const vec<u64>& boundaries,
                         size_t min_obs_size,
                         eigen::ArmaRng& rng);
vec<Segment> CircularBinarySegmentation(const Eigen::VectorXd& obvs,
                                        f64 max_p,
                                        size_t n_permutations,
                                        size_t min_obs_size,
                                        bool no_single_obvs_subsegments);
vec<Segment> MergeSingleObvsEndSegment(const vec<Segment>& sub_segments, size_t min_obs_size);
}  // namespace xoos::cnc::somatic
