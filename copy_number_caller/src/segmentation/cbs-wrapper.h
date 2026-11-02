#pragma once
#include <Eigen/Dense>

#include <xoos/types/vec.h>

#include "segmentation/genomic-segments.h"
#include "segmentation/segmentation-options.h"

namespace xoos::cnc {
vec<Segment> CircularBinarySegmentationWrapper(const Eigen::VectorXd& obvs, const SegmentationOptions& options);
}
