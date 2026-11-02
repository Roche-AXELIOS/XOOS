#include "segmentation/cbs-wrapper.h"

#include <Eigen/Dense>
#include <stdexcept>

#include <xoos/types/vec.h>

#include "segmentation/cbs-for-somatic.h"
#include "segmentation/cbs.h"
#include "segmentation/genomic-segments.h"
#include "segmentation/segmentation-mode.h"
#include "segmentation/segmentation-options.h"

namespace xoos::cnc {
vec<Segment> CircularBinarySegmentationWrapper(const Eigen::VectorXd& obvs, const SegmentationOptions& options) {
  if (options.segmentation_mode == SegmentationMode::kGermline) {
    return CircularBinarySegmentation(obvs, options);
  }
  if (options.segmentation_mode == SegmentationMode::kSomatic) {
    return somatic::CircularBinarySegmentation(
        obvs, options.max_p, options.n_permutations, options.min_obs_per_segment, options.no_single_obvs_subsegments);
  }
  throw std::runtime_error("invalid segmentation mode");
}
}  // namespace xoos::cnc
