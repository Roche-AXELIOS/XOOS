#include "segmentation/segmentation-mode.h"

namespace xoos::cnc {
const StrMap<SegmentationMode> kStringToSegmentationMode{
    {"germline", SegmentationMode::kGermline},
    {"somatic", SegmentationMode::kSomatic},
};
}  // namespace xoos::cnc
