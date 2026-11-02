#pragma once

#include <xoos/types/int.h>
#include <xoos/types/str-container.h>

namespace xoos::cnc {
enum class SegmentationMode : u8 {
  kGermline,
  kSomatic,
  kUnknown
};
extern const StrMap<SegmentationMode> kStringToSegmentationMode;
}  // namespace xoos::cnc
