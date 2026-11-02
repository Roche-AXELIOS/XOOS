#pragma once

#include <xoos/types/int.h>

namespace xoos::cnc {
enum class SegmentType : u8 {
  kUnknown,
  kSeed,
  kLogROnly,
  kBaf,
  kGermlineLikelihood,
  kSomaticWithBafLikelihood,
  kSomaticNoBafLikelihood,
  kSomaticTumorOnlyWithBafLikelihood
};

}  // namespace xoos::cnc
