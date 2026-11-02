#pragma once
#include <xoos/types/int.h>

namespace xoos::svc {

// Sample preparation type for tumor-normal model selection.
enum class SampleType : u8 {
  kAutoDetect,
  kFfpe,
  kCellLine
};

}  // namespace xoos::svc
