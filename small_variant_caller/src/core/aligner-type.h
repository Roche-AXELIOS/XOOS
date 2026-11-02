#pragma once
#include <xoos/types/int.h>

namespace xoos::svc {

// Aligner type for tumor-normal-wgs model selection.
// kCustom requires the user to explicitly provide --model and score thresholds.
enum class AlignerType : u8 {
  kBwa,
  kGiraffe,
  kCustom
};

}  // namespace xoos::svc
