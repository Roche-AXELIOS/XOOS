#pragma once

#include <xoos/types/int.h>

namespace xoos::cnc {

enum class TotalCopyNumberPrior : u8 {
  kDiploidBias,
  kUniform
};

}  // namespace xoos::cnc
