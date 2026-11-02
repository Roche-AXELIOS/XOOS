#pragma once

#include <xoos/types/int.h>

namespace xoos::cnc {
enum class LikelihoodModel : u8 {
  kSerialSummarized,
  kLogrSummarizedOnly
};
}  // namespace xoos::cnc
