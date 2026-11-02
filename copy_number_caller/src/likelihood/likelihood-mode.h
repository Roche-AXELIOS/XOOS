#pragma once

#include <xoos/types/int.h>

namespace xoos::cnc {
enum class LikelihoodMode : u8 {
  kUnknown,
  kSomatic,
  kGermline
};
}  // namespace xoos::cnc
