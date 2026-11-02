#pragma once
#include <xoos/types/int.h>

namespace xoos::svc {
// Enum for including or excluding `duplex_lowbq` in calculated BAM features `duplex_dp` and `duplex_af`
enum class DuplexLowbqMode : u8 {
  kInclude,
  kExclude
};
}  // namespace xoos::svc
