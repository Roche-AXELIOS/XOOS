#pragma once

#include <xoos/types/int.h>

namespace xoos::cnc {

enum class VcfPuritySource : u8 {
  kScna,
  kSomaticSnv,
  kUserInput,
  kNone
};

}  // namespace xoos::cnc
