#pragma once

#include <xoos/types/int.h>

namespace xoos::svc::phaser {

// A variant allele on a read or haplotype: unknown/ambiguous, reference, or alternate.
enum class Allele : u8 {
  kUnknown,
  kRef,
  kAlt
};

}  // namespace xoos::svc::phaser
