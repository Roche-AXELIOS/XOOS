#include "xoos/io/htslib-util/cigar-util.h"

#include <htslib/sam.h>

#include <xoos/types/int.h>

namespace xoos::io {

bool ConsumesRead(const u32 op) {
  static constexpr u32 kQueryBitMask = 0b01;
  // CIGAR ops: M, I, S, =, X
  return (bam_cigar_type(op) & kQueryBitMask) != 0;
}

bool ConsumesReference(const u32 op) {
  static constexpr u32 kReferenceBitMask = 0b10;
  // CIGAR ops: M, D, N, =, X
  return (bam_cigar_type(op) & kReferenceBitMask) != 0;
}

}  // namespace xoos::io
