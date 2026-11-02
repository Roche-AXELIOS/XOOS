#pragma once

#include <xoos/types/int.h>

namespace xoos::cnc {
enum class CopyNumberCallerModes : u8 {
  kSomaticTumorNormalWGS,
  kSomaticTumorOnlyWGS,
  kGermlineNormalWGS,
  kUnknown
};

// Modes that share the somatic VCF layout (allele-specific FORMAT fields, somatic copy-number status). Tumor-only mode
// reuses the somatic output apart from the additional GERMLINE INFO flag.
inline bool IsSomaticLikeMode(const CopyNumberCallerModes mode) {
  return (mode == CopyNumberCallerModes::kSomaticTumorNormalWGS) ||
         (mode == CopyNumberCallerModes::kSomaticTumorOnlyWGS);
}
}  // namespace xoos::cnc
