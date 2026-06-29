#pragma once

namespace xoos::demux::strand {

/**
 * @brief Defines the possible classification outcomes for a strand.
 */
enum class StrandType {
  kUnknown,     // Equal fw/rv counts or no informative observations
  kForward,     // More fw but not significant
  kReverse,     // More rv but not significant
  kForwardSig,  // Statistically significant evidence for forward
  kReverseSig   // Statistically significant evidence for reverse
};

};  // namespace xoos::demux::strand
