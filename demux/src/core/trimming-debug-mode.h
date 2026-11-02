#pragma once

#include <algorithm>
#include <vector>

namespace xoos::demux {

// Modes for the hidden --trimming-debug option. Multiple modes can be combined.
enum class TrimmingDebugMode {
  // Append sq:Z tag with raw trimmed sequences.
  kRawSeq,
};

// Check whether a specific debug mode is active.
inline bool HasTrimmingDebugMode(const std::vector<TrimmingDebugMode>& modes, const TrimmingDebugMode mode) {
  return std::ranges::find(modes, mode) != modes.end();
}

}  // namespace xoos::demux
