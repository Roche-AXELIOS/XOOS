#pragma once
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>

namespace xoos::cnc {

// A genomic interval using 0-based half-open coordinates [start, end).
struct Interval {
  // 0-based inclusive start position
  size_t start;
  // 0-based exclusive end position
  size_t end;
  auto operator<=>(const Interval&) const = default;
};

using ContigToIntervals = StrUnorderedMap<vec<Interval>>;

};  // namespace xoos::cnc
