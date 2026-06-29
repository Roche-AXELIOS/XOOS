#include "core/interval.h"

#include <xoos/types/int.h>

namespace xoos::alignment_metrics {

bool Interval::Overlaps(const Interval& other) const {
  return start < other.end && end > other.start;
}

bool Interval::Contains(const s64 position) const {
  return position >= start && position < end;
}

}  // namespace xoos::alignment_metrics
