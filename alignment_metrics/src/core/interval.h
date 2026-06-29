#pragma once

#include <compare>

#include <xoos/types/int.h>

namespace xoos::alignment_metrics {

/**
 * A genomic interval on a chromosome defined by a start and end position.
 * The interval is 0-based and half-open, meaning the start position is inclusive and the end position is exclusive.
 * For example, the interval [0, 5) includes positions 0, 1, 2, 3, and 4, but not position 5.
 */
struct Interval {
  // 0-based, inclusive start position of the genomic interval
  s64 start;
  // 0-based, exclusive end position of the genomic interval
  s64 end;

  Interval(const s64 in_start, const s64 in_end) : start(in_start), end(in_end) {
  }

  std::strong_ordering operator<=>(const Interval&) const = default;

  // Check if this interval overlaps with another interval
  bool Overlaps(const Interval& other) const;

  // Check if this interval contains a specific position
  bool Contains(s64 position) const;
};
}  // namespace xoos::alignment_metrics
