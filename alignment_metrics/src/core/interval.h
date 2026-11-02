#pragma once

#include <compare>
#include <tuple>

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
  // Bit vector of region labels this interval belongs to (one bit per label; see RegionLabelsSet).
  // A value of 0 means the interval is unlabeled (non-stratified mode). Ordering and equality below
  // deliberately ignore this field so that adding labels does not perturb the sorting and
  // binary-search invariants that downstream code (e.g. RegionLookupTable) relies on.
  u16 labels{0};

  Interval(const s64 in_start, const s64 in_end) : start(in_start), end(in_end) {
  }

  Interval(const s64 in_start, const s64 in_end, const u16 in_labels)
      : start(in_start), end(in_end), labels(in_labels) {
  }

  // Order and compare intervals by (start, end) only; `labels` is intentionally excluded so that
  // sorting and the sorted-order binary searches over intervals are unaffected by label membership.
  std::strong_ordering operator<=>(const Interval& other) const {
    return std::tie(start, end) <=> std::tie(other.start, other.end);
  }

  bool operator==(const Interval& other) const {
    return start == other.start && end == other.end;
  }

  // Check if this interval overlaps with another interval
  bool Overlaps(const Interval& other) const;

  // Check if this interval contains a specific position
  bool Contains(s64 position) const;
};
}  // namespace xoos::alignment_metrics
