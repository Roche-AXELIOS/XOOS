#pragma once

#include <IITree.h>

#include <string>

#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>

#include "observations/observations.h"
#include "segmentation/genomic-segments.h"

namespace xoos::cnc {
class IntervalTrees {
 public:
  IntervalTrees() = default;
  IntervalTrees(const IntervalTrees&) = default;
  IntervalTrees(IntervalTrees&&) = default;
  IntervalTrees& operator=(const IntervalTrees&) = default;
  IntervalTrees& operator=(IntervalTrees&&) = default;
  ~IntervalTrees() = default;
  explicit IntervalTrees(const vec<std::string>& regions);
  explicit IntervalTrees(const Observations& obvs);
  explicit IntervalTrees(const vec<GenomicSegment>& segments);
  /// Look up intervals overlapping [start, end) using 0-based half-open coordinates.
  vec<size_t> LookUp(std::string contig, size_t start, size_t end) const;

 private:
  // unordered map to store IITree (interval trees) for each contig
  // std::string: contig name
  // size_t index type for IITree
  // size_t data type for IITree
  using IntervalTreesT = StrUnorderedMap<IITree<size_t, size_t>>;
  IntervalTreesT _trees;
};
}  // namespace xoos::cnc
