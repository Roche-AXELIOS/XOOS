#include "segmentation/interval-trees.h"

#include <algorithm>
#include <cstddef>
#include <string>

#include <xoos/log/logging.h>
#include <xoos/types/vec.h>

#include "observations/observations.h"
#include "segmentation/genomic-segments.h"
#include "util/parse-util.h"

namespace xoos::cnc {

/**
 * @brief creates interval trees for list of region strigns
 * @param regions region strings of the format contig:start-end (1-based closed)
 * @return map, key=contig, item=interval tree for segments in the contig. data value for each node in the tree
 * corresponds to index w/i `segments`
 */
IntervalTrees::IntervalTrees(const vec<std::string>& regions) {
  for (const auto& region : regions) {
    const auto& [contig, start, end] = ParseRegionString(region);
    // ParseRegionString already converts to 0-based half-open
    _trees[contig].add(start, end, regions.size());
  }
  for (auto& [contig_name, interval_tree] : _trees) {
    interval_tree.index();
  }
}

/**
 * @brief creates interval trees for Observations
 * @param obvs input observations
 * @return map, key=contig, item=interval tree for segments in the contig. data value for each node in the tree
 * corresponds to index w/i `segments`
 */
IntervalTrees::IntervalTrees(const Observations& obvs) {
  for (size_t i = 0; i < obvs.contigs.size(); ++i) {
    const auto eidx = static_cast<Eigen::Index>(i);
    _trees[obvs.contigs[i]].add(obvs.starts[eidx], obvs.ends[eidx], i);
  }
  for (auto& [contig_name, interval_tree] : _trees) {
    interval_tree.index();
  }
}

/**
 * @brief creates interval trees for genomic segments
 * @param segments input segments
 * @return map, key=contig, item=interval tree for segments in the contigof corresponding nod. data value for each node
 * in the tree corresponds to index w/i `segments`
 */
IntervalTrees::IntervalTrees(const vec<GenomicSegment>& segments) {
  for (size_t i = 0; i < segments.size(); ++i) {
    _trees[segments[i].contig].add(segments[i].start, segments[i].end, i);
  }
  for (auto& pair : _trees) {
    pair.second.index();
  }
}

/**
 * @brief  looks up a <contig start end> pair inside _trees.
 * @param contig
 * @param start
 * @param end
 * @return idxs within original from which tree was constructed that represent overlaps to input
 */
vec<size_t> IntervalTrees::LookUp(std::string contig, size_t start, size_t end) const {
  const auto tree_ptr = _trees.find(contig);
  if (tree_ptr == _trees.end()) {
    Logging::Debug("contig {} is not observed in IntervalTrees", contig);
    return {};
  }
  const auto& tree = tree_ptr->second;
  // overlaps will contain the indexes in seed_segments that obvs[i] overlaps
  vec<size_t> overlaps;
  // IITree overlap() method fills in `overlaps` std::vector instead of returning a value
  tree.overlap(start, end, overlaps);
  if (overlaps.empty()) {
    Logging::Debug("observation {}:{}-{} does not appear in IntervalTrees. skipping...", contig, start + 1, end);
    return {};
  }
  vec<size_t> idxs(overlaps.size());
  for (size_t i = 0; i < overlaps.size(); ++i) {
    idxs[i] = tree.data(overlaps[i]);
  }
  std::ranges::sort(idxs);
  return idxs;
}
}  // namespace xoos::cnc
