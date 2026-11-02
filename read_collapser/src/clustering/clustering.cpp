#include "clustering/clustering.h"

#include <iterator>
#include <optional>
#include <ranges>
#include <stack>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <xoos/error/error.h>
#include <xoos/read-name/umi.h>
#include <xoos/util/hash.h>
#include <xoos/util/math.h>
#include <xoos/util/sequence-functions.h>

#include "core/read-collapser-options.h"
#include "io/alignment-io.h"
#include "io/alignment.h"
#include "metrics/metrics.h"
#include "util/duplex-util.h"

namespace xoos::read_collapser {

size_t ClusterIdHash::operator()(const ClusterId& cluster_id) const {
  return util::hash::Hash(cluster_id.super_region_id, cluster_id.count);
}

ClusterStrand DetermineClusterStrand(const AlignmentPtr& alignment, const bool cluster_by_strand) {
  using enum ClusterStrand;
  if (!cluster_by_strand) {
    return kBoth;
  }
  return (alignment->record->core.flag & BAM_FREVERSE) != 0 ? kRev : kFwd;
}

size_t UmiAndStrandHash::operator()(const UmiAndStrand& umi_strand) const {
  return xoos::util::hash::Hash(umi_strand.umi, static_cast<u8>(umi_strand.strand));
}

ClusterCoord CreateClusterCoord(const AlignmentPtr& alignment, const bool cluster_by_strand) {
  return {
      .start = alignment->StartPos(),
      .end = alignment->EndPos(),
      .strand = DetermineClusterStrand(alignment, cluster_by_strand),
  };
}

size_t ClusterCoordHash::operator()(const ClusterCoord& coord) const {
  return util::hash::Hash(coord.start, coord.end, coord.strand);
}

// Update duplicate read metrics based on the clusters provided.
namespace {

void UpdateDuplicateReadMetrics(const Clusters& clusters, Metrics& metrics) {
  for (const auto& cluster : clusters | std::views::values) {
    const auto cluster_size = cluster->alignments.size();
    // All but one read in a cluster are duplicates (pre-duplex deconvolution).
    metrics.duplicate_reads += cluster_size - 1;
  }
}

}  // namespace

// Update the cluster metrics based on the clusters provided.
void UpdateClusterMetrics(const Clusters& clusters, Metrics& metrics) {
  // update the histogram stats for each cluster
  for (const auto& cluster : clusters | std::views::values) {
    const auto cluster_size = cluster->alignments.size();
    ++metrics.total_clusters;
    metrics.UpdateClusterSizeHistogram(cluster_size, "total_clusters", 1);
    if (cluster_size == 1) {
      ++metrics.singleton_clusters;
    }
    bool cluster_is_full = false;
    bool cluster_is_partial = false;
    bool cluster_is_forward = false;
    bool cluster_is_reverse = false;
    for (const auto& read : cluster->alignments) {
      metrics.UpdateClusterSizeHistogram(cluster_size, "total_reads", 1);
      ++metrics.clustering_reads;
      if (read->IsPartial()) {
        cluster_is_partial = true;
        ++metrics.clustering_partial_reads;
      } else {
        cluster_is_full = true;
        ++metrics.clustering_full_reads;
      }
      if (read->IsForward()) {
        cluster_is_forward = true;
        metrics.UpdateClusterSizeHistogram(cluster_size, "forward_strand_reads", 1);
      } else {
        cluster_is_reverse = true;
      }
    }
    if (cluster_is_full && cluster_is_partial) {
      ++metrics.full_and_partial_read_clusters;
      metrics.UpdateClusterSizeHistogram(cluster_size, "full_and_partial_read_clusters", 1);
    } else if (cluster_is_full) {
      ++metrics.full_read_clusters;
      metrics.UpdateClusterSizeHistogram(cluster_size, "full_read_clusters", 1);
    } else if (cluster_is_partial) {
      ++metrics.partial_read_clusters;
      metrics.UpdateClusterSizeHistogram(cluster_size, "partial_read_clusters", 1);
    }
    if (cluster_is_forward && cluster_is_reverse) {
      ++metrics.mixed_strand_clusters;
      metrics.UpdateClusterSizeHistogram(cluster_size, "mixed_strand_clusters", 1);
    } else if (cluster_is_forward) {
      ++metrics.forward_strand_clusters;
      metrics.UpdateClusterSizeHistogram(cluster_size, "forward_strand_clusters", 1);
    } else if (cluster_is_reverse) {
      ++metrics.reverse_strand_clusters;
      metrics.UpdateClusterSizeHistogram(cluster_size, "reverse_strand_clusters", 1);
    }
  }
}

void MergePreliminaryClusters(const ClusterCoord& initial_coord,
                              const PreliminaryCluster& preliminary_cluster,
                              PreliminaryClusters& preliminary_clusters,
                              const u32 wiggle_room,
                              const u32 max_depth) {
  // Recursive depth first search is implemented using a stack rather than actual function recursion to reduce
  // the risk of stack overflow.
  std::stack<ClusterCoord> coords;
  coords.emplace(initial_coord);

  while (!coords.empty() && (max_depth == 0 || coords.size() < max_depth)) {
    auto [start, end, strand] = coords.top();
    coords.pop();

    const u32 neighbor_start_lower = math::SatSub(start, wiggle_room);
    const u32 neighbor_start_upper = start + wiggle_room;

    // Perform a walk of all potential neighbors of the current preliminary cluster. To do this we generate all
    // cluster coordinates which meet the constraint |a.start - b.start| + |a.end - b.end| <= wiggle_room, and
    // search for preliminary clusters at these coordinates. If any are found, we merge them into the current cluster
    // and add their coordinates to the stack for further searching.
    for (u32 neighbor_start = neighbor_start_lower, neighbor_end_lower = end, neighbor_end_upper = end;
         neighbor_start <= neighbor_start_upper;
         neighbor_start++) {
      for (u32 neighbor_end = neighbor_end_lower; neighbor_end <= neighbor_end_upper; neighbor_end++) {
        ClusterCoord neighbor_coord = {.start = neighbor_start, .end = neighbor_end, .strand = strand};

        // If there is no preliminary cluster at this coordinate, or if the preliminary cluster has already been merged
        // into a cluster, skip this neighbor.
        const auto neighbor_it = preliminary_clusters.find(neighbor_coord);
        if (neighbor_it == preliminary_clusters.end() || neighbor_it->second.cluster != nullptr) {
          continue;
        }

        // We have found a neighboring preliminary cluster, merge it into the current cluster.
        neighbor_it->second.cluster = preliminary_cluster.cluster;
        for (auto& alignment : neighbor_it->second.alignments) {
          alignment->cluster = preliminary_cluster.cluster.get();
          alignment->cluster->alignments.emplace_back(alignment);
        }

        // Extend the cluster further by searching for neighbors of this current neighbor.
        coords.emplace(neighbor_coord);
      }

      // Ensure the constraint |a.start - b.start| + |a.end - b.end| <= wiggle_room is maintained.
      // As neighbor_start increases, towards the range of neighbor_end_lower and neighbor_end_upper must increase,
      // and as neighbor_start passes start the range of neighbor_end_lower and neighbor_end_upper must decrease.
      if (neighbor_start < start) {
        --neighbor_end_lower;
        ++neighbor_end_upper;
      } else {
        ++neighbor_end_lower;
        --neighbor_end_upper;
      }
    }
  }
}

/**
 * Create clusters from a range of alignments.
 *
 * Note: The range defined by [@p cluster_start_index, @p cluster_end_index) is inclusive-exclusive.
 */
namespace {

void CreatePositionalClusters(const u32 wiggle_room,
                              const bool cluster_by_strand,
                              const CreateClusterId& create_cluster_id,
                              Clusters& clusters,
                              const std::forward_iterator auto& cluster_start,
                              const std::forward_iterator auto& cluster_end) {
  // Iterate through all the alignments that we should consider for clustering, if the ClusterCoord for
  // the alignment is not in the cluster_map, create a new cluster and add it to the cluster_map. Otherwise,
  // add the alignment to the existing cluster.
  // This process creates preliminary clusters. Each preliminary cluster contains all reads that have the same
  // alignment start and end positions as well as strand.
  PreliminaryClusters cluster_map;
  for (auto alignment_it = cluster_start; alignment_it != cluster_end; ++alignment_it) {
    ClusterCoord coord = CreateClusterCoord(*alignment_it, cluster_by_strand);
    if (auto cluster = cluster_map.find(coord); cluster == cluster_map.end()) {
      cluster_map.try_emplace(cluster, coord, PreliminaryCluster{nullptr, {*alignment_it}});
    } else {
      cluster->second.alignments.emplace_back(*alignment_it);
    }
  }

  for (auto& [coords, cluster_tmp] : cluster_map) {
    if (cluster_tmp.cluster != nullptr) {
      continue;
    }

    auto cluster_id = create_cluster_id();
    cluster_tmp.cluster = std::make_shared<Cluster>(cluster_id, vec<AlignmentPtr>{});
    for (auto& read : cluster_tmp.alignments) {
      read->cluster = cluster_tmp.cluster.get();
      read->cluster->alignments.emplace_back(read);
    }

    MergePreliminaryClusters(coords, cluster_tmp, cluster_map, wiggle_room, 0);
    clusters[cluster_id] = cluster_tmp.cluster;
  }
}

/**
 * Positional clustering of full reads, broken up into two phases:
 *   - Phase 1: Find all reads whose alignment start position is within the wiggle room of the previous alignment.
 *   - Phase 2: For each set of reads found in phase 1, determine the clusters they belong to based on their
 *   start, end, and strand.
 * Note: @p reads must be sorted by alignment start position.
 */
void ClusterFullReadsByPosition(const u32 wiggle_room,
                                const bool cluster_by_strand,
                                const CreateClusterId& create_cluster_id,
                                const vec<AlignmentPtr>& reads,
                                Clusters& clusters) {
  // The following code assumes there is at least one read, return early if there are none.
  if (reads.empty()) {
    return;
  }

  // The start of the first cluster is the first read, this will be updated for each new cluster we start.
  auto cluster_start_it = reads.begin();

  for (auto alignment_it = std::next(cluster_start_it); alignment_it != reads.end(); ++alignment_it) {
    // Check if the alignment is part of the current cluster by checking if it's start position is within the wiggle
    // room of the previous alignments start position. If it is, move onto the next alignment.
    auto prev_alignment = std::prev(alignment_it);
    if ((*alignment_it)->record->core.pos <= (*prev_alignment)->record->core.pos + wiggle_room) {
      // This alignment is part of the current cluster, move onto the next alignment.
      continue;
    }

    // The current alignment is not part of the cluster, so cluster the alignments from the start index up to the
    // current.
    CreatePositionalClusters(
        wiggle_room, cluster_by_strand, create_cluster_id, clusters, cluster_start_it, alignment_it);

    // Update the start of the next cluster to be the current alignment.
    cluster_start_it = alignment_it;
  }

  if (cluster_start_it != reads.end()) {
    // Handle the last cluster if it was not previously handled and there is at least one alignment.
    CreatePositionalClusters(
        wiggle_room, cluster_by_strand, create_cluster_id, clusters, cluster_start_it, reads.end());
  }
}

/**
 * @brief Clusters reads based on their UMI (Unique Molecular Identifier) pairs.
 *
 * This function processes a collection of reads and groups them into clusters
 * based on their UMI pairs. Each UMI pair is represented as a combination of 5' and 3' UMIs.
 * If a cluster for a given UMI pair does not exist, a new cluster is created. If a cluster
 * already exists for the UMI pair, the read is added to the existing cluster.
 *
 * @param create_cluster_id A callable object or function that generates unique cluster IDs.
 * @param reads A vector of reads to be clustered. Each read must
 *              contain valid 5' and 3' UMI values.
 *              The reads are modified in-place to assign cluster pointers to each read
 *              as they are clustered.
 * @param clusters A map that stores the resulting clusters, keyed by their unique cluster IDs.
 *                 This map is updated in-place as clusters are created and populated with reads.
 *
 * @note This function assumes that the `reads` vector is non-empty. If it is empty,
 *       the function returns early without performing any operations.
 * @note The `clusters` map is updated in-place with the newly created or updated clusters.
 */
void ClusterFullReadsByUmi(const CreateClusterId& create_cluster_id,
                           const vec<AlignmentPtr>& reads,
                           Clusters& clusters) {
  // The following code assumes there is at least one read, return early if there are none.
  if (reads.empty()) {
    return;
  }
  // Create a map of UMIs to clusters to determine if a cluster already exists for a given UMI pair.
  std::unordered_map<read_name::UmiPair, Cluster*, read_name::UmiPairHash> umi_cluster_map;
  for (const auto& read : reads) {
    const auto umi_pair = read_name::UmiPair{.umi5p = read->umi5p, .umi3p = read->umi3p};
    auto it = umi_cluster_map.find(umi_pair);
    // UMI pair does not exist in the map, create a new cluster.
    if (it == umi_cluster_map.end()) {
      auto new_cluster_id = create_cluster_id();
      auto new_cluster = std::make_shared<Cluster>(new_cluster_id, vec<AlignmentPtr>{});
      read->cluster = new_cluster.get();
      read->cluster->alignments.emplace_back(read);
      clusters[new_cluster_id] = new_cluster;
      umi_cluster_map[umi_pair] = new_cluster.get();
    } else {
      // Cluster exists for UMI pair and must be updated with read.
      // Assumes cluster already contains at least one read (from cluster creation step).
      read->cluster = it->second;
      read->cluster->alignments.emplace_back(read);
    }
  }
}

/**
 * @brief Assigns a partial read to the nearest cluster based on positional proximity.
 *
 * @param wiggle_room The maximum allowable positional deviation for a read to be considered part of a cluster.
 * @param cluster_by_strand Whether reads should be separated by strand.
 * @param candidate_clusters Vector of cluster statistics to search.
 * @param partial_read The partial read to be assigned to a cluster.
 *                     This read is modified in-place to update its cluster assignment if
 *                     a suitable cluster is found.
 */
void AssignPartialReadToCluster(const u32 wiggle_room,
                                const bool cluster_by_strand,
                                const vec<ClusterPtrAndStats>& candidate_clusters,
                                AlignmentPtr& partial_read) {
  f64 best_distance = std::numeric_limits<f64>::max();
  ClusterPtr best_cluster = nullptr;
  f64 second_best_distance = std::numeric_limits<f64>::max();

  // If a read is complete (untruncated) on the 5' end, we can use the start position to assign it to a cluster
  // Otherwise, if a read is complete on the 3' end, we can use the end position to assign it to a cluster.
  // We don't expect reads that are partial on both ends since they should have been discarded upstream
  const bool is_5p_complete = partial_read->IsFivePrimeComplete();
  const u32 pos = is_5p_complete ? partial_read->StartPos() : partial_read->EndPos();

  for (const auto& cluster_info : candidate_clusters) {
    if (const auto strand = DetermineClusterStrand(partial_read, cluster_by_strand);
        cluster_by_strand && cluster_info.strand != strand) {
      continue;
    }

    f64 distance = std::abs(static_cast<f64>(pos) - cluster_info.mean_pos);
    if (distance > wiggle_room) {
      continue;
    }

    // Check if the read is too long for the cluster
    if (is_5p_complete && partial_read->EndPos() > cluster_info.max_end + wiggle_room) {
      continue;
    }
    if (!is_5p_complete && partial_read->StartPos() < math::SatSub(cluster_info.min_start, wiggle_room)) {
      continue;
    }

    if (distance < best_distance) {
      second_best_distance = best_distance;
      best_distance = distance;
      best_cluster = cluster_info.cluster;
    } else if (distance < second_best_distance) {
      second_best_distance = distance;
    }
  }

  // Assign to the best cluster if unambiguous.
  if (best_cluster != nullptr && best_distance != second_best_distance) {
    best_cluster->alignments.emplace_back(partial_read);
    partial_read->cluster = best_cluster.get();
  }
}

std::pair<UnassignedPartialReads, UnassignedPartialReads> AssignPartialReadsToClusters(
    const ReadCollapserOptions& options,
    const bool cluster_by_strand,
    const UmiClusterPtrAndStats& cluster_stats_by_start,
    const UmiClusterPtrAndStats& cluster_stats_by_end,
    vec<AlignmentPtr>& unassigned_partials) {
  UnassignedPartialReads still_unassigned_3p_partials;
  UnassignedPartialReads still_unassigned_5p_partials;
  for (auto& partial_read : unassigned_partials) {
    const bool is_5p_complete = partial_read->IsFivePrimeComplete();
    const read_name::Umi umi = is_5p_complete ? partial_read->umi5p : partial_read->umi3p;
    if (cluster_stats_by_start.contains(umi) && is_5p_complete) {
      AssignPartialReadToCluster(
          options.wiggle_room_partial, cluster_by_strand, cluster_stats_by_start.at(umi), partial_read);
    } else if (cluster_stats_by_end.contains(umi) && !is_5p_complete) {
      AssignPartialReadToCluster(
          options.wiggle_room_partial, cluster_by_strand, cluster_stats_by_end.at(umi), partial_read);
    } else {
      // No candidate clusters with matching UMI, will consider it later
      // when we optionally cluster remaining unassigned partial reads into clusters based on position.
    }

    if (partial_read->cluster == nullptr) {
      const auto& strand = DetermineClusterStrand(partial_read, cluster_by_strand);
      if (is_5p_complete) {
        still_unassigned_3p_partials[{umi, strand}].emplace_back(partial_read);
      } else {
        still_unassigned_5p_partials[{umi, strand}].emplace_back(partial_read);
      }
    }
  }
  return {still_unassigned_3p_partials, still_unassigned_5p_partials};
}

/**
 * @brief Clusters partial reads based on alignment position. Reads
 *        should already be grouped by UMI and strand (if applicable).
 *
 * This function processes a collection of partial reads, where each alignment
 * contains only one UMI (either 5' or 3'). It clusters these reads based on
 * their positional proximity, considering a specified wiggle room. The reads
 * are grouped into clusters, which are stored in the provided `clusters` container.
 * The clustering process depends on whether the UMI is 5' or 3', as the alignments
 * are sorted differently in each case.
 *
 * @param wiggle_room The maximum positional difference allowed between alignments
 *                    to be included in the same cluster.
 * @param create_cluster_id A callable that generates unique cluster IDs.
 * @param partial_reads A vector of read pointers to be clustered. For
 *                      partial reads with 3' UMI only, these reads will be sorted by
 *                      their end positions; for partial reads with 5' UMI only, they are
 *                      already sorted by their start positions.
 * @param is_5p_complete A boolean indicating whether the partial read is complete and untruncated on the 5' end.
 *                       If false, the read is assumed to be complete on the 3' end.
 * @param clusters A reference to a container where the resulting clusters will
 *                 be stored. Each cluster is identified by its unique cluster ID.
 *                 It is updated in-place as clusters are created from the provided partial reads.
 *
 * @note This function assumes that the input reads are already grouped by UMI
 *       and strand (if applicable). If the input `partial_reads` is empty, the
 *       function returns early without performing any clustering.
 */
void ClusterUnassignedPartialReadsByPosition(const u32 wiggle_room,
                                             const CreateClusterId& create_cluster_id,
                                             vec<AlignmentPtr>& partial_reads,
                                             const bool is_5p_complete,
                                             Clusters& clusters) {
  // The following code assumes there is at least one read, return early if there are none.
  if (partial_reads.empty()) {
    return;
  }

  // The input reads must be sorted by EndPos() if it is incomplete (partial) on the 5' end
  if (!is_5p_complete) {
    std::ranges::sort(partial_reads, [](const auto& a, const auto& b) { return a->EndPos() < b->EndPos(); });
  }

  auto current_cluster = std::make_shared<Cluster>(create_cluster_id(), vec<AlignmentPtr>{});

  // Start with the first read.
  auto& first_read = partial_reads.front();
  first_read->cluster = current_cluster.get();
  current_cluster->alignments.emplace_back(first_read);
  u32 last_cluster_pos = is_5p_complete ? first_read->StartPos() : first_read->EndPos();

  // Iterate through the remaining reads and cluster them based on positional proximity.
  const auto remaining_reads = partial_reads | std::views::drop(1);

  for (const auto& partial_read : remaining_reads) {
    const u32 current_pos = is_5p_complete ? partial_read->StartPos() : partial_read->EndPos();

    // Check positional proximity to the last read added.
    if (current_pos <= last_cluster_pos + wiggle_room) {
      // Merge into current cluster
      partial_read->cluster = current_cluster.get();
      current_cluster->alignments.emplace_back(partial_read);
      // Update the pivot position to the new read's position.
      last_cluster_pos = current_pos;
    } else {
      // Finalize the current cluster, store it, and start a new one.
      clusters[current_cluster->cluster_id] = current_cluster;

      current_cluster = std::make_shared<Cluster>(create_cluster_id(), vec<AlignmentPtr>{});
      partial_read->cluster = current_cluster.get();
      current_cluster->alignments.emplace_back(partial_read);
      last_cluster_pos = current_pos;
    }
  }

  // Add the final cluster.
  if (!current_cluster->alignments.empty()) {
    clusters[current_cluster->cluster_id] = current_cluster;
  }
}

/**
 * @brief Calculates the mean, minimum, and maximum positions for each cluster (position-only, no UMI grouping).
 *
 * This function iterates through the provided clusters and computes the mean start and end positions, as well as the
 * minimum start and maximum end positions for each cluster. The results are stored in vectors keyed by the cluster's
 * relevant position (start for 3' partial assignment, end for 5' partial assignment).
 *
 * @param clusters A collection of clusters, where each cluster contains a set of alignments.
 * @param cluster_by_strand Whether to cluster by strand orientation.
 * @param cluster_stats_by_start Output vector of cluster statistics for assigning 3'-partial reads (reads with
 * truncated 3' ends) keyed by mean start.
 * @param cluster_stats_by_end Output vector of cluster statistics for assigning 5'-partial reads (reads with truncated
 * 5' ends) keyed by mean end.
 */
void CalculateClusterPositionStats(const Clusters& clusters,
                                   const bool cluster_by_strand,
                                   UmiClusterPtrAndStats& cluster_stats_by_start,
                                   UmiClusterPtrAndStats& cluster_stats_by_end) {
  for (const auto& [cluster_id, cluster_ptr] : clusters) {
    const read_name::Umi umi5p = cluster_ptr->alignments.front()->umi5p;
    const read_name::Umi umi3p = cluster_ptr->alignments.front()->umi3p;
    auto strand = DetermineClusterStrand(cluster_ptr->alignments.front(), cluster_by_strand);

    f64 mean_start = 0;
    f64 mean_end = 0;
    auto min_start = std::numeric_limits<u32>::max();
    auto max_end = std::numeric_limits<u32>::min();
    for (const auto& read : cluster_ptr->alignments) {
      mean_start += read->StartPos();
      mean_end += read->EndPos();
      min_start = std::min(min_start, read->StartPos());
      max_end = std::max(max_end, read->EndPos());
    }
    mean_start /= static_cast<f64>(cluster_ptr->alignments.size());
    mean_end /= static_cast<f64>(cluster_ptr->alignments.size());

    cluster_stats_by_start[umi5p].emplace_back(cluster_ptr, strand, mean_start, min_start, max_end);
    cluster_stats_by_end[umi3p].emplace_back(cluster_ptr, strand, mean_end, min_start, max_end);
  }
}

}  // namespace

/**
 * Cluster reads based on their alignment positions and optionally their UMIs. The clustering process consists of
 * several phases:
 *   - Phase 1: Separate reads into full, 5'-partial, and 3'-partial groups.
 *   - Phase 2: Cluster full reads by UMI if available, otherwise by position within the wiggle room distance.
 *   - Phase 3: If exclude_partial_reads, skip partial processing.
 *   - Phase 4: Calculate cluster position statistics of existing full clusters for partial assignment.
 *   - Phase 5: Assign partial reads to nearest full cluster within wiggle_room_partial.
 *   - Phase 6: Optionally cluster remaining unassigned partial reads into partial-only clusters.
 * Note: @p reads must be sorted by start position. UMIs must be pre-populated.
 * @param options Read collapser options.
 * @param create_cluster_id Callback function to create a unique cluster ID.
 * @param reads Vector of reads to cluster. This vector is modified in place to assign cluster pointers to each read.
 *              Although the reference to the vector itself is const, the alignments within the vector are modified by
 *              assigning cluster information to them. Any unwanted alignments (e.g. missing UMI,
 *              secondary/supplementary alignments) should be filtered out before calling this function.
 * @param clusters Output map of cluster IDs to cluster pointers. This map is modified in place to add the resulting
 *                 clusters from the clustering process.
 */
void ClusterReads(const ReadCollapserOptions& options,
                  const CreateClusterId& create_cluster_id,
                  const vec<AlignmentPtr>& reads,
                  Clusters& clusters,
                  Metrics& metrics) {
  if (reads.empty()) {
    return;
  }

  const bool cluster_by_umi = options.cluster_by_umi;

  // Phase 1: Separate reads into full and partial reads.
  vec<AlignmentPtr> full_reads;
  vec<AlignmentPtr> partial_reads;
  for (const auto& read : reads) {
    if (read->IsPartial()) {
      partial_reads.emplace_back(read);
    } else {
      full_reads.emplace_back(read);
    }
  }

  // Phase 2: Cluster full reads by UMI if available, otherwise by position.
  if (cluster_by_umi) {
    Clusters umi_clusters;
    ClusterFullReadsByUmi(create_cluster_id, full_reads, umi_clusters);
    // Cluster the resulting UMI clusters by position to further split clusters with the same UMI that are beyond the
    // wiggle room distance.
    for (const auto& [cluster_id, cluster_ptr] : umi_clusters) {
      ClusterFullReadsByPosition(
          options.wiggle_room, options.cluster_by_strand, create_cluster_id, cluster_ptr->alignments, clusters);
    }
  } else {
    ClusterFullReadsByPosition(options.wiggle_room, options.cluster_by_strand, create_cluster_id, full_reads, clusters);
  }

  // Phase 3: If partial reads are excluded, update metrics and return.
  if (options.exclude_partial_reads) {
    const size_t num_unclustered = partial_reads.size();
    metrics.unclustered_partial_reads += num_unclustered;
    metrics.clustering_input_reads -= num_unclustered;
    return;
  }

  // If there are no partial reads, we are done.
  if (partial_reads.empty()) {
    return;
  }

  // Phase 4: Calculate cluster position statistics for partial assignment.
  // These statistics are grouped by UMI (if any). If we are not clustering by UMI, all clusters
  // will have '*' as their UMI key so that all clusters are considered together for partial assignment.
  UmiClusterPtrAndStats cluster_stats_by_start_umi;
  UmiClusterPtrAndStats cluster_stats_by_end_umi;
  CalculateClusterPositionStats(
      clusters, options.cluster_by_strand, cluster_stats_by_start_umi, cluster_stats_by_end_umi);

  // Phase 5: Assign partial reads to nearest full cluster by position.
  // 5'-partial reads is complete on the 3' end, so we match by end position.
  // 3'-partial reads is complete on the 5' end, so we match by start position.
  auto [unassigned_3p_partials, unassigned_5p_partials] = AssignPartialReadsToClusters(
      options, options.cluster_by_strand, cluster_stats_by_start_umi, cluster_stats_by_end_umi, partial_reads);

  // Phase 6: Cluster remaining unassigned partial reads or leave them unclustered.
  if (!options.make_clusters_of_partial_reads_only) {
    size_t num_unclustered = 0;
    for (const auto& unassigned_partials : unassigned_3p_partials | std::views::values) {
      num_unclustered += unassigned_partials.size();
    }
    for (const auto& unassigned_partials : unassigned_5p_partials | std::views::values) {
      num_unclustered += unassigned_partials.size();
    }
    metrics.clustering_unclustered_partial_reads += num_unclustered;
    // Although these reads are not clustered, they are still considered during the clustering step,
    // so we include them in the clustering reads metric.
    metrics.clustering_reads += num_unclustered;
    metrics.clustering_partial_reads += num_unclustered;
    // Partial reads that are not made into clusters are still included in the vector, but are not clustered.
    return;
  }

  // Cluster remaining unassigned partial reads if specified by options.
  // These reads are already grouped by UMI and strand (if applicable).
  for (auto& [umi_strand, unassigned_partial_reads] : unassigned_3p_partials) {
    // 3' partial reads are complete on the 5' end, so we cluster by start position
    ClusterUnassignedPartialReadsByPosition(
        options.wiggle_room_partial, create_cluster_id, unassigned_partial_reads, true, clusters);
  }
  for (auto& [umi_strand, unassigned_partial_reads] : unassigned_5p_partials) {
    // 5' partial reads are complete on the 3' end, so we cluster by end position
    ClusterUnassignedPartialReadsByPosition(
        options.wiggle_room_partial, create_cluster_id, unassigned_partial_reads, false, clusters);
  }
}

/**
 * Cluster alignments and update clustering metrics.
 *
 * This function modifies the `alignments` vector by assigning each read in the vector
 * a pointer to the cluster it belongs to. Unclustered reads will have a `nullptr` cluster pointer.
 * Any supplementary alignments and non-rescued secondary alignments are ignored during clustering.
 *
 * Unmapped reads are not clustered either as they should be handled separately by
 * `ReadAndWriteUnmappedReads`.
 *
 * @param options Read collapser options.
 * @param alignments Vector of alignments to cluster. Primary alignments and rescued secondary alignments may be
 * clustered.
 * @param create_cluster_id Callback function to create a unique cluster ID.
 *
 * @return Clusters containing the clustered reads.
 */
Clusters ClusterAlignments(const ReadCollapserOptions& options,
                           vec<AlignmentPtr>& alignments,
                           const CreateClusterId& create_cluster_id,
                           Metrics& metrics) {
  vec<AlignmentPtr> primary_alignments;
  primary_alignments.reserve(alignments.size());
  for (const auto& alignment : alignments) {
    const auto flag = alignment->record->core.flag;
    if ((flag & BAM_FUNMAP) != 0) {
      alignment->cluster = nullptr;
      ++metrics.unmapped_reads;
    } else if ((flag & BAM_FSUPPLEMENTARY) != 0) {
      alignment->cluster = nullptr;
      ++metrics.unclustered_supplementary_alignments;
    } else if ((flag & BAM_FSECONDARY) != 0) {
      if (options.markdup_options.cluster_rescued_secondaries && HasRescuedSecondaryTag(alignment->record.get())) {
        primary_alignments.emplace_back(alignment);
        ++metrics.rescued_secondary_alignments;
      } else {
        alignment->cluster = nullptr;
        ++metrics.unclustered_secondary_alignments;
      }
    } else {
      primary_alignments.emplace_back(alignment);
    }
  }

  // Clustering input reads are all clustering-eligible alignments that were considered for clustering, which includes
  // clustered reads and unclustered partial reads (clustering_unclustered_partial_reads), but excludes supplementary
  // reads and non-rescued secondary reads, as well as discarded reads with missing UMIs.
  metrics.clustering_input_reads += primary_alignments.size();
  Clusters clusters;

  // Importantly, the alignments vector is modified in place to assign each read a pointer to the cluster it belongs to.
  // While primary_alignments is a filtered view of alignments, the AlignmentPtr objects within it are the same as those
  // in alignments, so modifying their cluster pointers also modifies those in alignments. This ensures that after
  // clustering, the original alignments vector, containing ALL reads (clustered, unclustered, supplementary,
  // secondary), has the correct cluster pointers assigned and is preserved for I/O operations.
  ClusterReads(options, create_cluster_id, primary_alignments, clusters, metrics);

  // Duplicate reads are based on pre-duplex deconvolution cluster sizes.
  UpdateDuplicateReadMetrics(clusters, metrics);

  // Deconvolve duplex reads if the option is enabled
  if (options.consensus_options.duplex_library_type != HDDeconvolutionType::kNone) {
    for (const auto& [cluster_id, cluster] : clusters) {
      const bool is_parent_parent =
          (options.consensus_options.duplex_library_type == HDDeconvolutionType::kParentParent);
      cluster->alignments = DeconvolveDuplexReads(cluster->alignments, is_parent_parent);
    }
  }

  UpdateClusterMetrics(clusters, metrics);
  return clusters;
}

}  // namespace xoos::read_collapser
