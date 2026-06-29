#include "clustering/clustering.h"

#include <iterator>
#include <ranges>
#include <stack>

#include <htslib/sam.h>

#include <xoos/error/error.h>
#include <xoos/util/hash.h>
#include <xoos/util/math.h>
#include <xoos/util/sequence-functions.h>

#include "metrics/metrics.h"
#include "util/duplex-util.h"

namespace xoos::read_collapser {

size_t ClusterIdHash::operator()(const ClusterId& cluster_id) const {
  return util::hash::Hash(cluster_id.super_region_id, cluster_id.count);
}

size_t UmiPairHash::operator()(const UmiPair& umi_pair) const {
  return util::hash::Hash(umi_pair.umi3p, umi_pair.umi5p);
}

ClusterStrand DetermineClusterStrand(const AlignmentPtr& alignment, const bool cluster_by_strand) {
  using enum ClusterStrand;
  if (!cluster_by_strand) {
    return kBoth;
  }
  return (alignment->record->core.flag & BAM_FREVERSE) != 0 ? kRev : kFwd;
}

size_t UmiStrandHash::operator()(const UmiStrand& umi_strand) const {
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
static void UpdateDuplicateReadMetrics(const Clusters& clusters) {
  auto& metrics = ConcurrentMetrics::Get();

  for (const auto& cluster : clusters | std::views::values) {
    const auto cluster_size = cluster->alignments.size();
    // All but one read in a cluster are duplicates (pre-duplex deconvolution).
    metrics.duplicate_reads += cluster_size - 1;
  }
}

// Update the cluster metrics based on the clusters provided.
void UpdateClusterMetrics(const Clusters& clusters) {
  auto& metrics = ConcurrentMetrics::Get();

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
    for (const auto& alignment : cluster->alignments) {
      metrics.UpdateClusterSizeHistogram(cluster_size, "total_reads", 1);
      ++metrics.clustering_reads;
      if (alignment->IsPartial()) {
        cluster_is_partial = true;
        ++metrics.clustering_partial_reads;
      } else {
        cluster_is_full = true;
        ++metrics.clustering_full_reads;
      }
      if (alignment->IsForward()) {
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

void DepthFirstClusterSearch(const ClusterCoord& initial_coord,
                             const PartialCluster& partial_cluster,
                             PartialClusters& partial_clusters,
                             u32 wiggle_room,
                             u32 max_depth) {
  // Recursive depth first search is implemented using a stack rather than actual function recursion to reduce
  // the risk of stack overflow.
  std::stack<ClusterCoord> coords;
  coords.emplace(initial_coord);

  while (!coords.empty() && (max_depth == 0 || coords.size() < max_depth)) {
    auto [start, end, strand] = coords.top();
    coords.pop();

    u32 neighbor_start_lower = math::SatSub(start, wiggle_room);
    u32 neighbor_start_upper = start + wiggle_room;

    // Perform a walk of all potential neighbors of the current partial cluster. To do this we generate all
    // cluster coordinates which meet the constraint |a.start - b.start| + |a.end - b.end| <= wiggle_room, and
    // search for partial clusters at these coordinates. If any are found, we merge them into the current cluster
    // and add their coordinates to the stack for further searching.
    for (u32 neighbor_start = neighbor_start_lower, neighbor_end_lower = end, neighbor_end_upper = end;
         neighbor_start <= neighbor_start_upper;
         neighbor_start++) {
      for (u32 neighbor_end = neighbor_end_lower; neighbor_end <= neighbor_end_upper; neighbor_end++) {
        ClusterCoord neighbor_coord = {neighbor_start, neighbor_end, strand};

        // If there is no partial cluster at this coordinate, or if the partial cluster has already been merged into
        // a cluster, skip this neighbor.
        auto neighbor_it = partial_clusters.find(neighbor_coord);
        if (neighbor_it == partial_clusters.end() || neighbor_it->second.cluster != nullptr) {
          continue;
        }

        // We have found a neighboring partial cluster, merge it into the current cluster.
        neighbor_it->second.cluster = partial_cluster.cluster;
        for (auto& alignment : neighbor_it->second.alignments) {
          alignment->cluster = partial_cluster.cluster.get();
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
 *
 * Note: The range defined by [@p cluster_start_index, @p cluster_end_index) is inclusive-exclusive.
 */
void ClusterAlignmentsByPosition(u32 wiggle_room,
                                 const bool cluster_by_strand,
                                 const CreateClusterId& create_cluster_id,
                                 Clusters& clusters,
                                 const std::forward_iterator auto& cluster_start,
                                 const std::forward_iterator auto& cluster_end) {
  // Iterate through all the alignments that we should  consider for clustering, if the ClusterCoord for
  // the alignment is not in the cluster_map, create a new cluster and add it to the cluster_map. Otherwise,
  // add the alignment to the existing cluster. This
  PartialClusters cluster_map;
  for (auto alignment_it = cluster_start; alignment_it != cluster_end; ++alignment_it) {
    ClusterCoord coord = CreateClusterCoord(*alignment_it, cluster_by_strand);
    if (auto cluster = cluster_map.find(coord); cluster == cluster_map.end()) {
      cluster_map.emplace_hint(cluster, coord, PartialCluster{nullptr, {*alignment_it}});
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
    for (auto& alignment : cluster_tmp.alignments) {
      alignment->cluster = cluster_tmp.cluster.get();
      alignment->cluster->alignments.emplace_back(alignment);
    }

    DepthFirstClusterSearch(coords, cluster_tmp, cluster_map, wiggle_room, 0);
    clusters[cluster_id] = cluster_tmp.cluster;
  }
}

/**
 * Positional clustering is broken up into two phases:
 *   - Phase 1: Find all alignments whose start position is within the wiggle room of the previous alignment.
 *   - Phase 2: For each set of alignments found in phase 1, determine the clusters they belong to based on their
 *   start, end, and strand.
 * Note: @p alignments must be sorted by start position.
 */
void ClusterAlignmentsByPosition(u32 wiggle_room,
                                 const bool cluster_by_strand,
                                 const CreateClusterId& create_cluster_id,
                                 const vec<AlignmentPtr>& alignments,
                                 Clusters& clusters) {
  // The following code assumes there is at least one alignment, return early if there are none.
  if (alignments.empty()) {
    return;
  }

  // The start of the first cluster is the first alignment, this will be updated for each new cluster we start.
  auto cluster_start_it = alignments.begin();

  for (auto alignment_it = std::next(cluster_start_it); alignment_it != alignments.end(); ++alignment_it) {
    // Check if the alignment is part of the current cluster by checking if it's start position is within the wiggle
    // room of the previous alignments start position. If it is, move onto the next alignment.
    auto prev_alignment = std::prev(alignment_it);
    if ((*alignment_it)->record->core.pos <= (*prev_alignment)->record->core.pos + wiggle_room) {
      // This alignment is part of the current cluster, move onto the next alignment.
      continue;
    }

    // The current alignment is not part of the cluster, so cluster the alignments from the start index up to the
    // current.
    ClusterAlignmentsByPosition(
        wiggle_room, cluster_by_strand, create_cluster_id, clusters, cluster_start_it, alignment_it);

    // Update the start of the next cluster to be the current alignment.
    cluster_start_it = alignment_it;
  }

  if (cluster_start_it != alignments.end()) {
    // Handle the last cluster if it was not previously handled and there is at least one alignment.
    ClusterAlignmentsByPosition(
        wiggle_room, cluster_by_strand, create_cluster_id, clusters, cluster_start_it, alignments.end());
  }
}

/**
 * Check if the UMI is a sequence containing only ACGT characters.
 * @param umi The UMI to check.
 */
static bool IsSequenceUmi(const std::optional<std::string>& umi) {
  return umi.has_value() && umi->find_first_not_of("ACGT") == std::string::npos;
}

/**
 * Determine if the UMI format used in the qname is the newer format, which uses ':' as a delimiter for UMIs, or the
 * legacy format, which uses '|' as a delimiter for UMIs.  Both use '|' to separate the read info and sid from the UMIs.
 *
 * @param qname The query name string to check for UMI format. It is expected to contain at least one '|' character
 * separating the read info and sid from the UMIs.
 * @return true if the qname uses the newer UMI format (':' delimiter), false if it uses the legacy UMI format ('|'
 * delimiter).
 * @throw  error::Error if the qname does not contain at least one '|' character, which is expected to separate the read
 * info and sid from the UMIs.
 */
static bool UmiFormatIsNew(const std::string_view& qname) {
  const size_t last_pipe_pos = qname.rfind('|');
  if (last_pipe_pos == std::string_view::npos) {
    throw error::Error("{} does not contain '|'", qname);
  }

  // Check if there's a second pipe before the last one
  const size_t second_last_pipe_pos = qname.rfind('|', last_pipe_pos - 1);

  // If there's no second pipe, it could be the new format (uses ':' as delimiter)
  // If there are 2+ pipes, it's the legacy format (uses '|' as delimiter)
  return (second_last_pipe_pos == std::string_view::npos);
}

/**
 * @brief Parses a query name (qname) string to extract two UMI (Unique Molecular Identifier) values.
 *
 * This function supports two formats:
 * - Legacy: <read_info>|<sid>|<umi5p>|<umi3p>
 * - Newer:  <read_info>|<sid>:<umi5p>:<umi3p>
 *
 * This function handles the following scenarios:
 * 1. UMIs Present: Extracts the 5' and 3' UMI sequences.
 * 2. UMIs Missing: If either UMI is missing (represented by `kNoUmi`), it is returned as `std::nullopt`.
 * 3. UMIs Omitted: If the dataset does not support UMIs and the fields are omitted entirely
 * (e.g., `<read_info>|<sid>`), this function will throw an error because the expected delimiters
 * are not found and the function expects to only be called for read names with UMIs.
 *
 * @param qname The input string containing the query name with embedded UMIs, separated by '|'.
 * @return A tuple containing two optional UMIs:
 * - The first UMI (5' UMI) as `std::optional<Umi>`.
 * - The second UMI (3' UMI) as `std::optional<Umi>`.
 * @throws error::Error If read name is ill-formed or does not contain the expected number of UMI delimiters.
 */
std::tuple<Umi, Umi> ParseUmi(const std::string_view& qname) {
  // Detect format by counting pipes
  // Legacy format: <read_info>|<sid>|<umi5p>|<umi3p>
  // Newer format: <read_info>|<sid>:<umi5p>:<umi3p>
  // Newer format no UMIs: <read_info>|<sid>
  if (const size_t last_pipe_pos = qname.rfind('|'); last_pipe_pos == std::string_view::npos) {
    throw error::Error("{} does not contain '|'", qname);
  }

  // if its the newer UMI format, then we expect the UMIs to be separated by ':', if its the legacy format we expect the
  // UMIs to be separated by '|'
  const char umi_delimiter = UmiFormatIsNew(qname) ? ':' : '|';

  // Find the last occurrence of the delimiter within the UMI section
  const size_t last_umi_delim_pos = qname.rfind(umi_delimiter);

  // If the delimiter is not found, it means there are no UMIs in the input
  if (last_umi_delim_pos == std::string_view::npos) {
    throw error::Error("'{}' does not contain additional expected UMI delimiters '{}'", qname, umi_delimiter);
  }

  // Find the second-to-last occurrence of the delimiter within the UMI section
  const size_t second_last_umi_delim_pos = qname.rfind(umi_delimiter, last_umi_delim_pos - 1);
  if (second_last_umi_delim_pos == std::string_view::npos) {
    throw error::Error("'{}' does not contain two UMI delimiters '{}'", qname, umi_delimiter);
  }

  // Extract UMI values
  auto umi5p =
      std::string(qname.substr(second_last_umi_delim_pos + 1, last_umi_delim_pos - second_last_umi_delim_pos - 1));
  auto umi3p = std::string(qname.substr(last_umi_delim_pos + 1));

  // Convert kNoUmi to nullopt
  return {
      (umi5p == kNoUmi) ? std::nullopt : std::make_optional(umi5p),
      (umi3p == kNoUmi) ? std::nullopt : std::make_optional(umi3p),
  };
}

/**
 * @brief Parses the UMI (Unique Molecular Identifier) from the read names of the given alignments
 *        and updates the alignment objects with the parsed UMI values.
 *
 * This function processes a collection of alignment objects, extracts the UMI from the read name
 * (expected to be at the end of the qname, separated by '|'), and updates the alignment objects
 * with the parsed UMI values. If the alignment is reverse, the 5' and 3' UMIs are swapped, and
 * sequence UMIs are reverse complemented. Alignments without a valid UMI are discarded, and
 * corresponding metrics are updated.
 *
 * @param alignments A vector of alignment pointers to be processed.
 * @param metrics A reference to a Metrics object for tracking discarded reads.
 *
 * @throws error::Error If the UMI cannot be parsed from the read name.
 */
static void ParseUmis(vec<AlignmentPtr>& alignments, Metrics& metrics) {
  vec<AlignmentPtr> umi_alignments;
  umi_alignments.reserve(alignments.size());
  for (const auto& alignment : alignments) {
    std::string_view qname{bam_get_qname(alignment->record.get())};
    // Parse the UMI from the read name. The UMI is expected to be at the end of the qname, separated by '|'.
    auto [umi5p, umi3p] = ParseUmi(qname);
    if (!umi5p.has_value() && !umi3p.has_value()) {
      // If both UMIs are missing, do not include this read in umi_alignments.
      ++metrics.discarded_missing_umi_reads;
      ++metrics.discarded_total_reads;
      continue;
    }
    // During alignment the direction of the read is determined, but the UMIs are not updated.
    // If the alignment is reverse, we swap the 5' and 3' UMIs. If the UMI are sequences,
    // we also need to reverse complement them.
    if (alignment->IsReverse()) {
      umi5p = IsSequenceUmi(umi5p) ? sequence::ReverseComplement(umi5p.value()) : umi5p;
      umi3p = IsSequenceUmi(umi3p) ? sequence::ReverseComplement(umi3p.value()) : umi3p;
      std::swap(umi5p, umi3p);
    }
    alignment->umi5p = umi5p;
    alignment->umi3p = umi3p;
    umi_alignments.emplace_back(alignment);
  }
  // Replace the original alignments vector with the valid UMI alignments.
  alignments = std::move(umi_alignments);
}

/**
 * @brief Clusters alignments based on their UMI (Unique Molecular Identifier) pairs.
 *
 * This function processes a collection of alignments and groups them into clusters
 * based on their UMI pairs. Each UMI pair is represented as a combination of 5' and 3' UMIs.
 * If a cluster for a given UMI pair does not exist, a new cluster is created. If a cluster
 * already exists for the UMI pair, the alignment is added to the existing cluster.
 *
 * @param create_cluster_id A callable object or function that generates unique cluster IDs.
 * @param alignments A vector of alignment pointers to be clustered. Each alignment must
 *                   contain valid 5' and 3' UMI values.
 * @param clusters A map that stores the resulting clusters, keyed by their unique cluster IDs.
 *
 * @note This function assumes that the `alignments` vector is non-empty. If it is empty,
 *       the function returns early without performing any operations.
 * @note The `clusters` map is updated in-place with the newly created or updated clusters.
 */
static void ClusterFullAlignmentsByUmi(const CreateClusterId& create_cluster_id,
                                       vec<AlignmentPtr>& alignments,
                                       Clusters& clusters) {
  // The following code assumes there is at least one alignment, return early if there are none.
  if (alignments.empty()) {
    return;
  }
  // Create a map of UMIs to clusters to determine if a cluster already exists for a given UMI pair.
  std::unordered_map<UmiPair, Cluster*, UmiPairHash> umi_cluster_map;
  for (const auto& alignment : alignments) {
    auto umi_pair = UmiPair{alignment->umi5p.value(), alignment->umi3p.value()};
    auto it = umi_cluster_map.find(umi_pair);
    // UMI pair does not exist in the map, create a new cluster.
    if (it == umi_cluster_map.end()) {
      auto new_cluster_id = create_cluster_id();
      auto new_cluster = std::make_shared<Cluster>(new_cluster_id, vec<AlignmentPtr>{});
      alignment->cluster = new_cluster.get();
      alignment->cluster->alignments.emplace_back(alignment);
      clusters[new_cluster_id] = new_cluster;
      umi_cluster_map[umi_pair] = new_cluster.get();
    } else {
      // Cluster exists for UMI pair and must be updated with alignment.
      // Assumes cluster already contains at least one alignment (from cluster creation step).
      alignment->cluster = it->second;
      alignment->cluster->alignments.emplace_back(alignment);
    }
  }
}

/**
 * @brief Assigns a partial alignment to a cluster based on matching UMI and positional proximity.
 *
 * This function attempts to assign an alignment to a cluster by comparing the alignment's UMI and position
 * to existing clusters. It considers strand separation, positional proximity, and cluster boundaries.
 *
 * @param wiggle_room The maximum allowable positional deviation for an alignment to be considered part of a cluster.
 * @param cluster_by_strand Specifies whether alignments should be separated by strand.
 * @param umi_to_cluster_info A map containing clusters and their associated positions.
 * @param alignment The alignment to be assigned to a cluster.
 * @param umi5p A boolean indicating whether to use the 5' UMI (true) or the 3' UMI (false) for clustering.
 */
static void AssignPartialAlignmentToCluster(const u32 wiggle_room,
                                            const bool cluster_by_strand,
                                            UmiClusterPtrAndStats& umi_to_cluster_info,
                                            AlignmentPtr& alignment,
                                            const bool umi5p) {
  f64 best_distance = std::numeric_limits<f64>::max();
  ClusterPtr best_cluster = nullptr;
  // Keep track of the second best distance to determine if any cluster assignments are ambiguous (i.e. equidistant).
  f64 second_best_distance = std::numeric_limits<f64>::max();
  auto umi = umi5p ? alignment->umi5p.value() : alignment->umi3p.value();
  // If the alignment's UMI is not in the map it cannot be assigned to a full UMI cluster.
  if (!umi_to_cluster_info.contains(umi)) {
    return;
  }

  const u32 pos = umi5p ? alignment->StartPos() : alignment->EndPos();

  for (const auto& cluster_info : umi_to_cluster_info[umi]) {
    // Skip alignments on different strands only when the cluster_by_strand is set to true.
    if (const auto strand = DetermineClusterStrand(alignment, cluster_by_strand);
        cluster_by_strand && cluster_info.strand != strand) {
      continue;
    }

    f64 distance = std::abs(static_cast<f64>(pos) - cluster_info.mean_pos);
    if (distance > wiggle_room) {
      continue;
    }

    // Check if the alignment is too long for the cluster.
    if ((umi5p && alignment->EndPos() > cluster_info.max_end + wiggle_room) ||
        (!umi5p && alignment->StartPos() < cluster_info.min_start - wiggle_room)) {
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

  // If a best cluster has been identified, assign the alignment to the cluster.
  // If the best distance is the same as the second best distance, the alignment is equidistant from two clusters and
  // is therefore ambiguous.
  if (best_cluster != nullptr && best_distance != second_best_distance) {
    best_cluster->alignments.emplace_back(alignment);
    alignment->cluster = best_cluster.get();
  }
}

/**
 * @brief Assigns partial UMI alignments to clusters or classifies them as unassigned.
 *
 * This function processes a collection of partial alignments and attempts to assign
 * each alignment to a cluster. If an alignment cannot be assigned to a cluster, it is added to the unassigned partial
 * alignments map for further processing.
 *
 * @param options The configuration options for the read collapser, including parameters
 *                for clustering and wiggle room for partial alignments.
 * @param partial_alignments A vector of pointers to partial alignments to be processed.
 * @param umi_to_cluster_info A reference to the UMI cluster information, which maps UMIs
 *                            to their corresponding clusters and position information.
 * @param unassigned_partial_alignments A map that stores unassigned partial alignments
 *                                       grouped by their UMI and strand.
 * @param umi5p A boolean flag indicating whether to use the 5' UMI (true) or the 3' UMI
 *              (false) for clustering.
 */
static void AssignPartialAlignmentsToClusters(const ReadCollapserOptions& options,
                                              vec<AlignmentPtr>& partial_alignments,
                                              UmiClusterPtrAndStats& umi_to_cluster_info,
                                              UnassignedPartialAlignments& unassigned_partial_alignments,
                                              const bool umi5p) {
  for (auto& alignment : partial_alignments) {
    AssignPartialAlignmentToCluster(
        options.wiggle_room_partial, options.cluster_by_strand, umi_to_cluster_info, alignment, umi5p);
    if (alignment->cluster == nullptr) {
      const auto& strand = DetermineClusterStrand(alignment, options.cluster_by_strand);
      const auto& umi = umi5p ? alignment->umi5p.value() : alignment->umi3p.value();
      unassigned_partial_alignments[{umi, strand}].emplace_back(alignment);
    }
  }
}

/**
 * @brief Calculates the mean, minimum, and maximum positions for each cluster
 *        and organizes them into UMI-based cluster statistics.
 *
 * This function iterates through the provided clusters, computes the mean start
 * and end positions, as well as the minimum start and maximum end positions for
 * each cluster. It then categorizes the clusters based on their 5' and 3' UMI
 * values, storing the computed statistics in the provided UMI cluster maps.
 *
 * @param clusters A collection of clusters, where each cluster contains a set
 *                 of alignments.
 * @param cluster_by_strand A flag indicating whether to cluster by
 *                           strand orientation.
 * @param umi5p_clusters A map to store the 5' UMI-based cluster statistics. Each
 *                       entry contains the cluster pointer, strand, mean start
 *                       position, minimum start position, and maximum end position.
 * @param umi3p_clusters A map to store the 3' UMI-based cluster statistics. Each
 *                       entry contains the cluster pointer, strand, mean end
 *                       position, minimum start position, and maximum end position.
 */
static void CalculateClusterMeanMinMaxPosition(const Clusters& clusters,
                                               const bool cluster_by_strand,
                                               UmiClusterPtrAndStats& umi5p_clusters,
                                               UmiClusterPtrAndStats& umi3p_clusters) {
  for (const auto& [cluster_id, cluster_ptr] : clusters) {
    auto umi5p = cluster_ptr->alignments.front()->umi5p.value();
    auto umi3p = cluster_ptr->alignments.front()->umi3p.value();
    auto strand = DetermineClusterStrand(cluster_ptr->alignments.front(), cluster_by_strand);

    f64 mean_start = 0;
    f64 mean_end = 0;
    auto min_start = std::numeric_limits<u32>::max();
    auto max_end = std::numeric_limits<u32>::min();
    for (const auto& alignment : cluster_ptr->alignments) {
      mean_start += alignment->StartPos();
      mean_end += alignment->EndPos();
      min_start = std::min(min_start, alignment->StartPos());
      max_end = std::max(max_end, alignment->EndPos());
    }
    mean_start /= static_cast<f64>(cluster_ptr->alignments.size());
    mean_end /= static_cast<f64>(cluster_ptr->alignments.size());

    umi5p_clusters[umi5p].emplace_back(cluster_ptr, strand, mean_start, min_start, max_end);
    umi3p_clusters[umi3p].emplace_back(cluster_ptr, strand, mean_end, min_start, max_end);
  }
}

/**
 * @brief Clusters partial alignments based on position. Alignments
 *        should already be grouped by UMI and strand (if applicable).
 *
 * This function processes a collection of partial alignments, where each alignment
 * contains only one UMI (either 5' or 3'). It clusters these alignments based on
 * their positional proximity, considering a specified wiggle room. The alignments
 * are grouped into clusters, which are stored in the provided `clusters` container.
 * The clustering process depends on whether the UMI is 5' or 3', as the alignments
 * are sorted differently in each case.
 *
 * @param wiggle_room The maximum positional difference allowed between alignments
 *                    to be included in the same cluster.
 * @param create_cluster_id A callable that generates unique cluster IDs.
 * @param partial_alignments A vector of alignment pointers to be clustered. For
 *                           partial reads with 3' UMI only, these alignments will be sorted by
 *                           their end positions; for partial reads with 5' UMI only, they are
 *                           already sorted by their start positions.
 * @param umi5p A boolean indicating whether the UMI is 5' (true) or 3' (false).
 * @param clusters A reference to a container where the resulting clusters will
 *                 be stored. Each cluster is identified by its unique cluster ID.
 *
 * @note This function assumes that the input alignments are already grouped by UMI
 *       and strand (if applicable). If the input `partial_alignments` is empty, the
 *       function returns early without performing any clustering.
 */
static void ClusterUnassignedPartialAlignmentsByPosition(const u32 wiggle_room,
                                                         const CreateClusterId& create_cluster_id,
                                                         vec<AlignmentPtr>& partial_alignments,
                                                         const bool umi5p,
                                                         Clusters& clusters) {
  // The following code assumes there is at least one alignment, return early if there are none.
  if (partial_alignments.empty()) {
    return;
  }

  // The input alignments must be sorted by EndPos() for 3' UMI reads, but are already sorted by StartPos() for 5' UMI
  // reads.
  if (!umi5p) {
    std::ranges::sort(partial_alignments, [](const auto& a, const auto& b) { return a->EndPos() < b->EndPos(); });
  }

  auto current_cluster = std::make_shared<Cluster>(create_cluster_id(), vec<AlignmentPtr>{});

  // Start with the first alignment.
  auto& first_alignment = partial_alignments.front();
  first_alignment->cluster = current_cluster.get();
  current_cluster->alignments.emplace_back(first_alignment);
  u32 last_cluster_pos = umi5p ? first_alignment->StartPos() : first_alignment->EndPos();

  // Iterate through the remaining alignments and cluster them based on positional proximity.
  auto remaining_alignments = partial_alignments | std::views::drop(1);
  auto non_supplementary_or_secondary_alignments =
      remaining_alignments | std::views::filter([](const AlignmentPtr& alignment) {
        return (alignment->record->core.flag & (BAM_FSUPPLEMENTARY | BAM_FSECONDARY)) == 0;
      });
  for (const auto& alignment : non_supplementary_or_secondary_alignments) {
    const u32 current_pos = umi5p ? alignment->StartPos() : alignment->EndPos();

    // Check positional proximity to the last read added.
    if (current_pos <= last_cluster_pos + wiggle_room) {
      // Merge into current cluster
      alignment->cluster = current_cluster.get();
      current_cluster->alignments.emplace_back(alignment);
      // Update the pivot position to the new read's position.
      last_cluster_pos = current_pos;
    } else {
      // Finalize the current cluster, store it, and start a new one.
      clusters[current_cluster->cluster_id] = current_cluster;

      current_cluster = std::make_shared<Cluster>(create_cluster_id(), vec<AlignmentPtr>{});
      alignment->cluster = current_cluster.get();
      current_cluster->alignments.emplace_back(alignment);
      last_cluster_pos = current_pos;
    }
  }

  // Add the final cluster.
  if (!current_cluster->alignments.empty()) {
    clusters[current_cluster->cluster_id] = current_cluster;
  }
}

/**
 * Positional and Fixed UMI clustering is broken up into six phases:
 *   - Phase 1: Parse UMIs from all alignments and separate full read, partial 5p UMI read, and partial 3p UMI read
 *              alignments
 *   - Phase 2: Perform UMI clustering of all full read alignments
 *   - Phase 3: Perform positional clustering on each UMI cluster to further split clusters by position (considering
 *              wiggle room and strand strategy)
 *   - Phase 4: Calculate the mean/min/max start and end positions of the clusters (used for assigning partial UMI
 *              reads to clusters)
 *   - Phase 5: Assign partial UMI alignments to clusters based on the nearest position and matching UMI.
 *   - Phase 6: Cluster the remaining unassigned partial alignments if specified by options.
 * Note: @p alignments must be sorted by start position.
 */
void ClusterAlignmentsByPositionAndUmi(const ReadCollapserOptions& options,
                                       const CreateClusterId& create_cluster_id,
                                       vec<AlignmentPtr>& alignments,
                                       Clusters& clusters) {
  auto& metrics = ConcurrentMetrics::Get();
  const auto initial_clustering_reads = alignments.size();
  // The following code assumes there is at least one alignment, return early if there are none.
  if (alignments.empty()) {
    return;
  }

  // Parse UMIs.
  ParseUmis(alignments, metrics);
  // Remove reads without valid UMIs from the clustering input read count
  metrics.clustering_input_reads -= (initial_clustering_reads - alignments.size());

  // Assign the alignment to the appropriate vector based on the presence of UMIs.
  vec<AlignmentPtr> full_umi_alignments;
  vec<AlignmentPtr> partial_umi5p_alignments;
  vec<AlignmentPtr> partial_umi3p_alignments;
  for (const auto& alignment : alignments) {
    if (alignment->umi5p && alignment->umi3p) {
      full_umi_alignments.emplace_back(alignment);
    } else if (alignment->umi5p) {
      partial_umi5p_alignments.emplace_back(alignment);
    } else if (alignment->umi3p) {
      partial_umi3p_alignments.emplace_back(alignment);
    }
  }

  // Perform UMI clustering to group full reads together based on 5p and 3p UMIs.
  Clusters umi_clusters;
  ClusterFullAlignmentsByUmi(create_cluster_id, full_umi_alignments, umi_clusters);

  // Further break down UMI clusters by position.
  for (const auto& [cluster_id, cluster_ptr] : umi_clusters) {
    // Perform positional clustering on the alignments in the UMI cluster.
    ClusterAlignmentsByPosition(
        options.wiggle_room, options.cluster_by_strand, create_cluster_id, cluster_ptr->alignments, clusters);
  }

  if (options.exclude_partial_reads) {
    // Note partial reads are still included in the alignments vector, but are not clustered.
    const size_t num_unclustered = partial_umi5p_alignments.size() + partial_umi3p_alignments.size();
    metrics.unclustered_partial_reads += num_unclustered;
    metrics.clustering_input_reads -= num_unclustered;
    return;
  }

  // Calculate the mean start and end positions of the clusters to determine the nearest cluster for partial UMI
  // alignments. The min start and max end positions are used to determine if a partial UMI alignment is too long for
  // a cluster.
  UmiClusterPtrAndStats umi5p_mean_min_max_pos_clusters;
  UmiClusterPtrAndStats umi3p_mean_min_max_pos_clusters;
  CalculateClusterMeanMinMaxPosition(
      clusters, options.cluster_by_strand, umi5p_mean_min_max_pos_clusters, umi3p_mean_min_max_pos_clusters);

  // Assign partial UMI alignments to clusters with matching UMI and nearest mean position.
  // Store unassigned partial UMI alignments for further clustering if specified by options.
  UnassignedPartialAlignments unassigned_partial_umi5p_alignments;
  UnassignedPartialAlignments unassigned_partial_umi3p_alignments;
  AssignPartialAlignmentsToClusters(
      options, partial_umi5p_alignments, umi5p_mean_min_max_pos_clusters, unassigned_partial_umi5p_alignments, true);
  AssignPartialAlignmentsToClusters(
      options, partial_umi3p_alignments, umi3p_mean_min_max_pos_clusters, unassigned_partial_umi3p_alignments, false);

  // If clustering of remaining unassigned partial alignments is not specified by options, leave them unclustered and
  // update the metrics.
  if (!options.make_clusters_of_partial_reads_only) {
    size_t num_unclustered = 0;
    for (const auto& partial_alignments : unassigned_partial_umi5p_alignments | std::views::values) {
      num_unclustered += partial_alignments.size();
    }
    for (const auto& partial_alignments : unassigned_partial_umi3p_alignments | std::views::values) {
      num_unclustered += partial_alignments.size();
    }
    metrics.clustering_unclustered_partial_reads += num_unclustered;
    // Since these reads were still considered for clustering, they must be counted towards clustering_reads and
    // clustering_partial_reads
    metrics.clustering_reads += num_unclustered;
    metrics.clustering_partial_reads += num_unclustered;
    // Note partial reads are still included in the alignments vector, but are not clustered.
    return;
  }

  // Cluster remaining unassigned partial alignments if specified by options.
  // These alignments are already grouped by UMI and strand (if applicable).
  for (auto& [umi_strand, partial_alignments] : unassigned_partial_umi5p_alignments) {
    ClusterUnassignedPartialAlignmentsByPosition(
        options.wiggle_room_partial, create_cluster_id, partial_alignments, true, clusters);
  }
  for (auto& [umi_strand, partial_alignments] : unassigned_partial_umi3p_alignments) {
    ClusterUnassignedPartialAlignmentsByPosition(
        options.wiggle_room_partial, create_cluster_id, partial_alignments, false, clusters);
  }
}

/**
 * Cluster alignments and update clustering metrics.
 *
 * This function modifies the `alignments` vector by assigning each read in the vector
 * a pointer to the cluster it belongs to. Unclustered reads will have a `nullptr` cluster pointer.
 * Any supplementary or secondary alignments are ignored during clustering.
 *
 * Unmapped reads are not clustered either as they should be handled separately by
 * `ReadAndWriteUnmappedAlignments`.
 *
 * @param options Read collapser options.
 * @param alignments Vector of alignments to cluster.
 * @param create_cluster_id Callback function to create a unique cluster ID.
 *
 * @return Clusters containing the clustered alignments.
 */
Clusters ClusterAlignments(const ReadCollapserOptions& options,
                           vec<AlignmentPtr>& alignments,
                           const CreateClusterId& create_cluster_id) {
  auto& metrics = ConcurrentMetrics::Get();
  vec<AlignmentPtr> primary_alignments;
  primary_alignments.reserve(alignments.size());
  for (const auto& alignment : alignments) {
    const auto flag = alignment->record->core.flag;
    if ((flag & BAM_FUNMAP) != 0) {
      alignment->cluster = nullptr;
      ++metrics.unmapped_reads;
    } else if ((flag & BAM_FSUPPLEMENTARY) != 0) {
      alignment->cluster = nullptr;
      ++metrics.unclustered_supplementary_reads;
    } else if ((flag & BAM_FSECONDARY) != 0) {
      alignment->cluster = nullptr;
      ++metrics.unclustered_secondary_reads;
    } else {
      primary_alignments.emplace_back(alignment);
    }
  }

  // Clustering input reads are all primary alignments that were considered for clustering, which includes clustered
  // reads and unclustered partial reads (clustering_unclustered_partial_reads), but excludes supplementary and
  // secondary reads, as well as discarded reads with missing UMIs and any partial reads excluded from the analysis with
  // --exclude-partial-reads (unclustered_partial_reads). UMI parsing is required to determine which reads are discarded
  // due to missing UMIs and which partial reads are unclustered due to the --exclude-partial-reads option, so this
  // metric must be decremented in ClusterAlignmentsByPositionAndUmi in these cases.
  metrics.clustering_input_reads += primary_alignments.size();
  Clusters clusters;

  // Importantly, the alignments vector is modified in place to assign each read a pointer to the cluster it belongs to.
  // While primary_alignments is a filtered view of alignments, the AlignmentPtr objects within it are the same as those
  // in alignments, so modifying their cluster pointers also modifies those in alignments. This ensures that after
  // clustering, the original alignments vector, containing ALL reads (clustered, unclustered, supplementary,
  // secondary), has the correct cluster pointers assigned and is preserved for I/O operations.
  if (options.cluster_by_umi) {
    ClusterAlignmentsByPositionAndUmi(options, create_cluster_id, primary_alignments, clusters);
  } else {
    ClusterAlignmentsByPosition(
        options.wiggle_room, options.cluster_by_strand, create_cluster_id, primary_alignments, clusters);
  }

  // Duplicate reads are based on pre-duplex deconvolution cluster sizes.
  UpdateDuplicateReadMetrics(clusters);

  // Deconvolve duplex reads if the option is enabled
  if (options.duplex_library_type != HDDeconvolutionType::kNone) {
    for (const auto& [cluster_id, cluster] : clusters) {
      const bool is_parent_parent = (options.duplex_library_type == HDDeconvolutionType::kParentParent);
      cluster->alignments = DeconvolveDuplexReads(cluster->alignments, is_parent_parent);
    }
  }

  UpdateClusterMetrics(clusters);
  return clusters;
}

}  // namespace xoos::read_collapser
