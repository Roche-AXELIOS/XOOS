#pragma once

#include <functional>
#include <memory>
#include <unordered_map>

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/read-name/umi.h>
#include <xoos/types/float.h>

#include "core/read-collapser-options.h"
#include "io/alignment.h"
#include "metrics/metrics.h"

namespace xoos::read_collapser {

/// Unique and reproducible identifier for a cluster.
struct ClusterId {
  /// The id of the super region this cluster belongs to.
  u32 super_region_id;
  /// The count of the cluster within the super region.
  u32 count;

  auto operator<=>(const ClusterId& cluster_id) const = default;
};

struct ClusterIdHash {
  size_t operator()(const ClusterId& cluster_id) const;
};

struct Cluster {
  ClusterId cluster_id{};
  vec<AlignmentPtr> alignments;
};

using ClusterPtr = std::shared_ptr<Cluster>;

/// Produce a unique cluster id
using CreateClusterId = std::function<ClusterId()>;

using Clusters = std::unordered_map<ClusterId, ClusterPtr, ClusterIdHash>;

/**
 * Cluster reads based on their alignment start and end positions and strand information.
 * If cluster_by_umi is enabled, the reads are first clustered by UMI and then by position.
 * This clustering method allows for something we call "drift" which is something we have observed where reads from the
 * same cluster can have slightly different start and end positions. To address this we allow for a "wiggle room" and
 * further reads are added to the cluster in recursive manner within the wiggle room threshold.
 *
 * Full reads are first clustered by position using wiggle_room. Partial reads (IsFivePrimeComplete/
 * IsThreePrimeComplete) are then assigned to the nearest full-read cluster within
 * wiggle_room_partial. Remaining unassigned partial reads are optionally clustered into partial-only
 * clusters if make_clusters_of_partial_reads_only is enabled.
 *
 * @param [in] options ReadCollapserOptions containing wiggle_room, wiggle_room_partial, cluster_by_strand,
 * exclude_partial_reads, and make_clusters_of_partial_reads_only parameters.
 * @param [in] create_cluster_id A function that produces a unique cluster id on each invocation.
 * @param [in] reads The reads to cluster. Note that this is a const reference to a vector of AlignmentPtrs,
 * but the alignments themselves are modified by assigning cluster information to them.
 * @param [out] clusters The clusters produced by clustering the reads.
 * @param [in,out] metrics The metrics to accumulate clustering metrics into.
 **/
void ClusterReads(const ReadCollapserOptions& options,
                  const CreateClusterId& create_cluster_id,
                  const vec<AlignmentPtr>& reads,
                  Clusters& clusters,
                  Metrics& metrics);

/// Represents the strand of a cluster. A cluster can contain reads aligned to the forward strand, reverse strand, or
/// both strands.
enum class ClusterStrand : u8 {
  kFwd,
  kRev,
  kBoth,
};

/// Determine the strand of the cluster for this read based on the read strand and whether clusters are
/// separated by strand or not.
ClusterStrand DetermineClusterStrand(const AlignmentPtr& alignment, bool cluster_by_strand);

// Update the cluster metrics based on the clusters provided (duplicate read metrics are handled separately).
void UpdateClusterMetrics(const Clusters& clusters, Metrics& metrics);

/// Represents a cluster coordinate, which is a unique identifier for a cluster during the positional
/// clustering process. ClusterCoord is used to lookup neighboring partial clusters during the clustering process to
/// merge them into a single cluster if they are within the wiggle room and strand constraints.
struct ClusterCoord {
  // 0-based inclusive start position
  u32 start;
  // 0-based exclusive end position
  u32 end;
  ClusterStrand strand;

  auto operator<=>(const ClusterCoord& coord) const = default;
};

ClusterCoord CreateClusterCoord(const AlignmentPtr& alignment, bool cluster_by_strand);

struct ClusterCoordHash {
  size_t operator()(const ClusterCoord& coord) const;
};

/// Represents a preliminary cluster, which is a cluster that has not yet been fully merged with its neighboring
/// clusters.
struct PreliminaryCluster {
  /**
   * The cluster that this preliminary cluster will ultimately be merged into. Can be nullptr if not yet determined.
   */
  ClusterPtr cluster;
  vec<AlignmentPtr> alignments;
};

using PreliminaryClusters = std::unordered_map<ClusterCoord, PreliminaryCluster, ClusterCoordHash>;

// Represents a cluster and its statistics, used to assign partial reads to clusters based on the nearest
// full read cluster.
struct ClusterPtrAndStats {
  ClusterPtr cluster;
  ClusterStrand strand{};
  f64 mean_pos{};
  u32 min_start{};
  u32 max_end{};
};

using UmiClusterPtrAndStats = std::unordered_map<read_name::Umi, vec<ClusterPtrAndStats>>;

// Represents a UMI and strand pair, used to group unassigned partial reads by UMI and strand for further
// clustering.
struct UmiAndStrand {
  read_name::Umi umi;
  ClusterStrand strand{};

  auto operator<=>(const UmiAndStrand&) const = default;
};

struct UmiAndStrandHash {
  size_t operator()(const UmiAndStrand& umi_strand) const;
};

using UnassignedPartialReads = std::unordered_map<UmiAndStrand, vec<AlignmentPtr>, UmiAndStrandHash>;

/**
 * Perform a depth-first search to merge neighboring preliminary clusters into a single cluster. This function starts
 * from
 * @ref initial_coord and recursively searches for neighboring preliminary clusters from @ref preliminary_clusters
 * within the wiggle room and with the same strand. The search continues until the maximum depth is reached or no more
 * neighboring preliminary clusters are found.
 */
void MergePreliminaryClusters(const ClusterCoord& initial_coord,
                              const PreliminaryCluster& preliminary_cluster,
                              PreliminaryClusters& preliminary_clusters,
                              u32 wiggle_room,
                              u32 max_depth);

/**
 * @brief Cluster alignments based on UMI (if enabled) and alignment position proximity.
 * This function takes a range of position-sorted alignments fetched from a given region
 * from the BAM file and clusters them based on their position and optionally their UMIs
 * and strand. Only primary alignments are clustered as they represent distinct reads.
 * Secondary and supplementary alignments are ignored for clustering purposes as they
 * represent the same read as the primary alignment.
 *
 * @param options ReadCollapserOptions containing clustering parameters such as wiggle_room, cluster_by_strand, and
 *                cluster_by_umi.
 * @param alignments A vector of AlignmentPtrs representing the alignments to be clustered.
 *                   These alignments are expected to be sorted by their alignment start position.
 *                   The alignments are modified in-place to assign cluster information to them as they are clustered.
 * @param create_cluster_id A function that produces a unique cluster id on each invocation.
 * @param metrics The metrics to accumulate clustering metrics into.
 * @return A map of ClusterId to ClusterPtr representing the resulting clusters after clustering the input alignments.
 *
 */
Clusters ClusterAlignments(const ReadCollapserOptions& options,
                           vec<AlignmentPtr>& alignments,
                           const CreateClusterId& create_cluster_id,
                           Metrics& metrics);

}  // namespace xoos::read_collapser
