#pragma once

#include <functional>
#include <memory>
#include <unordered_map>

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/float.h>

#include "core/read-collapser-options.h"
#include "io/alignment.h"

namespace xoos::read_collapser {

/// Unique and reproducible identifier for a cluster.
struct ClusterId {
  u32 super_region_id;  /// The id of the super region this cluster belongs to.
  u32 count;            /// The count of the cluster within the super region.

  auto operator<=>(const ClusterId& cluster_id) const = default;
};

struct ClusterIdHash {
  size_t operator()(const ClusterId& cluster_id) const;
};

struct UmiPair {
  Umi umi5p;
  Umi umi3p;

  auto operator<=>(const UmiPair&) const = default;
};

struct UmiPairHash {
  size_t operator()(const UmiPair& umi_pair) const;
};

struct Cluster {
  ClusterId cluster_id{};
  vec<AlignmentPtr> alignments{};
};

using ClusterPtr = std::shared_ptr<Cluster>;

/// Produce a unique cluster id
using CreateClusterId = std::function<ClusterId()>;

using Clusters = std::unordered_map<ClusterId, ClusterPtr, ClusterIdHash>;

/**
 * Cluster alignments based on their start and end positions and strand information. This clustering method
 * allows for something we call "drift" which is something we have observed where alignments from the same
 * cluster can have slightly different start and end positions. To address this we allow for a "wiggle room"
 * and further alignments are added to the cluster in recursive manner within the wiggle room threshold.
 * @param [in] wiggle_room The maximum |a1.start - a2.start| + |a1.end - a2.end| for two alignments to be in the same
 * cluster in each recursive step.
 * @param [in] cluster_by_strand If kTogether, alignments from different strands are clustered together. If kSeparate,
 * alignments from different strands are clustered separately.
 * @param [in] create_cluster_id A function that produces a unique cluster id on each invocation.
 * @param [in] alignments The alignments to cluster.
 * @param [out] clusters The clusters produced by clustering the alignments.
 **/
void ClusterAlignmentsByPosition(u32 wiggle_room,
                                 bool cluster_by_strand,
                                 const CreateClusterId& create_cluster_id,
                                 const vec<AlignmentPtr>& alignments,
                                 Clusters& clusters);

/**
 * Cluster alignments based on their position and UMI information. This clustering method
 * allows for something we call "drift" which is something we have observed where alignments from the same
 * cluster can have slightly different start and end positions. To address this we allow for a "wiggle room"
 * and further alignments are added to the cluster in recursive manner within the wiggle room threshold.
 * Full alignments are first clustered based on matching UMI and then position.
 * Partial UMI alignments are assigned to clusters based on the nearest position and matching UMI (if exclude partial
 *reads is not enabled). Any remaining unassigned partial UMI alignments are clustered together (if exclude partial
 *reads is not enabled and make clusters of partial reads only is enabled).
 * @param [in] options ReadCollapserOptions containing the cluster by strand, wiggle room, partial read wiggle room,
 * exclude partial reads, and make clusters of partial reads only parameters.
 * @param [in] create_cluster_id A function that produces a unique cluster id on each invocation.
 * @param [in] alignments The alignments to cluster.
 * @param [out] clusters The clusters produced by clustering the alignments.
 *
 * @note This function modifies the input alignments to remove any missing UMI alignments. Any unclustered partial UMI
 * alignments remain in the input alignments.
 **/
void ClusterAlignmentsByPositionAndUmi(const ReadCollapserOptions& options,
                                       const CreateClusterId& create_cluster_id,
                                       vec<AlignmentPtr>& alignments,
                                       Clusters& clusters);

/// Represents the strand of a cluster. A cluster can contain alignments from the forward strand, reverse strand, or
/// both strands.
enum class ClusterStrand {
  kFwd,
  kRev,
  kBoth,
};

/// Determine the strand of the cluster for this alignment based on the alignment strand and whether clusters are
/// separated by strand or not.
ClusterStrand DetermineClusterStrand(const AlignmentPtr& alignment, bool cluster_by_strand);

// Update the cluster metrics based on the clusters provided (duplicate read metrics are handled separately).
void UpdateClusterMetrics(const Clusters& clusters);

/// Represents a cluster coordinate, which is a unique identifier for a cluster during the positional
/// clustering process. ClusterCoord is used to lookup neighboring partial clusters during the clustering process to
/// merge them into a single cluster if they are within the wiggle room and strand constraints.
struct ClusterCoord {
  u32 start;
  u32 end;
  ClusterStrand strand;

  auto operator<=>(const ClusterCoord& coord) const = default;
};

ClusterCoord CreateClusterCoord(const AlignmentPtr& alignment, bool cluster_by_strand);

struct ClusterCoordHash {
  size_t operator()(const ClusterCoord& coord) const;
};

/// Represents a partial cluster, which is a cluster that has not yet been fully merged with its neighboring clusters.
struct PartialCluster {
  ClusterPtr cluster;  /// The cluster that this partial cluster will ultimately be merged into. Can be nullptr if not
                       /// yet determined.
  vec<AlignmentPtr> alignments;
};

using PartialClusters = std::unordered_map<ClusterCoord, PartialCluster, ClusterCoordHash>;

// Represents a cluster and its statistics, used to assign partial UMI alignments to clusters based on the nearest full
// UMI cluster.
struct ClusterPtrAndStats {
  ClusterPtr cluster{};
  ClusterStrand strand{};
  f64 mean_pos{};
  u32 min_start{};
  u32 max_end{};
};

using UmiClusterPtrAndStats = std::unordered_map<Umi, vec<ClusterPtrAndStats>>;

// Represents a UMI and strand pair, used to group unassigned partial alignments by UMI and strand for further
// clustering.
struct UmiStrand {
  Umi umi{};
  ClusterStrand strand{};

  auto operator<=>(const UmiStrand&) const = default;
};

struct UmiStrandHash {
  size_t operator()(const UmiStrand& umi_strand) const;
};

using UnassignedPartialAlignments = std::unordered_map<UmiStrand, vec<AlignmentPtr>, UmiStrandHash>;

/**
 * Perform a depth-first search to merge neighboring partial clusters into a single cluster. This function starts from
 * @ref initial_coord and recursively searches for neighboring partial clusters from @ref partial_clusters within the
 * wiggle room and with the same strand. The search continues until the maximum depth is reached or no more neighboring
 * partial clusters are found.
 */
void DepthFirstClusterSearch(const ClusterCoord& initial_coord,
                             const PartialCluster& partial_cluster,
                             PartialClusters& partial_clusters,
                             u32 wiggle_room,
                             u32 max_depth);

Clusters ClusterAlignments(const ReadCollapserOptions& options,
                           vec<AlignmentPtr>& alignments,
                           const CreateClusterId& create_cluster_id);

std::tuple<Umi, Umi> ParseUmi(const std::string_view& qname);

}  // namespace xoos::read_collapser
