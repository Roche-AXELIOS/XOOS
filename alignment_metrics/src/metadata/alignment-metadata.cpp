#include "metadata/alignment-metadata.h"

#include <exception>
#include <string>

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>
#include <xoos/util/string-functions.h>

namespace xoos::alignment_metrics {

/**
 * @brief Determines the cluster type based on the counts of full and partial reads.
 *
 * Classifies a read cluster as Full (only full-length reads), Partial (only partial reads),
 * Mixed (both full and partial reads), or Unknown (no reads) based on the provided counts.
 *
 * @param forward_partial_count Count of partial reads on the forward strand
 * @param reverse_partial_count Count of partial reads on the reverse strand
 * @param forward_full_count Count of full-length reads on the forward strand
 * @param reverse_full_count Count of full-length reads on the reverse strand
 * @return ClusterType indicating whether the cluster contains full, partial, mixed, or unknown reads
 */
ClusterType GetClusterType(const u32 forward_partial_count,
                           const u32 reverse_partial_count,
                           const u32 forward_full_count,
                           const u32 reverse_full_count) {
  auto cluster_type = ClusterType::kUnknown;
  if (forward_full_count + reverse_full_count > 0 &&
      (forward_partial_count + reverse_partial_count == 0)) {  // NOLINT - parentheses needed for clarity
    cluster_type = ClusterType::kFull;
  }
  if (forward_full_count + reverse_full_count == 0 &&
      (forward_partial_count + reverse_partial_count > 0)) {  // NOLINT - parentheses needed for clarity
    cluster_type = ClusterType::kPartial;
  }
  if (forward_full_count + reverse_full_count > 0 &&
      (forward_partial_count + reverse_partial_count > 0)) {  // NOLINT - parentheses needed for clarity
    cluster_type = ClusterType::kMixed;
  }
  return cluster_type;
}

/**
 * @brief Determines the strand composition of a read cluster.
 *
 * Classifies a cluster as ForwardCluster (only forward reads), ReverseCluster (only reverse reads),
 * MixedStrandCluster (both forward and reverse reads), or Unknown (no reads).
 *
 * @param forward_partial_count Count of partial reads on the forward strand
 * @param reverse_partial_count Count of partial reads on the reverse strand
 * @param forward_full_count Count of full-length reads on the forward strand
 * @param reverse_full_count Count of full-length reads on the reverse strand
 * @return ClusterStrand indicating the strand composition of the cluster
 */
ClusterStrand GetClusterStrand(const u32 forward_partial_count,
                               const u32 reverse_partial_count,
                               const u32 forward_full_count,
                               const u32 reverse_full_count) {
  auto cluster_strand = ClusterStrand::kUnknown;
  if (forward_full_count + forward_partial_count > 0 &&
      reverse_full_count + reverse_partial_count == 0) {  // NOLINT - parentheses needed for clarity
    cluster_strand = ClusterStrand::kForwardCluster;
  }
  if (forward_full_count + forward_partial_count == 0 &&
      (reverse_full_count + reverse_partial_count > 0)) {  // NOLINT - parentheses needed for clarity
    cluster_strand = ClusterStrand::kReverseCluster;
  }
  if (forward_full_count + forward_partial_count > 0 &&
      (reverse_full_count + reverse_partial_count > 0)) {  // NOLINT - parentheses needed for clarity
    cluster_strand = ClusterStrand::kMixedStrandCluster;
  }
  return cluster_strand;
}

/**
 * @brief Extracts alignment metadata from a BAM record's read name and flags.
 *
 * Parses the read name to extract cluster information (size, type, strand) from post-consensus data,
 * and determines read type (full/partial) from UMI presence in pre-consensus data. Also extracts
 * read strand from BAM flags. Post-consensus names follow the format:
 * {cluster_id}-{fwd_partial}-{rev_partial}-{fwd_full}-{rev_full}-{avg_family_size}.
 * Pre-consensus names with UMIs follow: @{id}|{sid_5p}|{sid_3p}|{umi_5p}|{umi_3p}.
 *
 * @param alignment Pointer to the BAM alignment record
 * @param ignore_family If true, skip parsing cluster/family information from the read name
 * @return AlignmentMetadata containing parsed cluster and read information
 */
AlignmentMetadata CreateAlignmentMetadata(const bam1_t* alignment, const bool ignore_family) {
  u32 cluster_size = 0;
  auto cluster_strand = ClusterStrand::kUnknown;
  auto cluster_type = ClusterType::kUnknown;
  auto read_type = ReadType::kUnknown;
  auto read_strand = ReadStrand::kForward;

  std::string read_name(bam_get_qname(alignment));
  if (!ignore_family) {
    /**
     * Post-consensus read names are formatted as
     * {cluster_id}-{forward_partial_count}-{reverse_partial_count}-{forward_full_count}-{reverse_full_count}-{average_family_size}
     */
    const auto read_name_parts = string::Split(read_name, "-");
    if (read_name_parts.size() == 6) {
      // Parse the read name for more information
      try {
        const s32 forward_partial_count = std::stoi(read_name_parts[1]);
        const s32 reverse_partial_count = std::stoi(read_name_parts[2]);
        const s32 forward_full_count = std::stoi(read_name_parts[3]);
        const s32 reverse_full_count = std::stoi(read_name_parts[4]);
        cluster_size = std::stoi(read_name_parts[5]);
        cluster_type =
            GetClusterType(forward_partial_count, reverse_partial_count, forward_full_count, reverse_full_count);
        cluster_strand =
            GetClusterStrand(forward_partial_count, reverse_partial_count, forward_full_count, reverse_full_count);
      } catch (const std::exception& e) {
        XOOS_LOG_DEBUG_ONCE(
            "Error parsing read name {}: {}. Assuming BAM is not post-consensus data.", read_name, e.what());
      }
    }
  }
  /**
   * Pre-consensus read names with UMIs are formatted as
   * @<original id>|<sid_5p>|<sid_3p>|<umi_5p>|<umi_3p>
   *
   * We currently check if both the 5p and 3p UMIs are present, and if so, we assume the read is a full read;
   * otherwise it is partial.
   *
   * TODO: This may change in the future (e.g. we may indicate full/partial status directly in the read name or
   * tag)
   */
  const auto read_name_parts = string::Split(read_name, "|");
  if (read_name_parts.size() >= 3) {
    if (read_name_parts[read_name_parts.size() - 1] == "*" || read_name_parts[read_name_parts.size() - 2] == "*") {
      read_type = ReadType::kPartial;
    } else {
      read_type = ReadType::kFull;
    }
  }
  // If the read is a duplex read with YC tag, then it should be treated as a full read
  // regardless of the UMIs
  if (bam_aux_get(alignment, "YC") != nullptr) {
    read_type = ReadType::kFull;
  }

  if ((alignment->core.flag & BAM_FREVERSE) != 0) {
    read_strand = ReadStrand::kReverse;
  }

  return AlignmentMetadata{.cluster_strand = cluster_strand,
                           .cluster_type = cluster_type,
                           .cluster_size = cluster_size,
                           .read_strand = read_strand,
                           .read_type = read_type};
}

}  // namespace xoos::alignment_metrics
