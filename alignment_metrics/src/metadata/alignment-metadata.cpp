#include "metadata/alignment-metadata.h"

#include <exception>
#include <string>
#include <string_view>

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/log/logging.h>
#include <xoos/read-name/consensus-read-name.h>
#include <xoos/read-name/read-name-parser.h>
#include <xoos/types/int.h>

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
 * Parses the read name using the shared `xoos::read_name` parser to extract cluster information
 * (size, type, strand) from post-consensus data and read type (full/partial) from the SID/UMI
 * presence bitfield in pre-consensus data. Also extracts read strand from BAM flags.
 *
 * Post-consensus names follow the format:
 * {cluster_id}-{fwd_partial}-{rev_partial}-{fwd_full}-{rev_full}-{avg_family_size}.
 *
 * Pre-consensus names follow the shared SBX read-name format. Full/partial status is derived from
 * the adapter bitfield reported by the parser rather than from UMI presence alone: for UMI-prep
 * datasets a read is full when both UMIs are classified, and for non-UMI datasets a read is full
 * when both SIDs are classified. This means partial/full can now be determined for reads without
 * UMIs (via the SID bits in the bitfield).
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

  const std::string_view read_name(bam_get_qname(alignment));
  if (!ignore_family) {
    // Post-consensus (inter-molecular consensus) read names encode the cluster composition. Parse
    // them with the shared parser; a non-consensus read name simply yields no cluster info.
    if (const auto consensus_info = read_name::ParseConsensusReadNameNoExcept(read_name)) {
      cluster_size = consensus_info->avg_family_size;
      cluster_type = GetClusterType(consensus_info->forward_partial,
                                    consensus_info->reverse_partial,
                                    consensus_info->forward_full,
                                    consensus_info->reverse_full);
      cluster_strand = GetClusterStrand(consensus_info->forward_partial,
                                        consensus_info->reverse_partial,
                                        consensus_info->forward_full,
                                        consensus_info->reverse_full);
    }
  }

  // Pre-consensus reads follow the shared SBX read-name format. The read-name format itself tells us
  // whether the dataset carries UMIs (UMI-prep emits UMI fields even on partial reads), so we decide
  // per read whether to interpret full/partial from the UMI bits or the SID bits of the adapter
  // bitfield. Only names the shared parser can decode (legacy names with UMI fields, or newer names
  // that carry an adapter bitfield) yield a full/partial classification; anything else — including
  // consensus names and names with no adapter/UMI information — leaves read_type as kUnknown.
  try {
    const bool expect_umis = read_name::HasUmiFields(read_name);
    const read_name::AdapterInfo adapter_info = read_name::ParseReadName(read_name, expect_umis);
    // A read is full when both ends are complete. When UMIs are expected, completeness is determined
    // by UMI presence on both ends; otherwise it is determined by SID presence on both ends. This
    // mirrors the completeness logic used by the Read Collapser and lets us classify full/partial for
    // reads without UMIs via the SID bits of the adapter bitfield.
    const bool is_5p_complete = expect_umis ? adapter_info.has_5p_umi : adapter_info.has_5p_sid;
    const bool is_3p_complete = expect_umis ? adapter_info.has_3p_umi : adapter_info.has_3p_sid;
    read_type = (is_5p_complete && is_3p_complete) ? ReadType::kFull : ReadType::kPartial;
  } catch (const std::exception& e) {
    // Read name does not carry decodable adapter/UMI information; leave read_type as kUnknown.
    XOOS_LOG_DEBUG_ONCE("Read name {} does not carry decodable adapter/UMI information: {}. Read type left as unknown.",
                        std::string(read_name),
                        e.what());
  }

  // If the read is a duplex read with YC tag, then it should be treated as a full read
  // regardless of the UMIs/SIDs, since reads with missing SID or hairpin adapter are discarded
  // upstream during demux.
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
