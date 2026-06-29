#pragma once

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/int.h>

namespace xoos::alignment_metrics {

// Strand composition of a read cluster (forward-only, reverse-only, mixed, or unknown).
enum class ClusterStrand : u8 {
  kMixedStrandCluster,
  kForwardCluster,
  kReverseCluster,
  kUnknown
};

// Completeness of reads in a cluster (full-length, partial, mixed, or unknown).
enum class ClusterType : u8 {
  kFull,
  kPartial,
  kMixed,
  kUnknown
};

// Strand orientation of a single read (forward, reverse, or unknown).
enum class ReadStrand : u8 {
  kForward,
  kReverse,
  kUnknown
};

// Read length type indicating whether a read is full-length, partial, or unknown.
enum class ReadType {
  kFull,
  kPartial,
  kUnknown
};

/**
 * @brief Metadata extracted from alignment records describing cluster and read characteristics.
 *
 * Contains information about cluster composition (strand, type, size) and individual read
 * properties (strand, type) parsed from BAM record names and flags. Used to stratify metrics
 * by cluster and read characteristics in post-consensus datasets.
 */
struct AlignmentMetadata {
  ClusterStrand cluster_strand{ClusterStrand::kUnknown};
  ClusterType cluster_type{ClusterType::kUnknown};
  u32 cluster_size{1};
  ReadStrand read_strand{ReadStrand::kUnknown};
  ReadType read_type{ReadType::kUnknown};
};

ClusterType GetClusterType(u32 forward_partial_count,
                           u32 reverse_partial_count,
                           u32 forward_full_count,
                           u32 reverse_full_count);

ClusterStrand GetClusterStrand(u32 forward_partial_count,
                               u32 reverse_partial_count,
                               u32 forward_full_count,
                               u32 reverse_full_count);

/**
 * Create alignment metadata object from the read name of a consensus read
 * that include information on the cluster from which the consensus read was generated.
 *
 * The read name of a consensus read is expected to be in the format:
 * `{cluster_id}-{forward_partial_count}-{reverse_partial_count}-{forward_full_count}-{reverse_full_count}-{average_family_size}`
 * where `cluster_id` is of the format `{N}:{M}` so that `N` is the batch ID and `M` is the cluster ID within the batch.
 *
 * @param alignment
 * @param ignore_family
 * @return AlignmentMetadata object
 */
AlignmentMetadata CreateAlignmentMetadata(const bam1_t* alignment, bool ignore_family);

}  // namespace xoos::alignment_metrics
