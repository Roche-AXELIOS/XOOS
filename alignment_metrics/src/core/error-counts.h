#pragma once

#include "metadata/alignment-metadata.h"

namespace xoos::alignment_metrics {

/**
 * Counts of an event (e.g. error, aligned base, etc.) categorized by read and/or cluster properties
 */
template <typename T>
  requires std::is_integral_v<T> && std::is_unsigned_v<T>
struct CountsByReadProperty {
  T mixed_strand_cluster{};
  T forward_cluster{};
  T reverse_cluster{};
  T full_cluster{};
  T partial_cluster{};
  T mixed_full_and_partial_cluster{};
  T forward_read{};
  T reverse_read{};
  T partial_read{};
  T full_read{};
  T total{};

  void Add(const AlignmentMetadata& alignment_metadata) {
    switch (alignment_metadata.cluster_type) {
      case ClusterType::kFull:
        ++full_cluster;
        break;
      case ClusterType::kPartial:
        ++partial_cluster;
        break;
      case ClusterType::kMixed:
        ++mixed_full_and_partial_cluster;
        break;
      default:
        break;
    }
    switch (alignment_metadata.cluster_strand) {
      case ClusterStrand::kMixedStrandCluster:
        ++mixed_strand_cluster;
        break;
      case ClusterStrand::kForwardCluster:
        ++forward_cluster;
        break;
      case ClusterStrand::kReverseCluster:
        ++reverse_cluster;
        break;
      default:
        break;
    }
    switch (alignment_metadata.read_strand) {
      case ReadStrand::kForward:
        ++forward_read;
        break;
      case ReadStrand::kReverse:
        ++reverse_read;
        break;
      default:
        break;
    }
    switch (alignment_metadata.read_type) {
      case ReadType::kFull:
        ++full_read;
        break;
      case ReadType::kPartial:
        ++partial_read;
        break;
      default:
        break;
    }
    ++total;
  }

  CountsByReadProperty& operator+=(const CountsByReadProperty& other) {
    mixed_strand_cluster += other.mixed_strand_cluster;
    forward_cluster += other.forward_cluster;
    reverse_cluster += other.reverse_cluster;
    full_cluster += other.full_cluster;
    partial_cluster += other.partial_cluster;
    mixed_full_and_partial_cluster += other.mixed_full_and_partial_cluster;
    forward_read += other.forward_read;
    reverse_read += other.reverse_read;
    partial_read += other.partial_read;
    full_read += other.full_read;
    total += other.total;
    return *this;
  }
};

/**
 * Counts of different events (e.g. aligned bases, insertions, deletions, substitutions)
 * categorized by the type of event
 */
template <typename T>
  requires std::is_integral_v<T> && std::is_unsigned_v<T>
struct CountsByErrorType {
  T depth{};
  T insertions{};
  T deletions{};
  T substitutions{};

  CountsByErrorType& operator+=(const CountsByErrorType& other) {
    depth += other.depth;
    insertions += other.insertions;
    deletions += other.deletions;
    substitutions += other.substitutions;
    return *this;
  }
};

/**
 * Counts of different events (e.g. aligned bases, insertions, deletions, substitutions)
 * categorized by the type of event and read/cluster properties
 */
template <typename T>
  requires std::is_integral_v<T> && std::is_unsigned_v<T>
struct CountsByErrorTypeAndReadProperty {
  CountsByReadProperty<T> substitutions;
  CountsByReadProperty<T> insertions;
  CountsByReadProperty<T> deletions;
  CountsByReadProperty<T> totals;
};

/**
 * Count of mismatches categorized by read and/or cluster properties
 */
template <typename T>
struct MismatchCounts {
  T total_bases{};
  T total_mismatches{};
  T forward_cluster_bases{};
  T forward_cluster_mismatches{};
  T reverse_cluster_bases{};
  T reverse_cluster_mismatches{};
  T mixed_strand_cluster_bases{};
  T mixed_strand_cluster_mismatches{};

  void Add(const AlignmentMetadata& alignment_metadata, bool is_mismatch) {
    if (is_mismatch) {
      switch (alignment_metadata.cluster_strand) {
        case ClusterStrand::kMixedStrandCluster:
          ++mixed_strand_cluster_mismatches;
          break;
        case ClusterStrand::kForwardCluster:
          ++forward_cluster_mismatches;
          break;
        case ClusterStrand::kReverseCluster:
          ++reverse_cluster_mismatches;
          break;
        default:
          break;
      }
      ++total_mismatches;
    }
    switch (alignment_metadata.cluster_strand) {
      case ClusterStrand::kMixedStrandCluster:
        ++mixed_strand_cluster_bases;
        break;
      case ClusterStrand::kForwardCluster:
        ++forward_cluster_bases;
        break;
      case ClusterStrand::kReverseCluster:
        ++reverse_cluster_bases;
        break;
      default:
        break;
    }
    ++total_bases;
  }
};

}  // namespace xoos::alignment_metrics
