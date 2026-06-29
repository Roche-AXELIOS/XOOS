#pragma once

#include <xoos/types/float.h>
#include <xoos/yc-decode/yc-decoder.h>

#include "io/alignment.h"

namespace xoos::read_collapser {

struct ClusterMetadata {
  u32 forward_partial = 0;
  u32 reverse_partial = 0;
  u32 forward_full = 0;
  u32 reverse_full = 0;
  u32 forward_total = 0;
  u32 reverse_total = 0;
};

struct HomopolymerRange {
  u64 start;
  u64 end;  // inclusive
};

using HomopolymerRanges = vec<HomopolymerRange>;

// Class to hold the consensus matrix, the strand information for each read
// in the cluster, and a metadata object that keeps track of the number of
// different types of reads (forward, reverse, partial, full) in the cluster.
class ConsensusMatrix {
 public:
  /**
   * Construct a ConsensusMatrix with the given read count and consensus length.
   * The matrix is initialized with gaps ('-') and the strand vector is initialized
   * with ReadStrand::kBoth.
   */
  explicit ConsensusMatrix(size_t read_count, size_t consensus_length);

  ConsensusMatrix() = default;

  // Set the strand for the `read_index` read in the cluster.
  void SetStrand(size_t read_index, ReadStrand new_strand);

  // Get the strand for the `read_index` read in the cluster.
  ReadStrand GetStrand(size_t read_index) const;

  // Set the duplex strand for the `read_index` read in the cluster.
  void SetDuplexStrand(size_t read_index, DuplexStrand new_duplex_strand);

  // Get the duplex strand for the `read_index` read in the cluster.
  DuplexStrand GetDuplexStrand(size_t read_index) const;

  // Check if the read at `read_index` is a deconvolved duplex read.
  bool IsDuplex(size_t read_index) const;

  // Set the base for the `read_index` read at the given position in the consensus matrix.
  void SetBase(size_t read_index, size_t pos, char base);

  // Fill the entire row for the `read_index` read with the given base.
  void FillRow(size_t read_index, char base);

  // Get the base for the `read_index` read at the given position in the consensus matrix.
  char GetBase(size_t read_index, size_t pos) const;

  // Set the base type for the `read_index` read at the given position in the consensus matrix.
  void SetBaseType(size_t read_index, size_t pos, yc_decode::BaseType base_type);

  // Get the base type for the `read_index` read at the given position in the consensus matrix.
  yc_decode::BaseType GetBaseType(size_t read_index, size_t pos) const;

  // Set the base types for all positions in the consensus matrix.
  void SetBaseTypes(size_t read_count, size_t consensus_length, yc_decode::BaseType base_type);

  // Check if the consensus matrix has base type information.
  bool HasBaseTypeInfo() const;

  /**
   * Update the counts of read types (forward-full, forward-partial, reverse-full, reverse-partial)
   * in the cluster metadata based on the strand and partiality of a read by 1.
   */
  void UpdateMetadata(bool is_reverse, bool is_partial);

  // Get the metadata for the cluster.
  ClusterMetadata GetMetadata() const;

  // Get the number of homopolymer ranges in the consensus matrix.
  size_t GetNumberOfHomopolymerRanges() const;

  // Get the sequence for the `read_index` read in the cluster as a string.
  std::string GetSequence(size_t read_index) const;

  // Get the length of the consensus matrix.
  size_t GetConsensusLength() const;

  // Get the number of reads in the cluster.
  size_t GetReadCount() const;

  /**
   * @brief Adds a new homopolymer range to the consensus matrix, merging with existing overlapping or adjacent ranges.
   *
   * This method inserts a new interval [start, end] into the `homopolymer_ranges` container. If the new interval
   * overlaps with to any existing intervals, those intervals are merged into a single interval that covers the union of
   * all overlapping ranges. If the new interval is completely contained within an existing interval, no changes are
   * made. The intervals in `homopolymer_ranges` are maintained in sorted order and are non-overlapping.
   *
   * @param start The inclusive start position of the homopolymer range.
   * @param end The inclusive end position of the homopolymer range. Must satisfy start < end.
   *
   * @throws error::Error if start >= end.
   */
  void AddHomopolymerRange(u64 start, u64 end);

  vec<HomopolymerRange> GetHomopolymerRanges() const {
    return _homopolymer_ranges;
  }

  // Check if a position is part of a homopolymer range.
  bool IsPartOfHomopolymer(u64 pos) const;

 private:
  vec<vec<char>> _msa_matrix{};
  vec<vec<yc_decode::BaseType>> _base_types{};
  vec<ReadStrand> _strand{};
  vec<DuplexStrand> _duplex_strand{};
  // Column ranges where homopolymers are present in the consensus matrix.
  HomopolymerRanges _homopolymer_ranges{};

  ClusterMetadata _metadata{};
};

}  // namespace xoos::read_collapser
