#include "consensus/consensus-matrix.h"

#include <algorithm>

#include <xoos/error/error.h>
#include <xoos/types/int.h>
#include <xoos/yc-decode/yc-decoder.h>

#include "consensus/base-encoder.h"
#include "io/alignment.h"

namespace xoos::read_collapser {

ConsensusMatrix::ConsensusMatrix(size_t read_count, size_t consensus_length)
    : _msa_matrix(read_count, vec<char>(consensus_length, kBaseGap)),
      _strand(read_count, ReadStrand::kFwd),
      _duplex_strand(read_count, DuplexStrand::kSimplex) {
}

void ConsensusMatrix::SetStrand(size_t read_index, ReadStrand new_strand) {
  if (read_index >= _strand.size()) {
    throw error::Error("Invalid read index {} for setting strand", read_index);
  }
  _strand[read_index] = new_strand;
}

ReadStrand ConsensusMatrix::GetStrand(size_t read_index) const {
  if (read_index >= _strand.size()) {
    throw error::Error("Invalid read index {} for getting strand", read_index);
  }
  return _strand[read_index];
}

void ConsensusMatrix::SetDuplexStrand(size_t read_index, DuplexStrand new_duplex_strand) {
  if (read_index >= _duplex_strand.size()) {
    throw error::Error("Invalid read index {} for setting duplex strand", read_index);
  }
  _duplex_strand[read_index] = new_duplex_strand;
}

DuplexStrand ConsensusMatrix::GetDuplexStrand(size_t read_index) const {
  if (read_index >= _duplex_strand.size()) {
    throw error::Error("Invalid read index {} for getting duplex strand", read_index);
  }
  return _duplex_strand[read_index];
}

bool ConsensusMatrix::IsDuplex(size_t read_index) const {
  if (read_index >= _duplex_strand.size()) {
    throw error::Error("Invalid read index {} for checking duplex", read_index);
  }
  return _duplex_strand[read_index] == DuplexStrand::kR1 || _duplex_strand[read_index] == DuplexStrand::kR2;
}

void ConsensusMatrix::SetBase(size_t read_index, size_t pos, char base) {
  if (read_index >= _msa_matrix.size() || pos >= _msa_matrix[read_index].size()) {
    throw error::Error("Invalid read index {} or position {} for setting base", read_index, pos);
  }
  _msa_matrix[read_index][pos] = base;
}

void ConsensusMatrix::FillRow(const size_t read_index, const char base) {
  if (read_index >= _msa_matrix.size()) {
    throw error::Error("Invalid read index {} for filling row", read_index);
  }
  std::ranges::fill(_msa_matrix[read_index], base);
}

char ConsensusMatrix::GetBase(size_t read_index, size_t pos) const {
  if (read_index >= _msa_matrix.size() || pos >= _msa_matrix[read_index].size()) {
    throw error::Error("Invalid read index {} or position {} for getting base", read_index, pos);
  }
  return _msa_matrix[read_index][pos];
}

void ConsensusMatrix::SetBaseType(size_t read_index, size_t pos, yc_decode::BaseType base_type) {
  if (read_index >= _base_types.size() || pos >= _base_types[read_index].size()) {
    throw error::Error("Invalid read index {} or position {} for setting base type", read_index, pos);
  }
  _base_types[read_index][pos] = base_type;
}

yc_decode::BaseType ConsensusMatrix::GetBaseType(size_t read_index, size_t pos) const {
  if (_base_types.empty()) {
    return yc_decode::BaseType::kSimplex;  // Default base type if no base types are set
  } else if (read_index >= _base_types.size() || pos >= _base_types[read_index].size()) {
    throw error::Error("Invalid read index {} or position {} for getting base type", read_index, pos);
  }
  return _base_types[read_index][pos];
}

void ConsensusMatrix::SetBaseTypes(const size_t read_count,
                                   const size_t consensus_length,
                                   const yc_decode::BaseType base_type) {
  _base_types = vec<vec<yc_decode::BaseType>>(read_count, vec<yc_decode::BaseType>(consensus_length, base_type));
}

bool ConsensusMatrix::HasBaseTypeInfo() const {
  return !_base_types.empty();
}

void ConsensusMatrix::UpdateMetadata(const bool is_reverse, const bool is_partial) {
  if (is_reverse) {
    if (is_partial) {
      ++_metadata.reverse_partial;
    } else {
      ++_metadata.reverse_full;
    }
    ++_metadata.reverse_total;
  } else {
    if (is_partial) {
      ++_metadata.forward_partial;
    } else {
      ++_metadata.forward_full;
    }
    ++_metadata.forward_total;
  }
}

ClusterMetadata ConsensusMatrix::GetMetadata() const {
  return _metadata;
}

size_t ConsensusMatrix::GetNumberOfHomopolymerRanges() const {
  return _homopolymer_ranges.size();
}

std::string ConsensusMatrix::GetSequence(size_t read_index) const {
  if (read_index >= _msa_matrix.size()) {
    throw error::Error("Invalid read index {} for getting sequence", read_index);
  }
  return {_msa_matrix[read_index].begin(), _msa_matrix[read_index].end()};
}

size_t ConsensusMatrix::GetConsensusLength() const {
  return _msa_matrix.empty() ? 0 : _msa_matrix[0].size();
}

size_t ConsensusMatrix::GetReadCount() const {
  return _msa_matrix.size();
}

void ConsensusMatrix::AddHomopolymerRange(u64 start, u64 end) {
  if (start >= end) {
    throw error::Error("Invalid homopolymer range {}-{}", start, end);
  }
  u64 current_start = start;
  u64 current_end = end;

  // Find the range of intervals to merge/remove
  std::optional<u32> merge_start_idx = std::nullopt;
  std::optional<u32> merge_end_idx = std::nullopt;  // Exclusive end index
  for (size_t i = 0; i < _homopolymer_ranges.size(); ++i) {
    const auto& existing_interval = _homopolymer_ranges[i];

    // No overlap, existing interval is completely before new/current merged interval
    if (existing_interval.end < current_start) {
      continue;
    }

    // No overlap, existing interval is completely after new/current merged interval
    if (existing_interval.start > current_end) {
      break;  // Since sorted, no further overlaps
    }

    // Overlap or complete containment
    if (existing_interval.start <= current_start && current_end <= existing_interval.end) {
      return;  // New interval is completely contained, do nothing.
    }

    // Merge
    current_start = std::min(current_start, existing_interval.start);
    current_end = std::max(current_end, existing_interval.end);

    if (!merge_start_idx.has_value()) {
      merge_start_idx = static_cast<u32>(i);
    }
    merge_end_idx = static_cast<u32>(i) + 1;
  }
  // Remove merged intervals
  if (merge_start_idx.has_value()) {
    _homopolymer_ranges.erase(_homopolymer_ranges.begin() + merge_start_idx.value(),
                              _homopolymer_ranges.begin() + merge_end_idx.value());
  }
  // Insert the new merged interval at the correct position
  // Find insertion point using binary search (std::lower_bound)
  const auto it_insert =
      std::ranges::lower_bound(_homopolymer_ranges,
                               HomopolymerRange{current_start, current_end},
                               [](const HomopolymerRange& a, const HomopolymerRange& b) { return a.start < b.start; });
  _homopolymer_ranges.insert(it_insert, HomopolymerRange{current_start, current_end});
}

bool ConsensusMatrix::IsPartOfHomopolymer(u64 pos) const {
  return std::ranges::any_of(_homopolymer_ranges,
                             [pos](const auto& range) { return pos >= range.start && pos <= range.end; });
}

}  // namespace xoos::read_collapser
