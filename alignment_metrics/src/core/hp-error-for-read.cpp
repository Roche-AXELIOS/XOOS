#include "hp-error-for-read.h"

#include <xoos/log/logging.h>

#include "alignment-metrics-options.h"

namespace xoos::alignment_metrics {

void HomopolymerErrorProfileForRead::AddSubstitution(u8 position) {
  if (position >= kHpMaxLength) {
    return;
  }
  _substitution_by_position |= (1u << position);
}

void HomopolymerErrorProfileForRead::AddInsertion(u8 position, const std::string& sequence) {
  if (position >= kHpMaxLength) {
    return;
  }
  _insertion_by_position |= (1u << position);
  _indel_by_position |= (1u << position);
  _insertion_sequences.emplace_back(sequence);
}

void HomopolymerErrorProfileForRead::AddDeletion(u8 position) {
  if (position >= kHpMaxLength) {
    return;
  }
  _deletion_by_position |= (1u << position);
  _indel_by_position |= (1u << position);
}

void HomopolymerErrorProfileForRead::UpdateMinQuality(u8 quality) {
  if (quality < _min_quality || _min_quality == kHomopolymerUninitializedBaseQuality) {
    _min_quality = quality;
  }
}

void HomopolymerErrorProfileForRead::UpdateMinBaseType(yc_decode::BaseType base_type) {
  if (base_type < _min_base_type) {
    _min_base_type = base_type;
  }
}

u8 HomopolymerErrorProfileForRead::GetMinQuality() const {
  return _min_quality;
}

yc_decode::BaseType HomopolymerErrorProfileForRead::GetMinBaseType() const {
  return _min_base_type;
}

HpErrorByPosition HomopolymerErrorProfileForRead::GetInsertionByPosition() const {
  return _insertion_by_position;
}

HpErrorByPosition HomopolymerErrorProfileForRead::GetDeletionByPosition() const {
  return _deletion_by_position;
}

HpErrorByPosition HomopolymerErrorProfileForRead::GetIndelByPosition() const {
  return _indel_by_position;
}

HpErrorByPosition HomopolymerErrorProfileForRead::GetSubstitutionByPosition() const {
  return _substitution_by_position;
}

bool HomopolymerErrorProfileForRead::AreInsertionsHomogeneous(char base) const {
  return std::ranges::all_of(_insertion_sequences, [base](const std::string_view seq) {
    return !seq.empty() && seq.find_first_not_of(base) == std::string::npos;
  });
}

};  // namespace xoos::alignment_metrics
