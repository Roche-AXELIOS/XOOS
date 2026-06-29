#pragma once

#include <limits>
#include <string>
#include <vector>

#include <xoos/types/int.h>
#include <xoos/yc-decode/yc-decoder.h>

namespace xoos::alignment_metrics {

// We use 255 as the uninitialized value for min_quality because a valid Phred-scale quality score
// can never exceed 93 due to ASCII encoding limits, and 255 is outside this range.
constexpr u8 kHomopolymerUninitializedBaseQuality = std::numeric_limits<u8>::max();

constexpr char kAllBase = 'X';
constexpr char kHpStatUninitializedBase = '\0';

// A 32-bit bit vector that indicates whether an error is present at a given position in a homopolymer.
using HpErrorByPosition = u32;

/**
 * Class to track homopolymer error profiles for a single read.
 * It is created for each position in a read that is the first base of a homopolymer.
 * It records insertions, deletions, and substitutions at each position within a reference
 * homopolymer region that the read overlaps.
 */
class HomopolymerErrorProfileForRead {
 public:
  /**
   * Record a substitution at the given relative position within the homopolymer.
   * @param position The relative position within the homopolymer (0-indexed).
   */
  void AddSubstitution(u8 position);
  /**
   * Record an insertion at the given relative position within the homopolymer.
   * @param position The relative position within the homopolymer (0-indexed).
   * @param sequence The sequence of bases that were inserted.
   */
  void AddInsertion(u8 position, const std::string& sequence);
  /**
   * Record a deletion at the given relative position within the homopolymer.
   * @param position The relative position within the homopolymer (0-indexed).
   */
  void AddDeletion(u8 position);
  /**
   * Update the minimum quality score for the homopolymer.
   * @param quality The new minimum quality score.
   */
  void UpdateMinQuality(u8 quality);
  /**
   * Update the minimum base type for the homopolymer.
   * @param base_type The new minimum base type.
   */
  void UpdateMinBaseType(yc_decode::BaseType base_type);
  /**
   * Get the minimum quality score recorded for the homopolymer.
   * @return The minimum quality score.
   */
  u8 GetMinQuality() const;
  /**
   * Get the minimum base type recorded for the homopolymer.
   * @return The minimum base type.
   */
  yc_decode::BaseType GetMinBaseType() const;
  /**
   * Get a bit vector indicating the positions of insertions within the homopolymer.
   * @return A bit vector where each bit represents whether an insertion occurred at that position.
   */
  HpErrorByPosition GetInsertionByPosition() const;
  /**
   * Get a bit vector indicating the positions of deletions within the homopolymer.
   * @return A bit vector where each bit represents whether a deletion occurred at that position.
   */
  HpErrorByPosition GetDeletionByPosition() const;
  /**
   * Get a bit vector indicating the positions of insertions and deletions within the homopolymer.
   * @return A bit vector where each bit represents whether an insertion or deletion occurred at that position.
   */
  HpErrorByPosition GetIndelByPosition() const;
  /**
   * Get a bit vector indicating the positions of substitutions within the homopolymer.
   * @return A bit vector where each bit represents whether a substitution occurred at that position.
   */
  HpErrorByPosition GetSubstitutionByPosition() const;
  /**
   * Check if all insertions in the homopolymer are of the same base as the given base.
   * @param base The base to compare against (e.g., 'A', 'C', 'G', 'T').
   * @return True if all insertions are of the same base, false otherwise.
   */
  bool AreInsertionsHomogeneous(char base) const;
  // Indicates that it is the first base of a homopolymer region
  bool valid{false};

 private:
  HpErrorByPosition _insertion_by_position{0};
  HpErrorByPosition _deletion_by_position{0};
  HpErrorByPosition _substitution_by_position{0};
  HpErrorByPosition _indel_by_position{0};
  u8 _min_quality{kHomopolymerUninitializedBaseQuality};
  yc_decode::BaseType _min_base_type{yc_decode::BaseType::kConcordant};
  std::vector<std::string> _insertion_sequences{};
};

}  // namespace xoos::alignment_metrics
