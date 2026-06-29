#pragma once
#include <map>
#include <string>
#include <vector>

#include <xoos/io/metadata-util.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "core/hp-error-for-read.h"

namespace xoos::alignment_metrics {

/**
 * @brief Statistics for a particular homopolymer region of a specific base and length.
 *
 * Tracks comprehensive error and coverage statistics for homopolymer sequences (e.g., AAAA, GGGGG).
 * Each HpStat instance represents a specific homopolymer type defined by its base (A/C/G/T) and length.
 * The struct maintains various coverage metrics (total, spanning, effective) and error counts (insertions,
 * deletions, substitutions) aggregated across all reads that align to homopolymers of this type.
 * Position-specific error tracking allows detailed analysis of where errors occur within the homopolymer.
 * This is stored on the first position of the homopolymer in the reference.
 */
struct HpStat {
  HpStat() = default;

  // statistics for the homopolymer
  // number of reads that cover the first position of an HP
  u64 total_reads{};

  // number of reads that cover the first position of an HP with at least anchor_len flanking bases on both ends of the
  // read; such read is called a spanning read
  u64 spanning_reads{};

  // number of spanning reads that covers an HP with at least one discordant base
  u64 discordant_reads{};

  // number of spanning reads that covers an HP with at least one base that does not pass the
  // base quality and base type filters
  u64 low_quality_reads{};

  // number of spanning reads that covers an HP such that bases aligned to the HP region
  // (1) all pass the base quality filter;
  // (2) has no substitution error or mismatch policy is not “skip”;
  // (3) has only insertions of the same base as the HP base if cognate_ins = True
  u64 effective_reads{};

  u64 count_ins{};    // number of reads that contain an insertion for this type of HP
  u64 count_del{};    // number of reads that contain a deletion for this type of HP
  u64 count_indel{};  // number of reads that contain an insertion or deletion for this type of HP
  u64 count_sub{};    // number of reads that contain a substitution for this type of HP

  // A vector of bit vector for each position in the homopolymer and for each read
  // if the bit is set, then the read has an insertion at that position
  std::vector<HpErrorByPosition> insertion_by_position;
  // A vector of bit vector for each position in the homopolymer and for each read
  // if the bit is set, then the read has a deletion at that position
  std::vector<HpErrorByPosition> deletion_by_position;
  // A vector of bit vector for each position in the homopolymer and for each read
  // if the bit is set, then the read has a substitution at that position
  std::vector<HpErrorByPosition> substitution_by_position;
  // A vector of bit vector for each position in the homopolymer and for each read
  // if the bit is set, then the read has an insertion or deletion at that position
  std::vector<HpErrorByPosition> indel_by_position;

  // Total number of homopolymers of this type
  u64 hp_count{};
};

// key is base + hp_length
// We track the base as well as the length of the homopolymer because we want to be able to track the statistics for
struct HpStatKey {
  char base;
  u8 hp_length;

  std::strong_ordering operator<=>(const HpStatKey& other) const {
    if (hp_length == other.hp_length) {
      if (base == other.base) {
        return std::strong_ordering::equal;
      }
      return base <=> other.base;
    }
    return hp_length <=> other.hp_length;
  }
};

// A container for the homopolymer statistics that tracks the homopolymer statistics for a range of homopolymer lengths
// and bases.
struct HpStats {
  explicit HpStats(u8 max_hp_length);

  std::map<HpStatKey, HpStat> hp_stats{};

  static std::vector<std::string> GetHeaders();

  // Sums all the hp stats for a given hp_length across all bases
  HpStat GetTotalHpStat(u8 hp_length) const;

  // Increments the HP count for a given key
  void IncrementHpCount(const HpStatKey& key);

  // Writes to the given tsv file
  void WriteTsv(const fs::path& output_path, const io::Comments& comments) const;

  friend HpStats operator+(const HpStats& lhs, const HpStats& rhs);

  HpStats& operator+=(const HpStats& other);
};
}  // namespace xoos::alignment_metrics
