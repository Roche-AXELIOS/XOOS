#pragma once
#include <optional>
#include <string_view>
#include <tuple>
#include <vector>

#include <xoos/types/int.h>
#include <xoos/yc-decode/yc-decoder.h>

#include "alignment-metrics-options.h"
#include "core/super-region.h"

namespace xoos::alignment_metrics {

// Start (0-inclusive) and end (0-inclusive) positions of a homopolymer along with its base
using Homopolymers = std::vector<std::tuple<u64, u64, char>>;

/**
 * Struct to hold both modified and unmodified base qualities and types for a read.
 */
struct ReadQualities {
  /**
   * Original base qualities from the read
   */
  std::vector<u8> qualities{};
  /**
   * Modified base qualities. For simplex reads, this is modified based on homopolymer masking.
   * For duplex reads, this is the same as the original base qualities unless the user specifies
   * to always mask homopolymers by base quality.
   */
  std::vector<u8> modified_qualities{};
  /**
   * Unmodified base types from the read. For reads without a YC tag, this will be all simplex.
   */
  std::vector<yc_decode::BaseType> base_types{};
  /**
   * Modified base types. For duplex reads, this is modified based on homopolymer masking.
   * For reads without a YC tag, this will be all simplex.
   */
  std::vector<yc_decode::BaseType> modified_base_types{};
};

/**
 * Finds homopolymers of length `min_hp_length` or greater in the given sequence, maximizing the size of each.
 * If the homopolymer lengths are greater than `max_hp_length`, they are not included in the result.
 *
 * The function iterates through the sequence between `alignment_start` and `alignment_end`, tracking starts
 * and ends of homopolymer regions, and adds them to the homopolymers collection. If the current base is the
 * same as the previous one, the end of the homopolymer is extended. If the current base is different from
 * the previous one, the homopolymer is added to the collection if it has a length greater than 1, and the
 * start and end positions are reset.
 */
Homopolymers FindHomopolymers(std::string_view sequence,
                              u16 min_hp_length,
                              u16 max_hp_length,
                              f64 hp_subsampling_fraction = kDefaultHpSubsamplingFraction,
                              std::optional<s32> hp_subsampling_seed = std::nullopt);

/**
 * Finds homopolymers of length `min_hp_length` or greater in the given sequence, maximizing the size of each.
 * If `max_hp_length` is specified and the homopolymer lengths are greater than `max_hp_length`, they are not
 * included in the result.
 *
 * The function iterates through the sequence between `alignment_start` and `alignment_end`, tracking starts
 * and ends of homopolymer regions, and adds them to the homopolymers collection. If the current base is the
 * same as the previous one, the end of the homopolymer is extended. If the current base is different from
 * the previous one, the homopolymer is added to the collection if it has a length greater than 1, and the
 * start and end positions are reset.
 */
Homopolymers FindReadHomopolymers(const bam1_t* bam1_ptr,
                                  u16 min_hp_length,
                                  std::optional<u16> max_hp_length = std::nullopt);
/**
 * Find homopolymers in the given sequence within the super region.
 */
Homopolymers FindReferenceHomopolymers(std::string_view sequence,
                                       u16 min_hp_length,
                                       u16 max_hp_length,
                                       const SuperRegion& region,
                                       f64 hp_subsampling_fraction = kDefaultHpSubsamplingFraction,
                                       std::optional<s32> hp_subsampling_seed = std::nullopt);

/**
 * Modify the base qualities for homopolymer regions in the read.
 * If any base in the homopolymer has a base quality of 5 or less, then all bases in the homopolymer
 * are set to the minimum base quality in that homopolymer.
 *
 * This is only applied to reads without a YC tag (i.e. simplex reads) or when the user specifies to always
 * mask homopolymers by base quality.
 */
void ModifyBaseQualityForHomopolymers(ReadQualities& read_qualities,
                                      const Homopolymers& homopolymers,
                                      u8 base_quality_threshold_for_hp_masking);

/**
 * Modify the base type for homopolymers regions in a duplex read.
 * If any base in the homopolymer has a base type of discordant, then all bases in the homopolymer
 * are set to discordant.
 *
 * This is only applied to reads with a YC tag.
 */
void ModifyBaseTypeForHomopolymers(ReadQualities& read_qualities, const Homopolymers& homopolymers);

/**
 * Mask the homopolymers in the given read by modifying the base qualities and/or types
 *
 * If the read has a YC tag, we will modify the base types. Otherwise, we will modify the base qualities.
 *
 * @param bam_ptr BAM record to mask homopolymers for
 * @param disable_hp_masking If true, do not mask homopolymer regions
 * @param disable_base_type_decoding If true, skip decoding the YC tag and treat the read as simplex.
 * @param base_quality_threshold_for_hp_masking Base quality threshold for masking homopolymers. If the minimum base
 * quality in the homopolymer is less than or equal to this threshold, we will mask the homopolymer by setting all bases
 * in the homopolymer to the minimum base quality in that homopolymer.
 * @return Modified and unmodified base qualities and types for the read
 */
ReadQualities MaskReadHomopolymers(const bam1_t* bam_ptr,
                                   bool disable_hp_masking = false,
                                   bool disable_base_type_decoding = false,
                                   u8 base_quality_threshold_for_hp_masking = kDefaultBaseQualityThresholdForHpMasking);

}  // namespace xoos::alignment_metrics
