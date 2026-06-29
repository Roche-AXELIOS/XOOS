#include "util/homopolymer-util.h"

#include <chrono>
#include <optional>
#include <random>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

#include <htslib/sam.h>

#include <xoos/error/error.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/int.h>
#include <xoos/util/sequence-functions.h>
#include <xoos/yc-decode/yc-decoder.h>

#include "alignment-metrics-options.h"
#include "core/super-region.h"

namespace xoos::alignment_metrics {

constexpr s32 kNBaseInt = 15;

Homopolymers FindHomopolymers(const std::string_view sequence,
                              const u16 min_hp_length,
                              const u16 max_hp_length,
                              const f64 hp_subsampling_fraction,
                              const std::optional<s32> hp_subsampling_seed) {
  Homopolymers homopolymers;
  u64 hp_start = 0;
  u64 hp_end = 0;

  const s32 seed = hp_subsampling_seed.value_or(static_cast<s32>(
      std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count()));
  std::mt19937 rng(seed);
  std::uniform_real_distribution dist(0.0, 1.0);

  const u64 alignment_end = sequence.size();
  // start at position 1
  for (u64 i = 1; i < alignment_end;) {
    // if we find a beginning of a homopolymer, then we need to find the end of it
    // We only support ACTG (no soft-masked bases)
    if (const char seq_at_i = sequence.at(i); seq_at_i == sequence.at(i - 1) && sequence::IsACGT(seq_at_i)) {
      // Current base is the same as the previous one; extend the end position of the homopolymer
      const char current_base = sequence.at(i);
      while (sequence[i] == current_base && i < alignment_end) {
        hp_end = i;
        ++i;
      }
      // both are inclusive, so we need to subtract 1
      if (hp_end - hp_start + 1 >= min_hp_length && hp_end - hp_start + 1 <= max_hp_length) {
        if (hp_subsampling_fraction < 1.0 && dist(rng) >= hp_subsampling_fraction) {
          continue;
        }
        homopolymers.emplace_back(hp_start, hp_end, current_base);
      }
    } else {
      hp_start = i;
      hp_end = i;
      ++i;
    }
  }
  return homopolymers;
}

Homopolymers FindReferenceHomopolymers(const std::string_view sequence,
                                       const u16 min_hp_length,
                                       const u16 max_hp_length,
                                       const SuperRegion& region,
                                       const f64 hp_subsampling_fraction,
                                       const std::optional<s32> hp_subsampling_seed) {
  Homopolymers homopolymers;
  for (const auto& subregion : region.subregions) {
    const auto subregion_sequence = sequence.substr(static_cast<size_t>(subregion.start - region.start),
                                                    static_cast<size_t>(subregion.end - subregion.start));
    const auto subregion_homopolymers = FindHomopolymers(
        subregion_sequence, min_hp_length, max_hp_length, hp_subsampling_fraction, hp_subsampling_seed);
    for (const auto& [start, end, base] : subregion_homopolymers) {
      homopolymers.emplace_back(ToUnsigned(subregion.start) + start, ToUnsigned(subregion.start) + end, base);
    }
  }
  return homopolymers;
}

Homopolymers FindReadHomopolymers(const bam1_t* bam1_ptr, u16 min_hp_length, std::optional<u16> max_hp_length) {
  Homopolymers homopolymers;
  u64 hp_start = 0;
  u64 hp_end = 0;
  const u64 alignment_end = bam1_ptr->core.l_qseq;
  const u8* seq = bam_get_seq(bam1_ptr);
  // start at position 1
  for (u64 i = 1; i < alignment_end;) {
    // if we find a beginning of a homopolymer, then we need to find the end of it
    const auto base0 = bam_seqi(seq, i);
    // exclude 'N' bases
    if (base0 == bam_seqi(seq, i - 1) && base0 != kNBaseInt) {
      // Current base is the same as the previous one; extend the end position of the homopolymer
      const auto current_base = base0;
      while (bam_seqi(seq, i) == current_base && i < alignment_end) {
        hp_end = i;
        ++i;
      }
      // both are inclusive, so we need to subtract 1
      // if max_hp_length is not set, then we don't need to check the length
      // otherwise we need to check if the length is less than or equal to the max_hp_length
      if (hp_end - hp_start + 1 >= min_hp_length &&
          (!max_hp_length.has_value() || hp_end - hp_start + 1 <= max_hp_length.value())) {
        homopolymers.emplace_back(hp_start, hp_end, seq_nt16_str[current_base]);
      }
    } else {
      hp_start = i;
      hp_end = i;
      ++i;
    }
  }
  return homopolymers;
}

void ModifyBaseQualityForHomopolymers(ReadQualities& read_qualities,
                                      const Homopolymers& homopolymers,
                                      const u8 base_quality_threshold_for_hp_masking) {
  for (const auto& [start, end, _] : homopolymers) {
    // Find the minimum quality in the homopolymer region
    const u8 min_quality = *std::min_element(read_qualities.qualities.begin() + ToSigned(start),
                                             read_qualities.qualities.begin() + ToSigned(end) + 1);
    if (min_quality <= base_quality_threshold_for_hp_masking && start <= end &&
        end < read_qualities.modified_qualities.size()) {
      std::fill(read_qualities.modified_qualities.begin() + ToSigned(start),
                read_qualities.modified_qualities.begin() + ToSigned(end) + 1,
                min_quality);
    }
  }
}

void ModifyBaseTypeForHomopolymers(ReadQualities& read_qualities, const Homopolymers& homopolymers) {
  for (const auto& [start, end, _] : homopolymers) {
    const yc_decode::BaseType min_base_type = *std::min_element(read_qualities.base_types.begin() + ToSigned(start),
                                                                read_qualities.base_types.begin() + ToSigned(end) + 1);
    if (min_base_type == yc_decode::BaseType::kDiscordant && start <= end &&
        end < read_qualities.modified_base_types.size()) {
      std::fill(read_qualities.modified_base_types.begin() + ToSigned(start),
                read_qualities.modified_base_types.begin() + ToSigned(end) + 1,
                min_base_type);
    }
  }
}

static ReadQualities CreateMappedQualities(const bam1_t* bam, bool disable_base_type_decoding) {
  ReadQualities read_qualities;
  // Create a vector of modified qualities and unadjusted qualities
  std::vector<u8> unadjusted_qualities(bam->core.l_qseq, 0);

  // Get the unadjusted qualities from the BAM record
  const u8* unadjusted_qual = bam_get_qual(bam);
  for (u32 i = 0; std::cmp_less(i, bam->core.l_qseq); i++) {
    unadjusted_qualities[i] = unadjusted_qual[i];
  }

  // If the read is duplex, we need to adjust the base type for bases.
  // We don't know if the bam record is duplex so we explicitly deconstruct the YC tag and if empty, we assume simplex.
  std::vector<yc_decode::BaseType> base_types(unadjusted_qualities.size(), yc_decode::BaseType::kSimplex);
  if (!disable_base_type_decoding) {
    try {
      const auto yc_tag = yc_decode::DeserializeYcTag(bam);
      base_types = yc_tag.GetBaseTypes();
    } catch (const std::runtime_error& e) {
      throw error::Error(
          "{}. "
          "To bypass YC tag deserialization, set the flag `--disable-base-type-decoding`.",
          e.what());
    }
    // If no base types were decoded, we assume simplex
    if (base_types.empty()) {
      base_types.resize(unadjusted_qualities.size(), yc_decode::BaseType::kSimplex);
    }
    // If some base types were decoded but the length doesn't match the read length, we throw an error
    // because it indicates a malformed YC tag
    try {
      yc_decode::RequireConsistentReadLength(bam, base_types);
    } catch (const std::runtime_error& e) {
      throw error::Error(
          "{}. "
          "To bypass YC tag deserialization, set the flag `--disable-base-type-decoding`.",
          e.what());
    }
  }
  read_qualities.base_types = base_types;
  read_qualities.modified_base_types = base_types;
  read_qualities.qualities = std::move(unadjusted_qualities);
  read_qualities.modified_qualities = read_qualities.qualities;
  return read_qualities;
}

ReadQualities MaskReadHomopolymers(const bam1_t* bam_ptr,
                                   const bool disable_hp_masking,
                                   const bool disable_base_type_decoding,
                                   const u8 base_quality_threshold_for_hp_masking) {
  auto read_qualities = CreateMappedQualities(bam_ptr, disable_base_type_decoding);
  if (disable_hp_masking) {
    return read_qualities;
  }
  // Make a pass through the sequence to find all homopolymer positions
  // We do not specify a max_hp_length here because we want to find all homopolymers
  // for the purpose of quality modification
  const Homopolymers homopolymers = FindReadHomopolymers(bam_ptr, kHpMinLength);
  if (!homopolymers.empty()) {
    // Only modify base quality if there is no base type information (i.e. no YC tag)
    // This is to avoid removing additional HPs with calibrated base qualities
    if (bam_aux_get(bam_ptr, "YC") == nullptr || disable_base_type_decoding) {
      ModifyBaseQualityForHomopolymers(read_qualities, homopolymers, base_quality_threshold_for_hp_masking);
    } else {
      ModifyBaseTypeForHomopolymers(read_qualities, homopolymers);
    }
  }
  return read_qualities;
}

}  // namespace xoos::alignment_metrics
