#include "read-processing.h"

#include <numeric>
#include <optional>

#include <xoos/io/htslib-util/htslib-util.h>

#include "core/config.h"
#include "core/filtering.h"
#include "read-name-utils.h"
#include "util/log-util.h"
#include "util/region-util.h"
#include "util/seq-util.h"

namespace xoos::svc {

using enum AlignOp;
using yc_decode::BaseType;

/**
 * @brief Determine if the provided read_id has been seen before, if not add it to the list of seen read_ids.
 * @return true if read_ids contains read_id, false otherwise and read_id is added to read_ids
 */
static bool ContainsReadIdInsertIfNot(ReadIds& read_ids, const ReadId read_id) {
  // Read ids are inserted at the back of the vector, therefore we iterate in reverse
  // as we are more likely to see recent read ids earlier.
  if (std::find(read_ids.crbegin(), read_ids.crend(), read_id) != read_ids.crend()) {
    return true;
  }
  read_ids.emplace_back(read_id);
  return false;
}

std::string Get1bpContext(const std::string& ref_seq, const u64 pos) {
  if (pos == 0 || pos + 1 >= ref_seq.size()) {
    // context is out of bounds for the reference sequence
    return "";
  }
  // extract the 1-bp context around the specified position
  std::string context;
  context.push_back(ref_seq[pos - 1]);
  context.push_back(ref_seq[pos + 1]);
  if (IsAnyNotACTG(context)) {
    return "";
  }
  return context;
}

void UpdateDerivedFeatureValues(UnifiedVariantFeature& feature) {
  if (feature.support > 0) {
    const auto support = static_cast<f64>(feature.support);
    feature.mapq_lt60_ratio = feature.mapq_lt60_count / support;
    feature.mapq_lt40_ratio = feature.mapq_lt40_count / support;
    feature.mapq_lt30_ratio = feature.mapq_lt30_count / support;
    feature.mapq_lt20_ratio = feature.mapq_lt20_count / support;
    feature.familysize_lt5_ratio = feature.familysize_lt5_count / support;
    feature.familysize_lt3_ratio = feature.familysize_lt3_count / support;
    feature.baseq_lt20_ratio = feature.baseq_lt20_count / support;
    feature.distance_mean = static_cast<f64>(feature.distance_sum) / support;
    feature.baseq_mean = feature.baseq_sum / support;
    feature.mapq_mean = feature.mapq_sum / support;
    feature.familysize_mean = feature.familysize_sum / support;
  }
  // for v9.2 chemistry, `duplex` consists of plus and minus reads. Thus, `duplex` is 2x whereas `nonduplex` is 1x
  feature.weighted_score = feature.nonduplex + (2.0 * feature.duplex);
  const auto tmp = (feature.duplex * 0.5 + feature.plusonly) / (feature.duplex + feature.plusonly + feature.minusonly);
  feature.strandbias = tmp * (1 - tmp);

#ifdef SOMATIC_ENABLE
  if (feature.tumor_support > 0) {
    feature.tumor_baseq_mean = feature.tumor_baseq_sum / feature.tumor_support;
    feature.tumor_mapq_mean = static_cast<f64>(feature.tumor_mapq_sum) / feature.tumor_support;
    feature.tumor_distance_mean = static_cast<f64>(feature.tumor_distance_sum) / feature.tumor_support;
  }
  if (feature.normal_support > 0) {
    feature.normal_baseq_mean = feature.normal_baseq_sum / feature.normal_support;
    feature.normal_mapq_mean = static_cast<f64>(feature.normal_mapq_sum) / feature.normal_support;
    feature.normal_distance_mean = static_cast<f64>(feature.normal_distance_sum) / feature.normal_support;
  }
#endif  // SOMATIC_ENABLE

  if (feature.duplex_lowbq > 0) {
    // `duplex_lowbq` is incremented by 0.5 for each `lowbq` read, so multiply it by 2 to get the actual count
    const auto duplex_lowbq = 2 * feature.duplex_lowbq;
    feature.mapq_mean_lowbq = feature.mapq_sum_lowbq / duplex_lowbq;
    feature.distance_mean_lowbq = static_cast<f64>(feature.distance_sum_lowbq) / duplex_lowbq;
  }
  if (feature.simplex > 0) {
    feature.mapq_mean_simplex = static_cast<f64>(feature.mapq_sum_simplex) / feature.simplex;
    feature.distance_mean_simplex = static_cast<f64>(feature.distance_sum_simplex) / feature.simplex;
  }
}

void UpdateDerivedFeatureValues(UnifiedReferenceFeature& feature) {
  if (feature.support > 0) {
    const auto support = static_cast<f64>(feature.support);
    feature.mapq_lt60_ratio = feature.mapq_lt60_count / support;
    feature.mapq_lt40_ratio = feature.mapq_lt40_count / support;
    feature.mapq_lt30_ratio = feature.mapq_lt30_count / support;
    feature.mapq_lt20_ratio = feature.mapq_lt20_count / support;
    feature.mapq_mean = feature.mapq_sum / support;
    feature.baseq_mean = feature.baseq_sum / support;
    feature.baseq_lt20_ratio = feature.baseq_lt20_count / support;
    feature.familysize_mean = feature.familysize_sum / support;
    feature.familysize_lt5_ratio = feature.familysize_lt5_count / support;
    feature.familysize_lt3_ratio = feature.familysize_lt3_count / support;
    feature.distance_mean = static_cast<f64>(feature.distance_sum) / support;
    feature.nonhomopolymer_mapq_lt60_ratio = feature.nonhomopolymer_mapq_lt60_count / support;
    feature.nonhomopolymer_mapq_lt40_ratio = feature.nonhomopolymer_mapq_lt40_count / support;
    feature.nonhomopolymer_mapq_lt30_ratio = feature.nonhomopolymer_mapq_lt30_count / support;
    feature.nonhomopolymer_mapq_lt20_ratio = feature.nonhomopolymer_mapq_lt20_count / support;
    feature.nonhomopolymer_baseq_mean = feature.nonhomopolymer_baseq_sum / support;
    feature.nonhomopolymer_mapq_mean = feature.nonhomopolymer_mapq_sum / support;
  }
  if (feature.duplex_lowbq > 0) {
    const auto duplex_lowbq = 2 * feature.duplex_lowbq;
    feature.mapq_mean_lowbq = feature.mapq_sum_lowbq / duplex_lowbq;
    feature.distance_mean_lowbq = static_cast<f64>(feature.distance_sum_lowbq) / duplex_lowbq;
  }
  if (feature.simplex > 0) {
    const auto simplex = static_cast<f64>(feature.simplex);
    feature.mapq_mean_simplex = feature.mapq_sum_simplex / simplex;
    feature.distance_mean_simplex = static_cast<f64>(feature.distance_sum_simplex) / simplex;
  }
}

// Feature thresholds for counting
static constexpr u32 kFamilySizeCountThreshold5 = 5;
static constexpr u32 kFamilySizeCountThreshold3 = 3;
static constexpr u32 kMapqCountThreshold60 = 60;
static constexpr u32 kMapqCountThreshold40 = 40;
static constexpr u32 kMapqCountThreshold30 = 30;
static constexpr u32 kMapqCountThreshold20 = 20;
static constexpr u32 kBaseqCountThreshold20 = 20;

void IncrementFeature(UnifiedVariantFeature& feature,
                      const AlignContext& align_ctx,
                      const f64 baseq,
                      const u64 distance,
                      const bool near_non_concordant) {
  if (ContainsReadIdInsertIfNot(feature.read_ids, align_ctx.read_id)) {
    // Skip previously processed reads to avoid double counting
    return;
  }

  const u8 mapq{align_ctx.bam_record->core.qual};
  if (align_ctx.read_type == ReadType::kSimplex) {
    ++feature.simplex;
    feature.mapq_sum_simplex += mapq;
    feature.distance_sum_simplex += distance;
    return;
  }

  if (near_non_concordant) {
    if (align_ctx.plus_counts > 0 && align_ctx.minus_counts > 0) {
      feature.duplex_lowbq += 0.5;
    }
    feature.mapq_sum_lowbq += mapq;
    feature.distance_sum_lowbq += distance;
    return;
  }

  if (align_ctx.plus_counts > 0 && align_ctx.minus_counts > 0) {
    ++feature.duplex;
  } else {
    ++feature.nonduplex;
    if (align_ctx.plus_counts > 0) {
      ++feature.plusonly;
    }
    if (align_ctx.minus_counts > 0) {
      ++feature.minusonly;
    }
  }

  const u32 familysize{align_ctx.plus_counts + align_ctx.minus_counts};
  feature.familysize_sum += familysize;
  if (familysize < kFamilySizeCountThreshold5) {
    ++feature.familysize_lt5_count;
  }
  if (familysize < kFamilySizeCountThreshold3) {
    ++feature.familysize_lt3_count;
  }

  feature.mapq_sum += mapq;
  if (mapq < kMapqCountThreshold60) {
    ++feature.mapq_lt60_count;
  }
  if (mapq < kMapqCountThreshold40) {
    ++feature.mapq_lt40_count;
  }
  if (mapq < kMapqCountThreshold30) {
    ++feature.mapq_lt30_count;
  }
  if (mapq < kMapqCountThreshold20) {
    ++feature.mapq_lt20_count;
  }

  feature.baseq_sum += baseq;
  if (baseq < kBaseqCountThreshold20) {
    ++feature.baseq_lt20_count;
  }

  feature.distance_sum += distance;

  if (feature.support == 0) {
    feature.mapq_min = mapq;
    feature.mapq_max = mapq;
    feature.baseq_min = baseq;
    feature.baseq_max = baseq;
    feature.distance_min = distance;
    feature.distance_max = distance;
  } else {
    feature.mapq_min = std::min(mapq, feature.mapq_min);
    feature.mapq_max = std::max(mapq, feature.mapq_max);
    feature.baseq_min = std::min(baseq, feature.baseq_min);
    feature.baseq_max = std::max(baseq, feature.baseq_max);
    feature.distance_min = std::min(distance, feature.distance_min);
    feature.distance_max = std::max(distance, feature.distance_max);
  }

  ++feature.support;
  if (bam_is_rev(align_ctx.bam_record)) {
    ++feature.support_reverse;
  }

  feature.weighted_depth += baseq / kMaxPossibleBaseQual * mapq / kMaxPossibleMapQual;
#ifdef SOMATIC_ENABLE
  if (align_ctx.is_tumor_sample.has_value()) {
    if (align_ctx.is_tumor_sample.value()) {
      ++feature.tumor_support;
      feature.tumor_baseq_sum += baseq;
      feature.tumor_mapq_sum += mapq;
      feature.tumor_distance_sum += distance;
      if (is_reverse) {
        ++feature.tumor_support_reverse;
      }
    } else {
      ++feature.normal_support;
      feature.normal_baseq_sum += baseq;
      feature.normal_mapq_sum += mapq;
      feature.normal_distance_sum += distance;
      if (is_reverse) {
        ++feature.normal_support_reverse;
      }
    }
  }
#endif  // SOMATIC_ENABLE
}

vec<AlignOpInfo> GetAlignOpInfos(const bam1_t* record, const std::string& target_seq) {
  vec<AlignOpInfo> align_infos;
  align_infos.reserve(record->core.n_cigar);

  // The 0-based position in the query sequence
  u64 q_pos = 0;
  // The 0-based position in the target sequence
  auto t_pos = static_cast<u64>(record->core.pos);

  const u32* cigar = bam_get_cigar(record);
  const u8* seq = bam_get_seq(record);

  // Lamda function to store the match or mismatch information in the align_infos vector and update the query and target
  // positions accordingly.
  // If the previous and the current base are both matches, we increment the length of the previous operation.
  // This is done to reduce the number of operations stored in the align_infos vector.
  // Note that mismatches are always stored as separate operations, even if they are consecutive.
  auto handle_match = [&seq, &q_pos, &target_seq, &t_pos, &align_infos](const u64 op_len) {
    for (u64 j = 0; j < op_len; ++j) {
      const char q_base = seq_nt16_str[bam_seqi(seq, q_pos + j)];
      const char t_base = target_seq.at(t_pos + j);
      auto align_info_op = q_base == t_base ? kMatch : kMismatch;
      if (!align_infos.empty() && align_infos.back().op == align_info_op && align_info_op == kMatch) {
        ++align_infos.back().op_len;
      } else {
        align_infos.emplace_back(align_info_op, 1, t_pos + j, q_pos + j);
      }
    }
    q_pos += op_len;
    t_pos += op_len;
  };

  for (u32 i = 0; i < record->core.n_cigar; ++i) {
    const auto op = bam_cigar_op(cigar[i]);
    const u64 op_len = bam_cigar_oplen(cigar[i]);
    const bool not_first_or_last_op = !align_infos.empty() && i < record->core.n_cigar - 1;
    switch (op) {
      case BAM_CMATCH: {
        handle_match(op_len);
        break;
      }
      case BAM_CINS: {
        // skip insertion if it is either the first or last op
        if (not_first_or_last_op) {
          // `reference_position - 1` because an insertion starts at the previous matching base
          align_infos.emplace_back(kInsertion, op_len, t_pos - 1, q_pos);
        }
        q_pos += op_len;
        break;
      }
      case BAM_CDEL: {
        // skip deletion if it is either the first or last op
        if (not_first_or_last_op) {
          align_infos.emplace_back(kDeletion, op_len, t_pos, q_pos);
        }
        t_pos += op_len;
        break;
      }
      case BAM_CREF_SKIP: {
        align_infos.emplace_back(kReferenceSkip, op_len, t_pos, q_pos);
        t_pos += op_len;
        break;
      }
      case BAM_CSOFT_CLIP: {
        // do not store information about the soft clip.
        q_pos += op_len;
        break;
      }
      case BAM_CEQUAL: {
        align_infos.emplace_back(kMatch, op_len, t_pos, q_pos);
        q_pos += op_len;
        t_pos += op_len;
        break;
      }
      case BAM_CDIFF: {
        align_infos.emplace_back(kMismatch, op_len, t_pos, q_pos);
        q_pos += op_len;
        t_pos += op_len;
        break;
      }
      default:
        // Hard clip and padding do not consume query or reference sequence bases.
        // So, query or reference position are not incremented.
        // Do not store information about them because they do not provide any useful information for
        // feature extraction.
        break;
    }
  }
  return align_infos;
}

// A duplex hairpin read has a family size of 2, i.e. 1 plus read and 1 minus read.
static constexpr u32 kDuplexPlusCounts{1};
static constexpr u32 kDuplexMinusCounts{1};

AlignContext::AlignContext(const ComputeBamFeaturesParams& params,
                           const vec<AlignOpInfo>& align_op_infos,
                           const bam1_t* bam_record,
                           const Region& region,
                           const std::string& ref_seq,
                           const ReadId read_id,
                           const PositionToVcfFeaturesMap& vcf_feats)
    : params(params),
      align_op_infos(align_op_infos),
      bam_record(bam_record),
      region(region),
      ref_seq(ref_seq),
      read_id(read_id),
      vcf_feats(vcf_feats),
      has_vcf_feats(!vcf_feats.empty()) {
  u32 family_size = 0;
  if (params.duplex) {
    plus_counts = kDuplexPlusCounts;
    minus_counts = kDuplexMinusCounts;
    family_size = kDuplexReadFamilySize;
  } else {
    ParseReadName(bam_get_qname(bam_record), plus_counts, minus_counts, family_size);
  }
  ref_end = align_op_infos.back().ref_pos + align_op_infos.back().op_len;
  // find the homopolymer overlapping the last aligned reference base position, which is `ref_end - 1`
  homopolymer_start =
      FindOverlappingHomopolymer(ref_seq, ref_end - 1, params.min_homopolymer_length, align_op_infos.front().ref_pos);

  if (params.decode_yc) {
    const auto& yc_tag = yc_decode::DeserializeYcTag(bam_record);
    base_types = yc_tag.GetBaseTypes();
    if (base_types.empty()) {
      read_type = ReadType::kSimplex;
    } else {
      read_type = ReadType::kDuplex;
      near_non_concordant = FindNearbyNonConcordantBase(bam_record, base_types);
    }
  } else {
    if (family_size < params.min_family_size) {
      skip_read = true;
      return;
    }
    if (params.duplex) {
      read_type = ReadType::kDuplex;
      near_non_concordant = FindNearbyNonConcordantBase(bam_record);
    } else {
      read_type = ReadType::kUmiConsensus;
    }
  }
  if (!base_types.empty() && params.min_base_type > BaseType::kDiscordant) {
    min_base_type = params.min_base_type;
  }
  if (params.tumor_read_group.has_value()) {
    is_tumor_sample = (params.tumor_read_group.value() == io::BamAuxGet<std::string>(bam_record, "RG"));
  }
}

bool AlignContext::CanProcessInsertion(const AlignOpInfo& info) const {
  if (skip_read) {
    return false;
  }
  if (info.ref_pos < region.start || info.ref_pos >= region.end) {
    // insertion out of bounds
    return false;
  }
  if (params.filter_homopolymer && homopolymer_start.has_value() && homopolymer_start.value() < info.ref_pos) {
    // insertion is in a homopolymer region
    return false;
  }
  if (has_vcf_feats && !vcf_feats.contains(info.ref_pos)) {
    // VCF features were extracted, but this position has no VCF features.
    return false;
  }
  if (min_base_type.has_value()) {
    // check if any of the bases in the insertion meet the minimum base type requirement
    return std::any_of(base_types.begin() + static_cast<s64>(info.read_pos),
                       base_types.begin() + static_cast<s64>(info.read_pos + info.op_len),
                       [this](const BaseType bt) { return bt >= min_base_type.value(); });
  }
  return true;
}

bool AlignContext::CanProcessDeletion(const AlignOpInfo& info) const {
  if (skip_read) {
    return false;
  }
  if (!IntervalOverlap(region.start, region.end, info.ref_pos, info.ref_pos + info.op_len)) {
    // deletion does not overlap with the region
    return false;
  }
  if (region.prev_interval.has_value() &&
      IntervalOverlap(region.prev_interval.value(), info.ref_pos, info.ref_pos + info.op_len)) {
    // deletion overlaps with the previous interval
    // it should be processed in that region instead of this region
    return false;
  }
  if (params.filter_homopolymer && homopolymer_start.has_value() && homopolymer_start.value() <= info.ref_pos) {
    // deletion is in a homopolymer region
    return false;
  }
  if (has_vcf_feats && !vcf_feats.contains(info.ref_pos - 1)) {
    // VCF features were extracted, but the previous position has no VCF features.
    // We use the previous position because a deletion starts at the previous matching base.
    return false;
  }
  // check if the base type of the base before the deletion meets the minimum base type requirement
  return !min_base_type.has_value() || base_types.at(info.read_pos - 1) >= min_base_type.value();
}

bool AlignContext::CanProcessMismatch(const AlignOpInfo& info) const {
  if (skip_read) {
    return false;
  }
  if (info.ref_pos < region.start || info.ref_pos >= region.end) {
    // mismatch out of bounds
    return false;
  }
  if (params.filter_homopolymer && homopolymer_start.has_value() && homopolymer_start.value() <= info.ref_pos) {
    // mismatch is in a homopolymer region
    return false;
  }
  if (has_vcf_feats && !vcf_feats.contains(info.ref_pos)) {
    // VCF features were extracted, but this position has no VCF features
    return false;
  }
  // check if the base type of the base at the mismatch position meets the minimum base type requirement
  return !min_base_type.has_value() || base_types.at(info.read_pos) >= min_base_type.value();
}

bool AlignContext::CanProcessMatch(const u64 read_pos, const u64 ref_pos) const {
  if (skip_read) {
    return false;
  }
  if (ref_pos > region.end) {
    // match out of bounds
    return false;
  }
  if (has_vcf_feats && !vcf_feats.contains(ref_pos) && !vcf_feats.contains(ref_pos - 1)) {
    // VCF features were extracted, but the current and previous positions have no VCF features.
    // We check the previous position because a deletion starts at the previous matching base.
    return false;
  }
  // check if the base type of the base at the match position meets the minimum base type requirement
  return !min_base_type.has_value() || base_types.at(read_pos) >= min_base_type.value();
}

bool AlignContext::IsInsertionNearNonConcordantBase(const AlignOpInfo& info) const {
  if (near_non_concordant.empty()) {
    return false;
  }
  for (u64 i = 0; i < info.op_len; ++i) {
    if (near_non_concordant.at(info.read_pos + i)) {
      return true;
    }
  }
  return false;
}

bool AlignContext::IsDeletionNearNonConcordantBase(const AlignOpInfo& info) const {
  return !near_non_concordant.empty() ? near_non_concordant.at(info.read_pos) : false;
}

bool AlignContext::IsMismatchNearNonConcordantBase(const AlignOpInfo& info) const {
  return !near_non_concordant.empty() ? near_non_concordant.at(info.read_pos) : false;
}

bool AlignContext::IsMatchNearNonConcordantBase(const u64 read_pos) const {
  return !near_non_concordant.empty() ? near_non_concordant.at(read_pos) : false;
}

void ProcessInsertion(UnifiedVariantFeatures& features, const AlignContext& align_ctx, const AlignOpInfo& info) {
  auto* begin = bam_get_qual(align_ctx.bam_record) + info.read_pos;
  auto* end = begin + info.op_len;
  const u32 baseq_sum = std::accumulate(begin, end, u32{0});
  const f32 baseq_mean = static_cast<f32>(baseq_sum) / static_cast<f32>(info.op_len);
  if (baseq_mean < static_cast<f32>(align_ctx.params.min_bq)) {
    return;
  }
  if (align_ctx.bam_record->core.pos < 0) {
    // Should not get here, check should be done in upstream call for this. Keep just incase call made without check
    WarnAsErrorIfSet("Found alignment with start position < 0: '{}'", bam_get_qname(align_ctx.bam_record));
    return;
  }
  const auto align_pos = static_cast<u64>(align_ctx.bam_record->core.pos);
  const auto min_distance = std::min(info.ref_pos - align_pos, align_ctx.ref_end - info.ref_pos - 1);
  if (static_cast<f32>(min_distance) < align_ctx.params.min_allowed_distance_from_end) {
    return;
  }

  // no `-1` in position here because insertion starts at the previous matching base
  const char ref_base{align_ctx.ref_seq[info.ref_pos]};
  if (IsNotACTG(ref_base)) {
    // Since we are not confident in variant calls at non-ACGT reference positions,
    // we do not extract variant features at this position.
    return;
  }

  auto* const begin_seq = bam_get_seq(align_ctx.bam_record);
  std::string alt_allele(info.op_len + 1, 'N');
  alt_allele[0] = ref_base;
  for (u64 i = 0; i < info.op_len; ++i) {
    const char base = seq_nt16_str[bam_seqi(begin_seq, info.read_pos + i)];
    if (IsNotACTG(base)) {
      return;
    }
    alt_allele[i + 1] = base;
  }

  const std::string ref_allele{ref_base};
  // no `-1` in position here because insertion starts at the previous matching base
  const VariantId id(align_ctx.region.chrom, info.ref_pos, ref_allele, alt_allele);

  UnifiedVariantFeature& feat = features[id];
  if (feat.support == 0) {
    feat.context = Get1bpContext(align_ctx.ref_seq, info.ref_pos);
    feat.context_index = UnifiedVariantFeature::ContextIndex(feat.context);
  }

  IncrementFeature(feat, align_ctx, baseq_mean, min_distance, align_ctx.IsInsertionNearNonConcordantBase(info));
}

void ProcessDeletion(UnifiedVariantFeatures& features, const AlignContext& align_ctx, const AlignOpInfo& info) {
  if (info.read_pos == 0) {
    return;
  }

  const u8 baseq = bam_get_qual(align_ctx.bam_record)[info.read_pos - 1];
  if (baseq < align_ctx.params.min_bq) {
    return;
  }

  const auto min_distance = std::min(info.ref_pos - static_cast<u32>(align_ctx.bam_record->core.pos),
                                     align_ctx.ref_end - (info.ref_pos + info.op_len));
  if (static_cast<f32>(min_distance) < align_ctx.params.min_allowed_distance_from_end) {
    return;
  }

  const std::string ref_allele{align_ctx.ref_seq.substr(info.ref_pos - 1, info.op_len + 1)};
  if (IsAnyNotACTG(ref_allele)) {
    // Since we are not confident in variant calls at non-ACGT reference positions,
    // we do not extract variant features at this position.
    return;
  }

  const std::string alt_allele{align_ctx.ref_seq[info.ref_pos - 1]};
  const VariantId id(align_ctx.region.chrom, info.ref_pos - 1, ref_allele, alt_allele);

  UnifiedVariantFeature& feat = features[id];
  if (feat.support == 0) {
    feat.context = Get1bpContext(align_ctx.ref_seq, info.ref_pos);
    feat.context_index = UnifiedVariantFeature::ContextIndex(feat.context);
  }

  IncrementFeature(feat, align_ctx, baseq, min_distance, align_ctx.IsDeletionNearNonConcordantBase(info));
}

void ProcessMismatch(UnifiedVariantFeatures& features, const AlignContext& align_ctx, const AlignOpInfo& info) {
  const u8 baseq = bam_get_qual(align_ctx.bam_record)[info.read_pos];
  if (baseq < align_ctx.params.min_bq) {
    return;
  }

  const auto min_distance =
      std::min(info.ref_pos - static_cast<u32>(align_ctx.bam_record->core.pos), align_ctx.ref_end - info.ref_pos - 1);
  if (static_cast<f32>(min_distance) < align_ctx.params.min_allowed_distance_from_end) {
    return;
  }

  const char ref_base{align_ctx.ref_seq[info.ref_pos]};
  if (IsNotACTG(ref_base)) {
    // Since we are not confident in variant calls at non-ACGT reference positions,
    // we do not extract variant features at this position.
    return;
  }

  auto* begin_seq = bam_get_seq(align_ctx.bam_record);
  const char base = seq_nt16_str[bam_seqi(begin_seq, info.read_pos)];
  if (IsNotACTG(base)) {
    return;
  }
  const std::string alt_allele{base};
  const std::string ref_allele{ref_base};
  const VariantId id(align_ctx.region.chrom, info.ref_pos, ref_allele, alt_allele);

  UnifiedVariantFeature& feat = features[id];
  if (feat.support == 0) {
    feat.context = Get1bpContext(align_ctx.ref_seq, info.ref_pos);
    feat.context_index = UnifiedVariantFeature::ContextIndex(feat.context);
  }

  IncrementFeature(feat, align_ctx, baseq, min_distance, align_ctx.IsMismatchNearNonConcordantBase(info));
}

/**
 * @brief Update the match features in the UnifiedReferenceFeature.
 * @param feat The UnifiedReferenceFeature to update.
 * @param align_ctx The AlignContext containing alignment information.
 * @param mapq The mapping quality of the read.
 * @param baseq The base quality of the read at the position.
 * @param min_distance The minimum distance from the end of the read to the reference position.
 * @param weighted_depth The weighted depth of the read at this position.
 */
static void UpdateMatchFeatures(UnifiedReferenceFeature& feat,
                                const AlignContext& align_ctx,
                                const u8 mapq,
                                const u8 baseq,
                                const u64 min_distance,
                                const f64 weighted_depth) {
  if (feat.support == 0) {
    feat.mapq_min = mapq;
    feat.mapq_max = mapq;
    feat.distance_min = min_distance;
    feat.distance_max = min_distance;
  } else {
    feat.mapq_min = std::min(feat.mapq_min, mapq);
    feat.mapq_max = std::max(feat.mapq_max, mapq);
    feat.distance_min = std::min(feat.distance_min, min_distance);
    feat.distance_max = std::max(feat.distance_max, min_distance);
  }
  ++feat.support;
  feat.distance_sum += min_distance;
  feat.weighted_depth += weighted_depth;

#ifdef SOMATIC_ENABLE
  if (is_tumor_sample.has_value()) {
    if (is_tumor_sample.value()) {
      ++feat.tumor_support;
    } else {
      ++feat.normal_support;
    }
  }
#endif  //  SOMATIC_ENABLE

  feat.mapq_sum += mapq;
  if (mapq < kMapqCountThreshold60) {
    ++feat.mapq_lt60_count;
  }
  if (mapq < kMapqCountThreshold40) {
    ++feat.mapq_lt40_count;
  }
  if (mapq < kMapqCountThreshold30) {
    ++feat.mapq_lt30_count;
  }
  if (mapq < kMapqCountThreshold20) {
    ++feat.mapq_lt20_count;
  }

  feat.baseq_sum += baseq;
  if (baseq < kBaseqCountThreshold20) {
    ++feat.baseq_lt20_count;
  }

  const u32 familysize{align_ctx.plus_counts + align_ctx.minus_counts};
  feat.familysize_sum += familysize;
  if (familysize < kFamilySizeCountThreshold5) {
    ++feat.familysize_lt5_count;
    if (familysize < kFamilySizeCountThreshold3) {
      ++feat.familysize_lt3_count;
    }
  }
}

/**
 * @brief Update the non-homopolymer match features in the UnifiedReferenceFeature.
 * @param feat The UnifiedReferenceFeature to update.
 * @param mapq The mapping quality of the read.
 * @param baseq The base quality of the read at the position.
 * @param weighted_depth The weighted depth of the read at this position.
 */
static void UpdateMatchNonHomopolymerFeatures(UnifiedReferenceFeature& feat,
                                              const u8 mapq,
                                              const u8 baseq,
                                              const f64 weighted_depth) {
  if (feat.nonhomopolymer_support == 0) {
    feat.nonhomopolymer_mapq_min = mapq;
    feat.nonhomopolymer_mapq_max = mapq;
    feat.nonhomopolymer_baseq_min = baseq;
    feat.nonhomopolymer_baseq_max = baseq;
  } else {
    feat.nonhomopolymer_mapq_min = std::min(feat.nonhomopolymer_mapq_min, mapq);
    feat.nonhomopolymer_mapq_max = std::max(feat.nonhomopolymer_mapq_max, mapq);
    feat.nonhomopolymer_baseq_min = std::min(feat.nonhomopolymer_baseq_min, baseq);
    feat.nonhomopolymer_baseq_max = std::max(feat.nonhomopolymer_baseq_max, baseq);
  }
  ++feat.nonhomopolymer_support;
  feat.nonhomopolymer_baseq_sum += baseq;
  feat.nonhomopolymer_weighted_depth += weighted_depth;
  feat.nonhomopolymer_mapq_sum += mapq;
  if (mapq < kMapqCountThreshold60) {
    ++feat.nonhomopolymer_mapq_lt60_count;
  }
  if (mapq < kMapqCountThreshold40) {
    ++feat.nonhomopolymer_mapq_lt40_count;
  }
  if (mapq < kMapqCountThreshold30) {
    ++feat.nonhomopolymer_mapq_lt30_count;
  }
  if (mapq < kMapqCountThreshold20) {
    ++feat.nonhomopolymer_mapq_lt20_count;
  }
}

void ProcessMatch(UnifiedReferenceFeatures& features,
                  const AlignContext& align_ctx,
                  const u64 read_pos,
                  const u64 ref_pos) {
  const u8 baseq = bam_get_qual(align_ctx.bam_record)[read_pos];
  if (baseq < align_ctx.params.min_bq) {
    return;
  }

  const auto min_distance =
      std::min(ref_pos - static_cast<u32>(align_ctx.bam_record->core.pos), align_ctx.ref_end - ref_pos - 1);
  if (static_cast<f32>(min_distance) < align_ctx.params.min_allowed_distance_from_end) {
    return;
  }

  UnifiedReferenceFeature& feat = features[ref_pos];
  if (ContainsReadIdInsertIfNot(feat.read_ids, align_ctx.read_id)) {
    return;
  }

  const u8 mapq{align_ctx.bam_record->core.qual};

  if (align_ctx.read_type == ReadType::kSimplex) {
    ++feat.simplex;
    feat.mapq_sum_simplex += mapq;
    feat.distance_sum_simplex += min_distance;
  } else if (align_ctx.IsMatchNearNonConcordantBase(read_pos)) {
    if (align_ctx.plus_counts > 0 && align_ctx.minus_counts > 0) {
      feat.duplex_lowbq += 0.5;
    }
    feat.mapq_sum_lowbq += mapq;
    feat.distance_sum_lowbq += min_distance;
  } else {
    const f64 weighted_depth{static_cast<f64>(baseq) / kMaxPossibleBaseQual * mapq / kMaxPossibleMapQual};
    UpdateMatchFeatures(feat, align_ctx, mapq, baseq, min_distance, weighted_depth);

    if (!align_ctx.homopolymer_start.has_value() || align_ctx.homopolymer_start > ref_pos) {
      // this matching base is not in a homopolymer
      UpdateMatchNonHomopolymerFeatures(feat, mapq, baseq, weighted_depth);
    }
  }
}

u32 CountVariantsInRead(const vec<AlignOpInfo>& infos) {
  u32 num_var_ops{0};
  for (const auto& info : infos) {
    switch (info.op) {
      case kInsertion:
      case kDeletion:
      case kMismatch:
        ++num_var_ops;
        break;
      default:
        break;
    }
  }
  return num_var_ops;
}

u32 CountVariantsInRead(const vec<AlignOpInfo>& infos,
                        const bam1_t* record,
                        const std::string& ref_seq,
                        const std::string& chrom,
                        const StrUnorderedSet& skip_variants) {
  using enum AlignOp;
  u32 num_var_ops{0};
  for (const auto& info : infos) {
    switch (info.op) {
      case kInsertion: {
        const char ref_base{ref_seq[info.ref_pos]};
        const std::string ref_allele{ref_base};
        auto* read = bam_get_seq(record);
        std::string alt_allele(info.op_len + 1, 'N');
        alt_allele[0] = ref_base;
        for (u64 i = 0; i < info.op_len; ++i) {
          alt_allele[i + 1] = seq_nt16_str[bam_seqi(read, info.read_pos + i)];
        }
        if (!skip_variants.contains(GetVariantCorrelationKey(chrom, info.ref_pos, ref_allele, alt_allele, false))) {
          ++num_var_ops;
        }
        break;
      }
      case kDeletion: {
        const std::string ref_allele{ref_seq.substr(info.ref_pos - 1, info.op_len + 1)};
        const std::string alt_allele{ref_seq[info.ref_pos - 1]};
        if (!skip_variants.contains(GetVariantCorrelationKey(chrom, info.ref_pos, ref_allele, alt_allele, false))) {
          ++num_var_ops;
        }
        break;
      }
      case kMismatch: {
        for (u64 i = 0; i < info.op_len; ++i) {
          const char ref_base{ref_seq[info.ref_pos + i]};
          const std::string ref_allele{ref_base};
          auto* read = bam_get_seq(record);
          const char alt_base = seq_nt16_str[bam_seqi(read, info.read_pos + i)];
          const std::string alt_allele{alt_base};
          if (!skip_variants.contains(GetVariantCorrelationKey(chrom, info.ref_pos, ref_allele, alt_allele, false))) {
            ++num_var_ops;
          }
        }
        break;
      }
      default:
        break;
    }
  }
  return num_var_ops;
}

/**
 * @brief Find nearby non-ACGT bases in the sequence and mark them in the boolean vector.
 * @param result Vector indicating whether each base is at or nearby a non-ACGT base
 * @param seq Sequence string
 * @note result is a boolean vector of the same length as seq, where each position is set to true if it is at or
 *       nearby a non-ACGT base, and false otherwise.
 */
static void NearNonACGT(const std::string& seq, vec<bool>& result) {
  const u64 seq_len = seq.length();
  for (u64 i = 0; i < seq_len; ++i) {
    if (IsACTG(seq[i])) {
      continue;
    }
    if (i + 1 < seq_len) {
      result.at(i + 1) = true;
    }
    if (i > 0) {
      u64 prev{i - 1};
      if (IsNotACTG(seq[prev])) {
        // if the previous base is also non-ACGT, we skip it because it is already marked
        // as non-ACGT and we do not want to mark it again
        continue;
      }
      result.at(prev) = true;
      while (prev > 0 && seq[prev - 1] == seq[prev]) {
        // the homopolymer before 'N' are treated as if they were 'N'
        --prev;
        result.at(prev) = true;
      }
      if (prev > 0) {
        // set one more base upstream, which can be ambiguous because of the 'N'
        // e.g. "CANTG", both 'C' and 'A' can be ambiguous
        result.at(prev - 1) = true;
      }
    }
  }
}

vec<bool> FindNearbyNonConcordantBase(const bam1_t* alignment) {
  const auto seq_len = static_cast<u64>(alignment->core.l_qseq);
  vec<bool> near_n(seq_len, false);
  std::string seq(seq_len, 'N');
  // find N using the read sequence and base quality
  auto* const begin_seq = bam_get_seq(alignment);
  auto* const qual_data = bam_get_qual(alignment);
  for (u64 i = 0; i < seq_len; ++i) {
    const char c{seq_nt16_str[bam_seqi(begin_seq, i)]};
    if (IsNotACTG(c)) {
      near_n[i] = true;
    } else {
      const u8 qual = qual_data[i];
      if (qual >= kMinConcordantBaseQuality) {
        seq[i] = c;
      } else {
        near_n[i] = true;
      }
    }
  }
  NearNonACGT(seq, near_n);
  return near_n;
}

vec<bool> FindNearbyNonConcordantBase(const bam1_t* alignment, const vec<BaseType>& base_types) {
  const auto seq_len = static_cast<u64>(alignment->core.l_qseq);
  vec<bool> result(seq_len, false);
  if (base_types.empty()) {
    return result;
  }
  std::string seq(seq_len, 'N');
  // find N using the read sequence and base types
  auto* const seq_start = bam_get_seq(alignment);
  for (u64 i = 0; i < seq_len; ++i) {
    const char c{seq_nt16_str[bam_seqi(seq_start, i)]};
    if (IsNotACTG(c)) {
      result[i] = true;
    } else {
      if (base_types.at(i) == BaseType::kConcordant) {
        seq[i] = c;
      } else {
        result[i] = true;
      }
    }
  }
  NearNonACGT(seq, result);
  return result;
}

u64 GetAlignmentLength(const vec<AlignOpInfo>& infos) {
  u64 read_length{0};
  for (const auto& info : infos) {
    switch (info.op) {
      case kMatch:
      case kMismatch:
      case kInsertion:
        read_length += info.op_len;
        break;
      default:
        break;
    }
  }
  return read_length;
}

void ProcessAlignment(const AlignContext& align_ctx,
                      UnifiedVariantFeatures& var_feats,
                      UnifiedReferenceFeatures& ref_feats) {
  const auto num_infos{align_ctx.align_op_infos.size()};
  for (size_t a = 0; a < num_infos; ++a) {
    const auto& info = align_ctx.align_op_infos.at(a);
    switch (info.op) {
      case kInsertion: {
        if (align_ctx.CanProcessInsertion(info)) {
          ProcessInsertion(var_feats, align_ctx, info);
        }
        break;
      }
      case kDeletion: {
        // The deletion must overlap this region, but it cannot overlap the previous region.
        // If the deletion also overlaps the previous region, then it should be processed there instead of this region.
        // This ensures that a deletion spanning across regions is only processed in its left-most overlapping region,
        // therefore preventing additional incomplete features for the same deletion.
        if (align_ctx.CanProcessDeletion(info)) {
          ProcessDeletion(var_feats, align_ctx, info);
        }
        break;
      }
      case kMismatch: {
        if (align_ctx.CanProcessMismatch(info)) {
          ProcessMismatch(var_feats, align_ctx, info);
        }
        break;
      }
      case kMatch: {
        // We skip reference feature extraction at a position if the following base is an insertion because
        // this position is part of the ALT allele, i.e. the insertion occurs at this position.
        const u64 end{a + 1 < num_infos && align_ctx.align_op_infos[a + 1].op == kInsertion ? info.op_len - 1
                                                                                            : info.op_len};
        for (u64 i = 0; i < end; ++i) {
          // Do not skip positions before `region_start` because these positions can be reference features for
          // deletions overlapping `region_start`.
          // Use `<=` instead of `<` to extract reference features for deletions that start at `region_end - 1`.
          const auto ref_pos = info.ref_pos + i;
          const auto read_pos = info.read_pos + i;
          if (align_ctx.CanProcessMatch(read_pos, ref_pos)) {
            ProcessMatch(ref_feats, align_ctx, read_pos, ref_pos);
          }
        }
        break;
      }
      default:
        // ignore all other operations
        break;
    }
  }
}

}  // namespace xoos::svc
