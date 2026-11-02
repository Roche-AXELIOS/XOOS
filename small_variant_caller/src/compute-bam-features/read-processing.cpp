#include "read-processing.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <numeric>
#include <optional>
#include <ranges>
#include <span>

#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/util/math.h>

#include "core/config.h"
#include "core/filtering.h"
#include "core/genotype.h"
#include "core/variant-info.h"
#include "read-name-utils.h"
#include "util/log-util.h"
#include "util/region-util.h"
#include "util/seq-util.h"

namespace xoos::svc {

using enum AlignOp;
using yc_decode::BaseType;

static constexpr f64 kAfqGtDistThreshold = 0.2;
static constexpr f64 kAfqGtMinDuplexDp = 4;
// Maximum possible Manhattan distance between two 3-component allele-frequency vectors whose components each sum
// to 1; used to convert a distance into a similarity score (`afq_gt* = kAfqGtMaxDistance - distance`).
static constexpr f64 kAfqGtMaxDistance = 2.0;

namespace {
// Weight applied to duplex reads when computing the weighted score; duplex reads consist of a plus and a
// minus read for v9.2 chemistry, so they count double relative to non-duplex reads.
constexpr f64 kDuplexWeight{2.0};
// Fractional weight applied to duplex reads when computing strand bias.
constexpr f64 kDuplexStrandWeight{0.5};
// Increment applied to the low-base-quality duplex feature when a low-quality duplex base is observed.
constexpr f64 kDuplexLowbqIncrement{0.5};
// Baseline alignment-bias value when reverse-strand support is balanced.
constexpr f64 kAlignmentBiasBaseline{0.25};

// Indel-length classes (bp) for per-length support buckets; longer indels fall into the ">3bp" bucket.
constexpr s64 kIndelLen1bp{1};
constexpr s64 kIndelLen2bp{2};
constexpr s64 kIndelLen3bp{3};

/**
 * @brief Determine if the provided read_id has been seen before, if not add it to the list of seen read_ids.
 * @return true if read_ids contains read_id, false otherwise and read_id is added to read_ids
 */
bool ContainsReadIdInsertIfNot(ReadIds& read_ids, const ReadId read_id) {
  // Read ids are inserted at the back of the vector, therefore we iterate in reverse
  // as we are more likely to see recent read ids earlier.
  if (std::ranges::find(std::ranges::reverse_view(read_ids), read_id) != read_ids.crend()) {
    return true;
  }
  read_ids.emplace_back(read_id);
  return false;
}

// Accumulate the per-variant duplex depth for an indel variant into length-class buckets on the reference feature.
// @param duplex_dp Per-variant duplex depth (integer duplex read count), matching the duplex_dp definition so the
//        resulting allele fractions are proper fractions.
// @pre vid.type must be kDeletion or kInsertion.
void AccumulateIndelClassSupport(UnifiedReferenceFeature& ref_feat, const VariantId& vid, const u32 duplex_dp) {
  const auto ref_len = static_cast<s64>(vid.ref.size());
  const auto alt_len = static_cast<s64>(vid.alt.size());
  const auto indel_len = std::abs(alt_len - ref_len);
  if (vid.type == VariantType::kDeletion) {
    if (indel_len == kIndelLen1bp) {
      ref_feat.del_1bp_support += duplex_dp;
    } else if (indel_len == kIndelLen2bp) {
      ref_feat.del_2bp_support += duplex_dp;
    } else if (indel_len == kIndelLen3bp) {
      ref_feat.del_3bp_support += duplex_dp;
    } else {
      ref_feat.del_gt3bp_support += duplex_dp;
    }
  } else {
    if (indel_len == kIndelLen1bp) {
      ref_feat.ins_1bp_support += duplex_dp;
    } else if (indel_len == kIndelLen2bp) {
      ref_feat.ins_2bp_support += duplex_dp;
    } else if (indel_len == kIndelLen3bp) {
      ref_feat.ins_3bp_support += duplex_dp;
    } else {
      ref_feat.ins_gt3bp_support += duplex_dp;
    }
  }
}

/**
 * @brief Compute the allele-frequency-based genotype-quality features for a variant.
 *
 * Scores how closely the duplex allele-frequency vector (variant AF, other-allele AF, reference AF) matches the ideal
 * vectors for the 0/1, 1/1, and 1/2 genotypes, and calls the closest genotype when the position has enough duplex depth
 * and one genotype is a clear, unique best match.
 *
 * @param var_feat Variant feature to populate (afq_gt01, afq_gt11, afq_gt12, duplex_gq, afq_gt).
 * @param var_af Variant-allele duplex allele frequency at the position.
 * @param var_other_af Other-allele (non-variant, non-reference) duplex allele frequency at the position.
 * @param ref_af Reference-allele duplex allele frequency at the position.
 * @param duplex_dp Total duplex depth at the position.
 */
// The observed AF components are passed as full-precision f64 (not read back from the now-f32 struct fields) so that
// genotype assignment is unaffected by the storage narrowing of `duplex_af`/`duplex_other_af` at exact thresholds.
void ComputeAfqGenotypeFeatures(
    UnifiedVariantFeature& var_feat, const f64 var_af, const f64 var_other_af, const f64 ref_af, const f64 duplex_dp) {
  // Number of components in an allele-frequency vector: variant AF, other-allele AF, and reference AF.
  constexpr size_t kAfVectorSize{3};
  using AfVector = std::array<f64, kAfVectorSize>;
  // Manhattan distance: sum of absolute per-component differences.
  const auto manhattan = [](const AfVector& a, const AfVector& b) {
    return std::inner_product(
        a.begin(), a.end(), b.begin(), 0.0, std::plus{}, [](const f64 x, const f64 y) { return std::abs(x - y); });
  };
  // Observed duplex AF vector (variant AF, other-allele AF, reference AF) and the ideal vectors for each genotype.
  const AfVector observed = {var_af, var_other_af, ref_af};
  const AfVector ideal_gt01 = {0.5, 0.0, 0.5};
  const AfVector ideal_gt11 = {1.0, 0.0, 0.0};
  const AfVector ideal_gt12 = {0.5, 0.5, 0.0};
  const f64 dist_gt01 = manhattan(observed, ideal_gt01);
  const f64 dist_gt11 = manhattan(observed, ideal_gt11);
  const f64 dist_gt12 = manhattan(observed, ideal_gt12);
  var_feat.afq_gt01 = kAfqGtMaxDistance - dist_gt01;
  var_feat.afq_gt11 = kAfqGtMaxDistance - dist_gt11;
  var_feat.afq_gt12 = kAfqGtMaxDistance - dist_gt12;
  const f64 min_dist = std::min({dist_gt01, dist_gt11, dist_gt12});
  var_feat.duplex_gq = (1.0 - min_dist) * duplex_dp;
  var_feat.afq_gt = Genotype::kGT00;
  if (duplex_dp < kAfqGtMinDuplexDp || min_dist > kAfqGtDistThreshold) {
    return;
  }
  const std::array<bool, 3> is_min = {math::IsEqualWithTolerance(dist_gt01, min_dist),
                                      math::IsEqualWithTolerance(dist_gt11, min_dist),
                                      math::IsEqualWithTolerance(dist_gt12, min_dist)};
  // Only call a genotype when exactly one candidate is the unique closest match.
  if (std::ranges::count(is_min, true) != 1) {
    return;
  }
  // NOSONAR(cpp:M23_002) — not always true: the `std::ranges::count(is_min, true) != 1` early return above
  // guarantees exactly one is_min element is true here, so this reaches the else branches for the kGT11/kGT12 cases.
  if (is_min[0]) {  // NOSONAR(cpp:M23_002)
    var_feat.afq_gt = Genotype::kGT01;
  } else if (is_min[1]) {
    var_feat.afq_gt = Genotype::kGT11;
  } else {
    var_feat.afq_gt = Genotype::kGT12;
  }
}

/**
 * @brief Compute allele-frequency-derived features for a single variant relative to its reference position.
 *
 * Populates the MAPQ, BASEQ, and duplex allele frequencies plus the duplex-other and afq genotype features, using the
 * position-level totals accumulated across all alleles.
 */
void UpdateVariantAlleleFrequencies(UnifiedVariantFeature& var_feat,
                                    UnifiedReferenceFeature& ref_feat,
                                    const u32 total_mapq_sum,
                                    const f64 total_baseq_sum,
                                    const bool include_duplex_lowbq) {
  var_feat.mq_af = 0;
  ref_feat.mq_af = 0;
  if (total_mapq_sum > 0) {
    const auto total = static_cast<f64>(total_mapq_sum);
    var_feat.mq_af = static_cast<f32>(var_feat.mapq_sum / total);
    ref_feat.mq_af = static_cast<f32>(ref_feat.nonhomopolymer_mapq_sum / total);
  }
  var_feat.bq_af = 0;
  ref_feat.bq_af = 0;
  if (total_baseq_sum > 0) {
    var_feat.bq_af = static_cast<f32>(var_feat.baseq_sum / total_baseq_sum);
    ref_feat.bq_af = static_cast<f32>(ref_feat.nonhomopolymer_baseq_sum / total_baseq_sum);
  }
  const f64 var_duplex = (include_duplex_lowbq ? var_feat.duplex_lowbq : 0) + var_feat.duplex;
  const f64 ref_duplex = (include_duplex_lowbq ? ref_feat.duplex_lowbq : 0) + ref_feat.support;
  // Compute duplex AFs in full precision; store as f32 but keep the f64 values for genotype assignment.
  f64 var_duplex_af = 0;
  f64 ref_duplex_af = 0;
  if (ref_feat.duplex_dp > 0) {
    var_duplex_af = var_duplex / ref_feat.duplex_dp;
    ref_duplex_af = ref_duplex / ref_feat.duplex_dp;
  }
  var_feat.duplex_af = static_cast<f32>(var_duplex_af);
  ref_feat.duplex_af = static_cast<f32>(ref_duplex_af);
  var_feat.duplex_other = std::max(0.0, ref_feat.duplex_dp - var_duplex - ref_duplex);
  const f64 var_duplex_other_af = std::max(0.0, 1.0 - var_duplex_af - ref_duplex_af);
  var_feat.duplex_other_af = static_cast<f32>(var_duplex_other_af);
  ComputeAfqGenotypeFeatures(var_feat, var_duplex_af, var_duplex_other_af, ref_duplex_af, ref_feat.duplex_dp);
}
}  // namespace

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
    feature.mapq_lt60_ratio = static_cast<f32>(feature.mapq_lt60_count / support);
    feature.mapq_lt40_ratio = static_cast<f32>(feature.mapq_lt40_count / support);
    feature.mapq_lt30_ratio = static_cast<f32>(feature.mapq_lt30_count / support);
    feature.mapq_lt20_ratio = static_cast<f32>(feature.mapq_lt20_count / support);
    feature.familysize_lt5_ratio = static_cast<f32>(feature.familysize_lt5_count / support);
    feature.familysize_lt3_ratio = static_cast<f32>(feature.familysize_lt3_count / support);
    feature.baseq_lt20_ratio = static_cast<f32>(feature.baseq_lt20_count / support);
    feature.distance_mean = static_cast<f64>(feature.distance_sum) / support;
    feature.baseq_mean = static_cast<f32>(feature.baseq_sum / support);
    feature.mapq_mean = static_cast<f32>(feature.mapq_sum / support);
    feature.familysize_mean = static_cast<f32>(feature.familysize_sum / support);
    const f64 f = static_cast<f64>(feature.support_reverse) / support;
    feature.alignmentbias = kAlignmentBiasBaseline - f * (1 - f);
  }
  // for v9.2 chemistry, `duplex` consists of plus and minus reads. Thus, `duplex` is 2x whereas `nonduplex` is 1x
  feature.weighted_score = feature.nonduplex + (kDuplexWeight * feature.duplex);
  const auto tmp = (feature.duplex * kDuplexStrandWeight + feature.plusonly) /
                   (feature.duplex + feature.plusonly + feature.minusonly);
  feature.strandbias = tmp * (1 - tmp);

  if (feature.duplex_lowbq > 0) {
    // `duplex_lowbq` is incremented by 0.5 for each `lowbq` read, so multiply it by 2 to get the actual count
    const auto duplex_lowbq = 2 * feature.duplex_lowbq;
    feature.mapq_mean_lowbq = static_cast<f32>(feature.mapq_sum_lowbq / duplex_lowbq);
    feature.distance_mean_lowbq = static_cast<f64>(feature.distance_sum_lowbq) / duplex_lowbq;
  }
  if (feature.simplex > 0) {
    const auto simplex = static_cast<f64>(feature.simplex);
    feature.mapq_mean_simplex = static_cast<f32>(feature.mapq_sum_simplex / simplex);
    feature.distance_mean_simplex = static_cast<f64>(feature.distance_sum_simplex) / simplex;
  }
}

void UpdateDerivedFeatureValues(UnifiedReferenceFeature& feature) {
  if (feature.support > 0) {
    const auto support = static_cast<f64>(feature.support);
    feature.mapq_lt60_ratio = static_cast<f32>(feature.mapq_lt60_count / support);
    feature.mapq_lt40_ratio = static_cast<f32>(feature.mapq_lt40_count / support);
    feature.mapq_lt30_ratio = static_cast<f32>(feature.mapq_lt30_count / support);
    feature.mapq_lt20_ratio = static_cast<f32>(feature.mapq_lt20_count / support);
    feature.mapq_mean = static_cast<f32>(feature.mapq_sum / support);
    feature.baseq_mean = static_cast<f32>(feature.baseq_sum / support);
    feature.baseq_lt20_ratio = static_cast<f32>(feature.baseq_lt20_count / support);
    feature.familysize_mean = static_cast<f32>(feature.familysize_sum / support);
    feature.familysize_lt5_ratio = static_cast<f32>(feature.familysize_lt5_count / support);
    feature.familysize_lt3_ratio = static_cast<f32>(feature.familysize_lt3_count / support);
    feature.distance_mean = static_cast<f64>(feature.distance_sum) / support;
    feature.nonhomopolymer_mapq_lt60_ratio = static_cast<f32>(feature.nonhomopolymer_mapq_lt60_count / support);
    feature.nonhomopolymer_mapq_lt40_ratio = static_cast<f32>(feature.nonhomopolymer_mapq_lt40_count / support);
    feature.nonhomopolymer_mapq_lt30_ratio = static_cast<f32>(feature.nonhomopolymer_mapq_lt30_count / support);
    feature.nonhomopolymer_mapq_lt20_ratio = static_cast<f32>(feature.nonhomopolymer_mapq_lt20_count / support);
    feature.nonhomopolymer_baseq_mean = static_cast<f32>(feature.nonhomopolymer_baseq_sum / support);
    feature.nonhomopolymer_mapq_mean = static_cast<f32>(feature.nonhomopolymer_mapq_sum / support);
  }
  if (feature.duplex_lowbq > 0) {
    const auto duplex_lowbq = 2 * feature.duplex_lowbq;
    feature.mapq_mean_lowbq = static_cast<f32>(feature.mapq_sum_lowbq / duplex_lowbq);
    feature.distance_mean_lowbq = static_cast<f64>(feature.distance_sum_lowbq) / duplex_lowbq;
  }
  if (feature.simplex > 0) {
    const auto simplex = static_cast<f64>(feature.simplex);
    feature.mapq_mean_simplex = static_cast<f32>(feature.mapq_sum_simplex / simplex);
    feature.distance_mean_simplex = static_cast<f64>(feature.distance_sum_simplex) / simplex;
  }
}

void UpdateDerivedFeatureValues(VarIdToVarBamFeatures& variant_features,
                                PosToRefBamFeatures& reference_features,
                                const vec<VariantId>& vids,
                                const bool include_duplex_lowbq) {
  auto& ref_feat = reference_features[vids[0].GetRefFeaturePos()];
  UpdateDerivedFeatureValues(ref_feat);

  // Tally total MAPQ sum, total BASEQ sum, duplex DP, and number of variants at this position.
  u32 total_mapq_sum = ref_feat.nonhomopolymer_mapq_sum;
  f64 total_baseq_sum = ref_feat.nonhomopolymer_baseq_sum;
  ref_feat.duplex_dp = (include_duplex_lowbq ? ref_feat.duplex_lowbq : 0) + ref_feat.support;
  ref_feat.num_alt = 0;
  for (const auto& vid : vids) {
    const auto itr = variant_features.find(vid);
    if (itr != variant_features.end()) {
      ++ref_feat.num_alt;
      const auto& var_feat = itr->second;
      total_mapq_sum += var_feat.mapq_sum;
      total_baseq_sum += var_feat.baseq_sum;
      // Compute the total duplex dp here over all variant and ref positions with regular and lowbq
      ref_feat.duplex_dp += (include_duplex_lowbq ? var_feat.duplex_lowbq : 0) + var_feat.duplex;
    }
  }
  // Note: per-position aggregate features (indel/SNV support, AF, and log-ratios) are computed separately in
  // UpdateAggregateFeatures, keyed by the raw variant position so that co-located SNVs, insertions, and deletions
  // share the same values.
  // Update features for each variant and calculate allele frequencies and other values.
  for (const auto& vid : vids) {
    const auto itr = variant_features.find(vid);
    if (itr == variant_features.end()) {
      continue;
    }

    auto& var_feat = itr->second;
    UpdateDerivedFeatureValues(var_feat);

    // Update features that depend on both variant-specific and reference position specific values and totals.
    UpdateVariantAlleleFrequencies(var_feat, ref_feat, total_mapq_sum, total_baseq_sum, include_duplex_lowbq);
    const u32 duplex_sum = var_feat.duplex + ref_feat.nonhomopolymer_support;
    if (duplex_sum > 0) {
      const auto duplex = static_cast<f64>(var_feat.duplex);
      var_feat.adt = std::abs(duplex - ref_feat.nonhomopolymer_support) / duplex_sum;
      var_feat.adtl = std::log10(duplex_sum) * var_feat.adt;
      var_feat.indel_af = static_cast<f32>(duplex / duplex_sum);
    }
  }
}

void UpdateAggregateFeatures(VarIdToVarBamFeatures& variant_features,
                             PosToRefBamFeatures& reference_features,
                             const vec<VariantId>& vids) {
  // These per-position aggregate features are keyed by the raw variant position (GetAggFeaturePos), not the
  // reference-feature position (GetRefFeaturePos, which offsets deletions by 1). Keying by the raw position means all
  // variants reported at the same locus — SNVs, insertions, and deletions — share the same aggregate values.
  auto& ref_feat = reference_features[vids[0].GetAggFeaturePos()];

  // Duplex depth denominator: total duplex depth at this raw position = reference-matching duplex support plus the
  // summed duplex support of every non-deletion alt allele reported here. Like the original ref_feat.duplex_dp, the
  // reference term is ref_feat.support. These multi-allelic aggregate features intentionally use the duplex count only
  // and never the duplex_lowbq reads, regardless of the --include-duplex-lowbq CLI option (unlike ref_feat.duplex_dp /
  // duplex_af, which the option still governs), so no ref_feat.duplex_lowbq term is added.
  //
  // Deletion reads carry a CIGAR match at the anchor base P (the deletion itself is at P+1), so they are already
  // counted in ref_feat.support@P. Because these aggregate features are keyed to the raw position P, a deletion read
  // would otherwise be counted twice in the denominator: once as reference support and once as its own alt duplex
  // support. Insertion reads do not have this issue — the match op preceding the insertion is skipped, so insertion
  // reads never increment ref_feat.support@P. To avoid double-counting, deletion alt support is not added to the
  // denominator; deletion reads contribute their single count via the reference term already in ref_feat.support.
  // Accumulates integer duplex read counts (ref_feat.support and per-variant duplex counts are both u32), so u32.
  u32 total_duplex_dp = ref_feat.support;
  for (const auto& vid : vids) {
    const auto itr = variant_features.find(vid);
    if (itr == variant_features.end()) {
      continue;
    }
    const auto& var_feat = itr->second;
    // Per-variant duplex depth. Feeds the per-class support numerator (always) and the summed denominator. Uses the
    // duplex count only, never the duplex_lowbq reads (see denominator note above). Integer read count (u32).
    const u32 var_duplex_dp = var_feat.duplex;
    // Deletion reads are already inside ref_feat.support (anchor match at P), so adding them here too would
    // double-count them in the denominator. Add only non-deletion alt duplex support; deletion reads keep their
    // single contribution via the reference term. Deletions still count toward their per-class support numerator
    // below. (Insertion/SNV reads are not in ref_feat.support at P, so they must be added here.)
    if (vid.type != VariantType::kDeletion) {
      total_duplex_dp += var_duplex_dp;
    }
    // Aggregate duplex support across all alt alleles at this position, bucketed by indel class.
    if (vid.type == VariantType::kDeletion || vid.type == VariantType::kInsertion) {
      AccumulateIndelClassSupport(ref_feat, vid, var_duplex_dp);
    } else if (vid.type == VariantType::kSNV) {
      // Aggregate duplex support across all SNV alt alleles at this position. Models SNVs co-occurring at the same
      // locus as an indel (and equals the SNV's own support when the variant itself is a lone SNV at the position).
      ref_feat.snv_support += var_duplex_dp;
    }
  }
  // Per-position allele frequencies, normalized by total duplex depth at this position. dp is f64, so the u32 support
  // numerators promote to floating point before division; results are narrowed to the f32 af fields.
  if (total_duplex_dp > 0) {
    const auto dp = static_cast<f64>(total_duplex_dp);
    ref_feat.del_1bp_af = static_cast<f32>(ref_feat.del_1bp_support / dp);
    ref_feat.del_2bp_af = static_cast<f32>(ref_feat.del_2bp_support / dp);
    ref_feat.del_3bp_af = static_cast<f32>(ref_feat.del_3bp_support / dp);
    ref_feat.del_gt3bp_af = static_cast<f32>(ref_feat.del_gt3bp_support / dp);
    ref_feat.ins_1bp_af = static_cast<f32>(ref_feat.ins_1bp_support / dp);
    ref_feat.ins_2bp_af = static_cast<f32>(ref_feat.ins_2bp_support / dp);
    ref_feat.ins_3bp_af = static_cast<f32>(ref_feat.ins_3bp_support / dp);
    ref_feat.ins_gt3bp_af = static_cast<f32>(ref_feat.ins_gt3bp_support / dp);
    ref_feat.snv_af = static_cast<f32>(ref_feat.snv_support / dp);
  }
  // Natural-log ratios of indel-class duplex depths. Adding 1 to both numerator and denominator avoids
  // division by zero and log(0). The support fields are u32, so the +1.0 literal forces floating-point division
  // (integer division here would truncate the ratio to 0 or 1 and corrupt the feature); results narrow to f32.
  ref_feat.del_1bp_ins_1bp_af_logratio =
      static_cast<f32>(std::log((ref_feat.del_1bp_support + 1.0) / (ref_feat.ins_1bp_support + 1.0)));
  ref_feat.del_1bp_2bp_af_logratio =
      static_cast<f32>(std::log((ref_feat.del_1bp_support + 1.0) / (ref_feat.del_2bp_support + 1.0)));
  ref_feat.del_2bp_3bp_af_logratio =
      static_cast<f32>(std::log((ref_feat.del_2bp_support + 1.0) / (ref_feat.del_3bp_support + 1.0)));
  ref_feat.del_3bp_gt3bp_af_logratio =
      static_cast<f32>(std::log((ref_feat.del_3bp_support + 1.0) / (ref_feat.del_gt3bp_support + 1.0)));
  ref_feat.ins_1bp_2bp_af_logratio =
      static_cast<f32>(std::log((ref_feat.ins_1bp_support + 1.0) / (ref_feat.ins_2bp_support + 1.0)));
  ref_feat.ins_2bp_3bp_af_logratio =
      static_cast<f32>(std::log((ref_feat.ins_2bp_support + 1.0) / (ref_feat.ins_3bp_support + 1.0)));
  ref_feat.ins_3bp_gt3bp_af_logratio =
      static_cast<f32>(std::log((ref_feat.ins_3bp_support + 1.0) / (ref_feat.ins_gt3bp_support + 1.0)));
}

void UpdateTumorNormalFeatures(BamRegionFeatureCollection& bam_feats,
                               const vec<VariantId>& vids,
                               const bool include_duplex_lowbq) {
  // Update derived feature values for tumor sample's variant and reference features
  UpdateDerivedFeatureValues(bam_feats.tumor_var_features, bam_feats.tumor_ref_features, vids, include_duplex_lowbq);

  // Update derived feature values for normal sample's variant and reference features
  UpdateDerivedFeatureValues(bam_feats.normal_var_features, bam_feats.normal_ref_features, vids, include_duplex_lowbq);

  // Calculate tumor-normal AF ratio for each variant at this position
  for (const auto& vid : vids) {
    const auto tumor_itr = bam_feats.tumor_var_features.find(vid);
    if (tumor_itr == bam_feats.tumor_var_features.end()) {
      continue;
    }
    const f64 tumor_af = tumor_itr->second.duplex_af;
    f64 af_sum = tumor_af;

    bool has_normal = false;
    const auto normal_itr = bam_feats.normal_var_features.find(vid);
    if (normal_itr != bam_feats.normal_var_features.end()) {
      has_normal = true;
      af_sum += normal_itr->second.duplex_af;
    }

    if (af_sum > 0) {
      const auto ratio = static_cast<f32>(tumor_af / af_sum);
      // combined, tumor, and normal features share the same value
      // it can be serialized as either "bam_rat", "tumor_bam_rat", or "normal_bam_rat"
      tumor_itr->second.tn_af_ratio = ratio;
      if (has_normal) {
        normal_itr->second.tn_af_ratio = ratio;
      }
      const auto itr = bam_feats.var_features.find(vid);
      if (itr != bam_feats.var_features.end()) {
        itr->second.tn_af_ratio = ratio;
      }
    }
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

namespace {
/**
 * @brief Increment the duplex base-type counter for a variant position.
 *
 * When the YC tag is available, the per-position BaseType determines which counter
 * (concordant / simplex / discordant) is incremented. When the tag is absent (legacy
 * mode), duplex reads default to concordant.
 */
template <typename Feature>
void IncrementDuplexBaseType(Feature& feature, const AlignContext& align_ctx, const u64 read_pos) {
  using enum BaseType;
  if (!align_ctx.base_types.empty() && read_pos < align_ctx.base_types.size()) {
    switch (align_ctx.base_types.at(read_pos)) {
      case kConcordant:
        ++feature.duplex_concordant;
        break;
      case kSimplex:
        ++feature.duplex_simplex;
        break;
      case kDiscordant:
        ++feature.duplex_discordant;
        break;
      default:
        break;
    }
  } else if (align_ctx.read_type == ReadType::kDuplex) {
    ++feature.duplex_concordant;
  } else {
    // Non-duplex read without base_types — no base-type counter to increment
  }
}

/**
 * @brief Increment the variant support counters for the read's sequencing class.
 * @param feature Variant feature counters to update.
 * @param align_ctx Alignment context containing the sequencing protocol and strand counts.
 */
void IncrementVariantReadClassCounters(UnifiedVariantFeature& feature, const AlignContext& align_ctx) {
  if (align_ctx.params->sequencing_protocol == SequencingProtocol::kSimplex) {
    // SBX-S reads have no plus/minus family counts and are non-duplex support.
    ++feature.nonduplex;
  } else if (align_ctx.plus_counts > 0 && align_ctx.minus_counts > 0) {
    // Support from both strands identifies a duplex read.
    ++feature.duplex;
  } else {
    // Support from at most one strand identifies a non-duplex read.
    ++feature.nonduplex;
    if (align_ctx.plus_counts > 0) {
      // Record support from only the plus strand.
      ++feature.plusonly;
    }
    if (align_ctx.minus_counts > 0) {
      // Record support from only the minus strand.
      ++feature.minusonly;
    }
  }
}

}  // namespace

void IncrementFeature(UnifiedVariantFeature& feature,
                      const AlignContext& align_ctx,
                      const f64 baseq,
                      const u64 distance,
                      const bool near_non_concordant,
                      const u64 read_pos) {
  if (ContainsReadIdInsertIfNot(feature.read_ids, align_ctx.read_id)) {
    // Skip previously processed reads to avoid double counting
    return;
  }

  const u8 mapq{align_ctx.bam_record->core.qual};
  const bool is_simplex_protocol{align_ctx.params->sequencing_protocol == SequencingProtocol::kSimplex};
  if (align_ctx.read_type == ReadType::kSimplex) {
    ++feature.simplex;
    feature.mapq_sum_simplex += mapq;
    feature.distance_sum_simplex += distance;
    if (!is_simplex_protocol) {
      return;
    }
  }

  if (!is_simplex_protocol && near_non_concordant) {
    if (align_ctx.plus_counts > 0 && align_ctx.minus_counts > 0) {
      feature.duplex_lowbq += kDuplexLowbqIncrement;
    }
    feature.mapq_sum_lowbq += mapq;
    feature.distance_sum_lowbq += distance;
    return;
  }

  if (!is_simplex_protocol) {
    IncrementDuplexBaseType(feature, align_ctx, read_pos);
  }

  IncrementVariantReadClassCounters(feature, align_ctx);

  const u32 familysize{is_simplex_protocol ? kSimplexReadFamilySize : align_ctx.plus_counts + align_ctx.minus_counts};
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
  const bool is_reverse = bam_is_rev(align_ctx.bam_record);
  if (is_reverse) {
    ++feature.support_reverse;
  }

  feature.weighted_depth += baseq / kMaxPossibleBaseQual * mapq / kMaxPossibleMapQual;
}

vec<AlignOpInfo> GetAlignOpInfos(const bam1_t* const record, const std::string& target_seq) {
  vec<AlignOpInfo> align_infos;
  align_infos.reserve(record->core.n_cigar);

  // The 0-based position in the query sequence
  u64 q_pos = 0;
  // The 0-based position in the target sequence
  auto t_pos = static_cast<u64>(record->core.pos);

  const std::span<const u32> cigar{bam_get_cigar(record), record->core.n_cigar};
  const u8* const seq = bam_get_seq(record);

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
                           const bam1_t* const bam_record,
                           const Region& region,
                           const std::string& ref_seq,
                           const ReadId read_id,
                           const std::set<u64>& vcf_feat_positions)
    : params(&params),
      align_op_infos(&align_op_infos),
      bam_record(bam_record),
      region(&region),
      ref_seq(&ref_seq),
      read_id(read_id),
      vcf_feat_positions(&vcf_feat_positions),
      has_vcf_feats(!vcf_feat_positions.empty()) {
  if (!InitFamilySize()) {
    return;
  }
  InitReadType();
  InitRefContext();
  ClassifyTumorNormal();
}

bool AlignContext::InitFamilySize() {
  u32 family_size = 0;
  if (params->sequencing_protocol == SequencingProtocol::kSimplex) {
    family_size = kSimplexReadFamilySize;
  } else if (IsDuplexProtocol(params->sequencing_protocol)) {
    plus_counts = kDuplexPlusCounts;
    minus_counts = kDuplexMinusCounts;
    family_size = kDuplexReadFamilySize;
  } else {
    ParseReadName(bam_get_qname(bam_record), plus_counts, minus_counts, family_size);
  }
  if (family_size < params->min_family_size) {
    skip_read = true;
    return false;
  }
  return true;
}

void AlignContext::InitReadType() {
  using enum ReadType;
  if (params->sequencing_protocol == SequencingProtocol::kSimplex) {
    read_type = kSimplex;
    return;
  }
  const bool is_duplex = IsDuplexProtocol(params->sequencing_protocol);
  if (!is_duplex) {
    read_type = kUmiConsensus;
    return;
  }

  read_type = kDuplex;
  if (params->decode_yc == YcDecodeMethod::kNone) {
    near_non_concordant = FindNearbyNonConcordantBase(bam_record);
  } else {
    // After YC split, the record is one of the two constituent reads.
    // The read type stays "Duplex" for updating relevant duplex features.
    const auto& yc_tag = yc_decode::DeserializeYcTag(bam_record);
    base_types = yc_tag.GetBaseTypes();
    if (base_types.empty() && params->sequencing_protocol == SequencingProtocol::kDuplexSimplex) {
      // Absence of YC tag indicates a simplex read in DuplexSimplex protocol
      read_type = kSimplex;
    } else {
      near_non_concordant = FindNearbyNonConcordantBase(bam_record, base_types);
    }
  }

  if (!base_types.empty() && params->min_base_type > BaseType::kDiscordant) {
    min_base_type = params->min_base_type;
  }
}

void AlignContext::InitRefContext() {
  ref_end = align_op_infos->back().ref_pos + align_op_infos->back().op_len;
  // find the homopolymer overlapping the last aligned reference base position, which is `ref_end - 1`
  homopolymer_start = FindOverlappingHomopolymer(
      *ref_seq, ref_end - 1, params->min_homopolymer_length, align_op_infos->front().ref_pos);
}

void AlignContext::ClassifyTumorNormal() {
  if (params->tumor_rg_ids.has_value() && !params->tumor_rg_ids->empty()) {
    const auto rg_id = io::BamAuxGet<std::string>(bam_record, "RG");
    is_tumor_sample = rg_id.has_value() && params->tumor_rg_ids->contains(rg_id.value());
  } else if (params->tumor_sample_name.has_value()) {
    const auto rg_id = io::BamAuxGet<std::string>(bam_record, "RG");
    is_tumor_sample = rg_id.has_value() && (params->tumor_sample_name.value() == rg_id.value());
  } else {
    // No tumor read-group or sample-name configured: leave the existing classification unchanged.
  }
}

bool AlignContext::CanProcessInsertion(const AlignOpInfo& info) const {
  if (skip_read) {
    return false;
  }
  if (info.ref_pos < region->start || info.ref_pos >= region->end) {
    // insertion out of bounds
    return false;
  }
  if (params->filter_homopolymer == HomopolymerFilter::kAlignmentEnd && homopolymer_start.has_value() &&
      homopolymer_start.value() < info.ref_pos) {
    // insertion is in a homopolymer region
    return false;
  }
  if (has_vcf_feats && !vcf_feat_positions->contains(info.ref_pos)) {
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
  if (!IntervalOverlap(region->start, region->end, info.ref_pos, info.ref_pos + info.op_len)) {
    // deletion does not overlap with the region
    return false;
  }
  if (region->prev_interval.has_value() &&
      IntervalOverlap(region->prev_interval.value(), info.ref_pos, info.ref_pos + info.op_len)) {
    // deletion overlaps with the previous interval
    // it should be processed in that region instead of this region
    return false;
  }
  if (params->filter_homopolymer == HomopolymerFilter::kAlignmentEnd && homopolymer_start.has_value() &&
      homopolymer_start.value() <= info.ref_pos) {
    // deletion is in a homopolymer region
    return false;
  }
  if (has_vcf_feats && !vcf_feat_positions->contains(info.ref_pos - 1)) {
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
  if (info.ref_pos < region->start || info.ref_pos >= region->end) {
    // mismatch out of bounds
    return false;
  }
  if (params->filter_homopolymer == HomopolymerFilter::kAlignmentEnd && homopolymer_start.has_value() &&
      homopolymer_start.value() <= info.ref_pos) {
    // mismatch is in a homopolymer region
    return false;
  }
  if (has_vcf_feats && !vcf_feat_positions->contains(info.ref_pos)) {
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
  if (ref_pos > region->end) {
    // match out of bounds
    return false;
  }
  if (has_vcf_feats && !vcf_feat_positions->contains(ref_pos) && !vcf_feat_positions->contains(ref_pos - 1)) {
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
  return !near_non_concordant.empty() && near_non_concordant.at(info.read_pos);
}

bool AlignContext::IsMismatchNearNonConcordantBase(const AlignOpInfo& info) const {
  return !near_non_concordant.empty() && near_non_concordant.at(info.read_pos);
}

bool AlignContext::IsMatchNearNonConcordantBase(const u64 read_pos) const {
  return !near_non_concordant.empty() && near_non_concordant.at(read_pos);
}

void ProcessInsertion(BamRegionFeatureCollection& features, const AlignContext& align_ctx, const AlignOpInfo& info) {
  const std::span<const u8> quals{bam_get_qual(align_ctx.bam_record),
                                  static_cast<size_t>(align_ctx.bam_record->core.l_qseq)};
  const auto op_quals = quals.subspan(info.read_pos, info.op_len);
  const u32 baseq_sum = std::accumulate(op_quals.begin(), op_quals.end(), u32{0});
  const f32 baseq_mean = static_cast<f32>(baseq_sum) / static_cast<f32>(info.op_len);
  if (baseq_mean < static_cast<f32>(align_ctx.params->min_bq)) {
    return;
  }
  if (align_ctx.bam_record->core.pos < 0) {
    // Should not get here, check should be done in upstream call for this. Keep just incase call made without check
    WarnAsErrorIfSet("Found alignment with start position < 0: '{}'", bam_get_qname(align_ctx.bam_record));
    return;
  }
  const auto align_pos = static_cast<u64>(align_ctx.bam_record->core.pos);
  const auto min_distance = std::min(info.ref_pos - align_pos, align_ctx.ref_end - info.ref_pos - 1);
  if (static_cast<f32>(min_distance) < align_ctx.params->min_allowed_distance_from_end) {
    return;
  }

  // no `-1` in position here because insertion starts at the previous matching base
  const char ref_base{(*align_ctx.ref_seq)[info.ref_pos]};
  if (IsNotACTG(ref_base)) {
    // Since we are not confident in variant calls at non-ACGT reference positions,
    // we do not extract variant features at this position.
    return;
  }

  const auto* const begin_seq = bam_get_seq(align_ctx.bam_record);
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
  const VariantId vid(align_ctx.region->chrom, info.ref_pos, ref_allele, alt_allele);

  // update general features
  UnifiedVariantFeature& feat = features.var_features[vid];
  if (feat.support == 0) {
    feat.context = Get1bpContext((*align_ctx.ref_seq), info.ref_pos);
    feat.context_index = UnifiedVariantFeature::ContextIndex(feat.context);
  }
  IncrementFeature(
      feat, align_ctx, baseq_mean, min_distance, align_ctx.IsInsertionNearNonConcordantBase(info), info.read_pos);

  if (align_ctx.is_tumor_sample.has_value()) {
    UnifiedVariantFeature& tn_feat =
        align_ctx.is_tumor_sample.value() ? features.tumor_var_features[vid] : features.normal_var_features[vid];
    if (tn_feat.support == 0) {
      tn_feat.context = Get1bpContext((*align_ctx.ref_seq), info.ref_pos);
      tn_feat.context_index = UnifiedVariantFeature::ContextIndex(tn_feat.context);
    }
    IncrementFeature(
        tn_feat, align_ctx, baseq_mean, min_distance, align_ctx.IsInsertionNearNonConcordantBase(info), info.read_pos);
  }
}

void ProcessDeletion(BamRegionFeatureCollection& features, const AlignContext& align_ctx, const AlignOpInfo& info) {
  if (info.read_pos == 0) {
    return;
  }

  const u8 baseq = bam_get_qual(align_ctx.bam_record)[info.read_pos - 1];
  if (baseq < align_ctx.params->min_bq) {
    return;
  }

  const auto min_distance = std::min(info.ref_pos - static_cast<u32>(align_ctx.bam_record->core.pos),
                                     align_ctx.ref_end - (info.ref_pos + info.op_len));
  if (static_cast<f32>(min_distance) < align_ctx.params->min_allowed_distance_from_end) {
    return;
  }

  const std::string ref_allele{(*align_ctx.ref_seq).substr(info.ref_pos - 1, info.op_len + 1)};
  if (IsAnyNotACTG(ref_allele)) {
    // Since we are not confident in variant calls at non-ACGT reference positions,
    // we do not extract variant features at this position.
    return;
  }

  const std::string alt_allele{(*align_ctx.ref_seq)[info.ref_pos - 1]};
  const VariantId vid(align_ctx.region->chrom, info.ref_pos - 1, ref_allele, alt_allele);

  // update generate features
  UnifiedVariantFeature& feat = features.var_features[vid];
  if (feat.support == 0) {
    feat.context = Get1bpContext((*align_ctx.ref_seq), info.ref_pos);
    feat.context_index = UnifiedVariantFeature::ContextIndex(feat.context);
  }
  IncrementFeature(
      feat, align_ctx, baseq, min_distance, align_ctx.IsDeletionNearNonConcordantBase(info), info.read_pos - 1);

  if (align_ctx.is_tumor_sample.has_value()) {
    UnifiedVariantFeature& tn_feat =
        align_ctx.is_tumor_sample.value() ? features.tumor_var_features[vid] : features.normal_var_features[vid];
    if (tn_feat.support == 0) {
      tn_feat.context = Get1bpContext((*align_ctx.ref_seq), info.ref_pos);
      tn_feat.context_index = UnifiedVariantFeature::ContextIndex(tn_feat.context);
    }
    IncrementFeature(
        tn_feat, align_ctx, baseq, min_distance, align_ctx.IsDeletionNearNonConcordantBase(info), info.read_pos - 1);
  }
}

void ProcessMismatch(BamRegionFeatureCollection& features, const AlignContext& align_ctx, const AlignOpInfo& info) {
  const u8 baseq = bam_get_qual(align_ctx.bam_record)[info.read_pos];
  if (baseq < align_ctx.params->min_bq) {
    return;
  }

  const auto min_distance =
      std::min(info.ref_pos - static_cast<u32>(align_ctx.bam_record->core.pos), align_ctx.ref_end - info.ref_pos - 1);
  if (static_cast<f32>(min_distance) < align_ctx.params->min_allowed_distance_from_end) {
    return;
  }

  const char ref_base{(*align_ctx.ref_seq)[info.ref_pos]};
  if (IsNotACTG(ref_base)) {
    // Since we are not confident in variant calls at non-ACGT reference positions,
    // we do not extract variant features at this position.
    return;
  }

  const auto* const begin_seq = bam_get_seq(align_ctx.bam_record);
  const char base = seq_nt16_str[bam_seqi(begin_seq, info.read_pos)];
  if (IsNotACTG(base)) {
    return;
  }
  const std::string alt_allele{base};
  const std::string ref_allele{ref_base};
  const VariantId vid(align_ctx.region->chrom, info.ref_pos, ref_allele, alt_allele);

  // update general features
  UnifiedVariantFeature& feat = features.var_features[vid];
  if (feat.support == 0) {
    feat.context = Get1bpContext((*align_ctx.ref_seq), info.ref_pos);
    feat.context_index = UnifiedVariantFeature::ContextIndex(feat.context);
  }
  IncrementFeature(
      feat, align_ctx, baseq, min_distance, align_ctx.IsMismatchNearNonConcordantBase(info), info.read_pos);

  if (align_ctx.is_tumor_sample.has_value()) {
    UnifiedVariantFeature& tn_feat =
        align_ctx.is_tumor_sample.value() ? features.tumor_var_features[vid] : features.normal_var_features[vid];
    if (tn_feat.support == 0) {
      tn_feat.context = Get1bpContext((*align_ctx.ref_seq), info.ref_pos);
      tn_feat.context_index = UnifiedVariantFeature::ContextIndex(tn_feat.context);
    }
    IncrementFeature(
        tn_feat, align_ctx, baseq, min_distance, align_ctx.IsMismatchNearNonConcordantBase(info), info.read_pos);
  }
}

namespace {
/**
 * @brief Update the non-homopolymer match features in the UnifiedReferenceFeature.
 * @param feat The UnifiedReferenceFeature to update.
 * @param mapq The mapping quality of the read.
 * @param baseq The base quality of the read at the position.
 * @param weighted_depth The weighted depth of the read at this position.
 */
void UpdateMatchNonHomopolymerFeatures(UnifiedReferenceFeature& feat,
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

/**
 * @brief Update the match features in the UnifiedReferenceFeature.
 * @param feat The UnifiedReferenceFeature to update.
 * @param align_ctx The AlignContext containing alignment information.
 * @param baseq The base quality of the match.
 * @param min_distance The minimum distance from the match to the alignment end's reference position.
 * @param near_non_concordant Whether the match is near a non-concordant base.
 * @param not_in_homopolymer Whether the match is not in a homopolymer region.
 */
void UpdateMatchFeatures(UnifiedReferenceFeature& feat,
                         const AlignContext& align_ctx,
                         const u8 baseq,
                         const u64 min_distance,
                         const bool near_non_concordant,
                         const bool not_in_homopolymer,
                         const u64 read_pos) {
  const u8 mapq{align_ctx.bam_record->core.qual};
  const bool is_simplex_protocol{align_ctx.params->sequencing_protocol == SequencingProtocol::kSimplex};

  if (align_ctx.read_type == ReadType::kSimplex) {
    ++feat.simplex;
    feat.mapq_sum_simplex += mapq;
    feat.distance_sum_simplex += min_distance;
    if (!is_simplex_protocol) {
      return;
    }
  }

  if (!is_simplex_protocol && near_non_concordant) {
    if (align_ctx.plus_counts > 0 && align_ctx.minus_counts > 0) {
      feat.duplex_lowbq += kDuplexLowbqIncrement;
    }
    feat.mapq_sum_lowbq += mapq;
    feat.distance_sum_lowbq += min_distance;
    return;
  }

  if (!is_simplex_protocol) {
    IncrementDuplexBaseType(feat, align_ctx, read_pos);
  }

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

  const f64 weighted_depth{static_cast<f64>(baseq) / kMaxPossibleBaseQual * mapq / kMaxPossibleMapQual};
  feat.weighted_depth += weighted_depth;

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

  const u32 familysize{is_simplex_protocol ? kSimplexReadFamilySize : align_ctx.plus_counts + align_ctx.minus_counts};
  feat.familysize_sum += familysize;
  if (familysize < kFamilySizeCountThreshold5) {
    ++feat.familysize_lt5_count;
    if (familysize < kFamilySizeCountThreshold3) {
      ++feat.familysize_lt3_count;
    }
  }

  if (not_in_homopolymer) {
    // this matching base is not in a homopolymer
    UpdateMatchNonHomopolymerFeatures(feat, mapq, baseq, weighted_depth);
  }
}
}  // namespace

void ProcessMatch(BamRegionFeatureCollection& features,
                  const AlignContext& align_ctx,
                  const u64 read_pos,
                  const u64 ref_pos) {
  const u8 baseq = bam_get_qual(align_ctx.bam_record)[read_pos];
  if (baseq < align_ctx.params->min_bq) {
    return;
  }

  const auto min_distance =
      std::min(ref_pos - static_cast<u32>(align_ctx.bam_record->core.pos), align_ctx.ref_end - ref_pos - 1);
  if (static_cast<f32>(min_distance) < align_ctx.params->min_allowed_distance_from_end) {
    return;
  }

  UnifiedReferenceFeature& feat = features.ref_features[ref_pos];
  if (ContainsReadIdInsertIfNot(feat.read_ids, align_ctx.read_id)) {
    return;
  }

  const bool near_non_concordant = align_ctx.IsMatchNearNonConcordantBase(read_pos);
  const bool not_in_homopolymer = !align_ctx.homopolymer_start.has_value() || align_ctx.homopolymer_start > ref_pos;

  // update general features
  UpdateMatchFeatures(feat, align_ctx, baseq, min_distance, near_non_concordant, not_in_homopolymer, read_pos);

  if (align_ctx.is_tumor_sample.has_value()) {
    UnifiedReferenceFeature& tn_feat = align_ctx.is_tumor_sample.value() ? features.tumor_ref_features[ref_pos]
                                                                         : features.normal_ref_features[ref_pos];
    UpdateMatchFeatures(tn_feat, align_ctx, baseq, min_distance, near_non_concordant, not_in_homopolymer, read_pos);
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
                        const bam1_t* const record,
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
        const auto* const read = bam_get_seq(record);
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
          const auto* const read = bam_get_seq(record);
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

namespace {
/**
 * @brief Find nearby non-ACGT bases in the sequence and mark them in the boolean vector.
 * @param result Vector indicating whether each base is at or nearby a non-ACGT base
 * @param seq Sequence string
 * @note result is a boolean vector of the same length as seq, where each position is set to true if it is at or
 *       nearby a non-ACGT base, and false otherwise.
 */
void NearNonACGT(const std::string& seq, vec<bool>& result) {  // NOSONAR: vec<bool> is intentional
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
}  // namespace

vec<bool> FindNearbyNonConcordantBase(  // NOSONAR(cpp:M23_357) — vec<bool> intentional for cache locality
    const bam1_t* const alignment) {
  const auto seq_len = static_cast<u64>(alignment->core.l_qseq);
  vec<bool> near_n(seq_len, false);  // NOSONAR(cpp:M23_357) — vec<bool> intentional for cache locality
  std::string seq(seq_len, 'N');
  // find N using the read sequence and base quality
  const auto* const begin_seq = bam_get_seq(alignment);
  const std::span<const u8> qual_data{bam_get_qual(alignment), seq_len};
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

vec<bool> FindNearbyNonConcordantBase(  // NOSONAR(cpp:M23_357) — vec<bool> intentional for cache locality
    const bam1_t* const alignment,
    const vec<BaseType>& base_types) {
  if (base_types.empty()) {
    // no base types available, fall back to using base quality
    return FindNearbyNonConcordantBase(alignment);
  }
  const auto seq_len = static_cast<u64>(alignment->core.l_qseq);
  vec<bool> result(seq_len, false);  // NOSONAR: vec<bool> retained for cache locality on read-length vector
  std::string seq(seq_len, 'N');
  // find N using the read sequence and base types
  const auto* const seq_start = bam_get_seq(alignment);
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

void ProcessAlignment(const AlignContext& align_ctx, BamRegionFeatureCollection& bam_feats) {
  const auto num_infos{(*align_ctx.align_op_infos).size()};
  for (size_t a = 0; a < num_infos; ++a) {
    const auto& info = (*align_ctx.align_op_infos).at(a);
    switch (info.op) {
      case kInsertion: {
        if (align_ctx.CanProcessInsertion(info)) {
          ProcessInsertion(bam_feats, align_ctx, info);
        }
        break;
      }
      case kDeletion: {
        // The deletion must overlap this region, but it cannot overlap the previous region.
        // If the deletion also overlaps the previous region, then it should be processed there instead of this region.
        // This ensures that a deletion spanning across regions is only processed in its left-most overlapping region,
        // therefore preventing additional incomplete features for the same deletion.
        if (align_ctx.CanProcessDeletion(info)) {
          ProcessDeletion(bam_feats, align_ctx, info);
        }
        break;
      }
      case kMismatch: {
        if (align_ctx.CanProcessMismatch(info)) {
          ProcessMismatch(bam_feats, align_ctx, info);
        }
        break;
      }
      case kMatch: {
        // We skip reference feature extraction at a position if the following base is an insertion because
        // this position is part of the ALT allele, i.e. the insertion occurs at this position.
        const u64 end{a + 1 < num_infos && (*align_ctx.align_op_infos)[a + 1].op == kInsertion ? info.op_len - 1
                                                                                               : info.op_len};
        for (u64 i = 0; i < end; ++i) {
          // Do not skip positions before `region_start` because these positions can be reference features for
          // deletions overlapping `region_start`.
          // Use `<=` instead of `<` to extract reference features for deletions that start at `region_end - 1`.
          const auto ref_pos = info.ref_pos + i;
          const auto read_pos = info.read_pos + i;
          if (align_ctx.CanProcessMatch(read_pos, ref_pos)) {
            ProcessMatch(bam_feats, align_ctx, read_pos, ref_pos);
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
