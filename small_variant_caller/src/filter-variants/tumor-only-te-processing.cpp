#include "tumor-only-te-processing.h"

#include <algorithm>
#include <limits>

#include <xoos/util/container-functions.h>
#include <xoos/util/string-functions.h>

#include "core/sequencing-protocol.h"
#include "core/variant-feature-extraction.h"
#include "core/vcf-fields.h"
#include "shap-value-tsv.h"

namespace xoos::svc {

namespace {
/**
 * @brief Convert regions to a map of chromosome to set of positions.
 * @param regions Map of chromosome to vector of intervals
 * @return Map of chromosome to set of positions
 */
StrUnorderedMap<std::unordered_set<u64>> GetChromToPositionSet(const StrMap<vec<Interval>>& regions) {
  StrUnorderedMap<std::unordered_set<u64>> result;
  for (const auto& [chrom, intervals] : regions) {
    for (const auto& [start, end] : intervals) {
      for (auto pos = start; pos < end; ++pos) {
        result[chrom].emplace(pos);
      }
    }
  }
  return result;
}

/**
 * @brief Check if a chromosome and position exists in the given map.
 * @param chrom_to_positions Map of chromosome to set of positions
 * @param chrom Chromosome name
 * @param pos Position to check
 * @return True if the chromosome and position exists, false otherwise
 */
bool ContainsChromPosition(const StrUnorderedMap<std::unordered_set<u64>>& chrom_to_positions,
                           const std::string& chrom,
                           const u64 pos) {
  const auto it = chrom_to_positions.find(chrom);
  if (it == chrom_to_positions.end()) {
    return false;
  }
  return it->second.contains(pos);
}

/// Options for UpdateRecord / UpdatePhasedRecord to avoid adjacent bool parameters.
struct RecordUpdateOptions {
  bool is_hotspot{false};
  bool is_duplex_protocol{false};
};

/**
 * @brief Update a VCF record with selected feature information. Intended for `somatic` workflow only.
 * @param record VCF record to be updated
 * @param bam_feat BamFeatureTuple containing variant and reference features
 * @param opts Update options (hotspot flag, duplex protocol flag)
 */
void UpdateRecord(const io::VcfRecordPtr& record, const BamFeatureTuple& bam_feat, const RecordUpdateOptions& opts) {
  const auto& var_feat = bam_feat.var_feat;
  const auto& ref_feat = bam_feat.ref_feat;
  record->SetFormatField(kWeightedCountsId, vec<f32>{static_cast<f32>(var_feat.weighted_score)});
  record->SetFormatField(kMachineLearningId, vec<f32>{static_cast<f32>(var_feat.ml_score)});
  record->SetFormatField(kNonDuplexCountsId, vec<s32>{static_cast<s32>(var_feat.nonduplex)});
  record->SetFormatField(kDuplexCountsId, vec<s32>{static_cast<s32>(var_feat.duplex)});
  record->SetFormatField(kStrandBiasMetricId, vec<f32>{static_cast<f32>(var_feat.strandbias)});
  record->SetFormatField(kPlusOnlyCountsId, vec<s32>{static_cast<s32>(var_feat.plusonly)});
  record->SetFormatField(kMinusOnlyCountsId, vec<s32>{static_cast<s32>(var_feat.minusonly)});
  record->SetFormatField(kSequenceContextId, vec<std::string>{var_feat.context});
  // update the existing AD field, since we use different read filtering than Mutect2
  record->SetFormatField(
      kAlleleDepthId,
      vec<s32>{static_cast<s32>(ref_feat.support), static_cast<s32>(var_feat.duplex + var_feat.nonduplex)});
  if (opts.is_hotspot) {
    record->SetFormatField(kHotspotId, vec<s32>{1});
  }
  if (opts.is_duplex_protocol) {
    record->SetFormatField(
        kDuplexConcordantCountsId,
        vec<s32>{static_cast<s32>(ref_feat.duplex_concordant), static_cast<s32>(var_feat.duplex_concordant)});
    record->SetFormatField(
        kDuplexSimplexCountsId,
        vec<s32>{static_cast<s32>(ref_feat.duplex_simplex), static_cast<s32>(var_feat.duplex_simplex)});
    record->SetFormatField(
        kDuplexDiscordantCountsId,
        vec<s32>{static_cast<s32>(ref_feat.duplex_discordant), static_cast<s32>(var_feat.duplex_discordant)});
    record->SetFormatField(kDuplexLowbqCountsId,
                           vec<s32>{static_cast<s32>(std::round(ref_feat.duplex_lowbq * kDuplexLowbqToCountFactor)),
                                    static_cast<s32>(std::round(var_feat.duplex_lowbq * kDuplexLowbqToCountFactor))});
  }
}

/**
 * @brief Update a phased VCF record with empty values. Intended for `somatic` workflow only.
 * @param record VCF record to be updated
 * @param opts Update options (duplex protocol flag)
 */
void UpdatePhasedRecord(const io::VcfRecordPtr& record, const RecordUpdateOptions& opts) {
  // Phased calls don't use ML based filtering and do not have values for these fields
  record->SetFormatField(kWeightedCountsId, vec<f32>{static_cast<f32>(0.0)});
  record->SetFormatField(kMachineLearningId, vec<f32>{static_cast<f32>(0.0)});
  record->SetFormatField(kNonDuplexCountsId, vec<s32>{0});
  record->SetFormatField(kDuplexCountsId, vec<s32>{0});
  record->SetFormatField(kStrandBiasMetricId, vec<f32>{static_cast<f32>(0.0)});
  record->SetFormatField(kPlusOnlyCountsId, vec<s32>{0});
  record->SetFormatField(kMinusOnlyCountsId, vec<s32>{0});
  record->SetFormatField(kSequenceContextId, vec<std::string>{"NA"});
  if (opts.is_duplex_protocol) {
    record->SetFormatField(kDuplexConcordantCountsId, vec<s32>{0, 0});
    record->SetFormatField(kDuplexSimplexCountsId, vec<s32>{0, 0});
    record->SetFormatField(kDuplexDiscordantCountsId, vec<s32>{0, 0});
    record->SetFormatField(kDuplexLowbqCountsId, vec<s32>{0, 0});
  }
}

/**
 * @brief Set up a VCF record for forced variant calls.
 * @param record VCF record to be updated
 * @param variant_id VariantId struct
 */
void SetupForceRecord(const io::VcfRecordPtr& record, const VariantId& variant_id) {
  const vec<std::string> alleles = {variant_id.ref, variant_id.alt};
  record->SetAlleles(alleles);
  record->SetFilter(kFilteringForcedId);
}
}  // namespace

TumorOnlyTeProcessing::TumorOnlyTeProcessing(const GlobalContext& global_ctx,
                                             const FilterSettings& filter_settings,
                                             const PhasedFilterSettings& phased_filter_settings,
                                             const u64 prev_pos,
                                             VarIdToVcfFeatures vcf_features,
                                             BamRegionFeatureCollection bam_features)
    : _filter_settings(&filter_settings),
      _phased_filter_settings(&phased_filter_settings),
      _prev_pos(prev_pos),
      _vcf_features(std::move(vcf_features)),
      _bam_features(std::move(bam_features)),
      _hdr(global_ctx.hdr),
      _phased(global_ctx.phased),
      _is_duplex_protocol(IsDuplexProtocol(global_ctx.bam_feat_params.sequencing_protocol)) {
  // Set forced call positions
  _forced_call_chrom_pos = global_ctx.force_calls.has_value() ? GetChromToPositionSet(global_ctx.force_calls.value())
                                                              : StrUnorderedMap<std::unordered_set<u64>>{};
  if (!_forced_call_chrom_pos.empty()) {
    // Build map of chrom->pos->VariantId for fast lookup when creating forced calls
    for (const auto& [vid, feat] : _bam_features.var_features) {
      if (ContainsChromPosition(_forced_call_chrom_pos, vid.chrom, vid.pos)) {
        // only need to track variants at forced call positions
        _forced_call_chrom_pos_to_vids[vid.chrom][vid.pos].emplace_back(vid);
      }
    }
  }
}

std::optional<BamFeatureTuple> TumorOnlyTeProcessing::ScoreVariant(const VariantId& vid,
                                                                   const ChromMedianDepth& normalize_targets,
                                                                   const vec<FeatureColumn>& scoring_cols,
                                                                   const ScoreCalculator& somatic_calculator,
                                                                   vec<vec<std::string>>& shap_value_rows) {
  // Not using `TumorNormalBamFeatureRow` here because there is no distinction for tumor/normal in tumor-only TE
  std::optional<BamFeatureTuple> bam_feat = _bam_features.GetBamFeatureTuple(vid);

  if (bam_feat.has_value()) {
    const auto vcf_features_it = _vcf_features.find(vid);
    const VcfFeature& vcf_feat = vcf_features_it == _vcf_features.end() ? kZeroVcfFeature : vcf_features_it->second;

    const DepthTuple& normalize_target = normalize_targets.GetValue(vid.chrom);

    const auto feature_vec = GetFeatureVec(scoring_cols, vid, bam_feat.value(), vcf_feat, normalize_target);
    const auto score = somatic_calculator.CalculateScore(feature_vec);
    bam_feat->var_feat.ml_score = score.probability;
    if (!score.shap_values.empty()) {
      shap_value_rows.emplace_back(AssembleShapValueTsvRow(vid, score));
    }
    bam_feat->var_feat.filter_status = FilterVariant(vid, bam_feat->var_feat, bam_feat->ref_feat, *_filter_settings);
    return bam_feat;
  }
  return std::nullopt;
}

void TumorOnlyTeProcessing::CreateNewForceRecordsForInbetweenPositions(const u64 pos,
                                                                       const std::string& chrom,
                                                                       vec<io::VcfRecordPtr>& out_records) {
  // For forced calling we must check every position from the last seen variant position up to the current record's
  // position. Any forced calls that need to be made should be done so before dealing with the current record to keep
  // the sorted order intact.
  const u64 max_pos_to_check = pos;

  const auto chrom_it = _forced_call_chrom_pos_to_vids.find(chrom);
  if (chrom_it == _forced_call_chrom_pos_to_vids.end()) {
    _prev_pos = pos;
    return;
  }

  const auto& pos_to_vids = chrom_it->second;
  for (u64 i = _prev_pos; i < max_pos_to_check; ++i) {
    if (ContainsChromPosition(_seen_forced_call_chrom_pos, chrom, i)) {
      // already seen this forced call position
      continue;
    }
    const auto& pos_it = pos_to_vids.find(i);
    if (pos_it == pos_to_vids.end()) {
      // no variants at this forced call position
      continue;
    }
    // create forced call records for all variants at this position
    for (const auto& vid : pos_it->second) {
      const std::optional<BamFeatureTuple> bam_feat = _bam_features.GetBamFeatureTuple(vid);
      if (!bam_feat.has_value()) {
        // this should not happen because we built the map from existing features, but just in case
        continue;
      }
      const io::VcfRecordPtr new_record = io::VcfRecord::CreateFromHeader(_hdr);
      new_record->SetChromosome(chrom);
      new_record->SetPosition(static_cast<s32>(i));
      SetupForceRecord(new_record, vid);
      const auto key_to_check = GetVariantCorrelationKey(chrom, i, vid.ref, vid.alt, false);
      const auto is_hotspot = _filter_settings->hotspots.contains(key_to_check);
      UpdateRecord(new_record, bam_feat.value(), {.is_hotspot = is_hotspot, .is_duplex_protocol = _is_duplex_protocol});
      out_records.emplace_back(new_record);
    }
    _seen_forced_call_chrom_pos[chrom].insert(i);
  }
  _prev_pos = pos;
}

void TumorOnlyTeProcessing::ProcessRecord(const VariantId& vid,
                                          const io::VcfRecordPtr& record,
                                          const std::optional<BamFeatureTuple>& bam_feature,
                                          vec<io::VcfRecordPtr>& out_records,
                                          const bool make_forcedcalls,
                                          const std::string& key,
                                          const std::string& key_unpadded) {
  bool phased_record = false;
  if (_phased && record->HasFormatFieldNoCheck("PID")) {
    const auto pid_values = record->GetFormatField<std::string>("PID");
    if (!pid_values.empty()) {
      phased_record = true;
    }
  }
  auto record_copy = record->Clone(_hdr);
  if (bam_feature.has_value()) {
    // Calculate variant quality based on the ML score alone because this is a binary classification.
    record_copy->SetQuality(static_cast<s32>(CalculatePhredQualityScore(bam_feature->var_feat.ml_score)));
  }
  // Set GQ to 0 in VCF record because the genotype is unclear.
  record_copy->SetFormatField(kFieldGq, vec<s32>{0});

  bool found = false;
  if (!phased_record && bam_feature.has_value()) {
    for (const auto& item : bam_feature->var_feat.filter_status) {
      record_copy->AddFilter(item);
    }
    const auto is_hotspot = _filter_settings->hotspots.contains(key_unpadded);
    UpdateRecord(
        record_copy, bam_feature.value(), {.is_hotspot = is_hotspot, .is_duplex_protocol = _is_duplex_protocol});
    out_records.emplace_back(record_copy);
    found = true;
  } else if (phased_record) {
    // Phased variant
    found = true;
    vec<s32> n_ad = record->GetFormatField<s32>("AD");
    vec<f32> n_af = record->GetFormatField<f32>("AF");
    vec<s32> n_mq = record->GetInfoField<s32>("MMQ");
    vec<s32> n_bq = record->GetInfoField<s32>("MBQ");
    const vec<std::string> filters =
        FilterPhasedVariant(key,
                            (*_phased_filter_settings),
                            static_cast<u32>(n_ad[1]),
                            n_af[0],
                            static_cast<u8>(std::clamp(
                                *std::ranges::max_element(n_mq), 0, static_cast<s32>(std::numeric_limits<u8>::max()))),
                            static_cast<u8>(std::clamp(
                                *std::ranges::max_element(n_bq), 0, static_cast<s32>(std::numeric_limits<u8>::max()))));
    for (const auto& item : filters) {
      record_copy->AddFilter(item);
    }
    UpdatePhasedRecord(record_copy, {.is_duplex_protocol = _is_duplex_protocol});
    out_records.emplace_back(record_copy);
  } else {
    // variant has no BAM features and is not phased; no output
  }
  if (make_forcedcalls && found && ContainsChromPosition(_forced_call_chrom_pos, vid.chrom, vid.pos)) {
    _seen_forced_call_chrom_pos[vid.chrom].insert(vid.pos);
  }
}

}  // namespace xoos::svc
