#include "tumor-normal-processing.h"

#include "core/filtering.h"
#include "core/sequencing-protocol.h"
#include "core/variant-feature-extraction.h"
#include "core/vcf-fields.h"
#include "shap-value-tsv.h"
#include "util/seq-util.h"

namespace xoos::svc {

TumorNormalProcessing::TumorNormalProcessing(const GlobalContext& global_ctx,
                                             VarIdToVcfFeatures vcf_features,
                                             BamRegionFeatureCollection bam_features)
    : _vcf_features(std::move(vcf_features)),
      _bam_features(std::move(bam_features)),
      _hdr(global_ctx.hdr),
      _snv_min_ml_score(global_ctx.snv_min_ml_score),
      _indel_min_ml_score(global_ctx.indel_min_ml_score),
      _min_tumor_support(global_ctx.min_tumor_support),
      _max_normal_support(global_ctx.max_normal_support),
      _min_tumor_af(global_ctx.min_tumor_af),
      _catalog_min_tumor_af(global_ctx.catalog_min_tumor_af),
      _catalog_snv_min_ml_score(global_ctx.catalog_snv_min_ml_score),
      _catalog_indel_min_ml_score(global_ctx.catalog_indel_min_ml_score),
      _min_dp_ratio(global_ctx.min_dp_ratio),
      _max_indel_size(global_ctx.max_indel_size),
      _normalize_scoring_features(global_ctx.normalize_scoring_features),
      _scoring_cols(global_ctx.model_config.scoring_cols),
      _is_duplex_protocol(IsDuplexProtocol(global_ctx.bam_feat_params.sequencing_protocol)) {
}

namespace {
/**
 * Function sets the required INFO and FORMAT values for tumor-normal-wgs VCF records before they are output. A set
 * genotype is assigned and values in the VCF record updated to reflect the values observed in both the VCF record and
 * support alignments from the BAM file. The Tumor-Normal somatic model score is added as well.
 * @param record A VCFRecordPtr reference
 * @param variant_id A VariantId struct reference
 * @param bam_feat A TumorNormalBamFeatureRow struct reference
 * @param vcf_feat A VcfFeature struct reference
 * @param tumor_first Boolean, True if the tumor is sample is the first sample in the VCF. Assumes only one tumor and
 * one normal in VCF.
 * @param pred_score The predicted score from the tumor-normal somatic model for the variant.
 */
void SetSomaticTNValues(const io::VcfRecordPtr& record,
                        const VariantId& variant_id,
                        const TumorNormalBamFeatureTuple& bam_feat,
                        const VcfFeature& vcf_feat,
                        const bool tumor_first,
                        const PredictionScore& pred_score,
                        const bool is_duplex_protocol) {
  // Setup the format and info fields for tumor normal vcf records. Assumes only two samples in the VCF, one tumor, one
  // normal
  const auto tumor_support = static_cast<s32>(bam_feat.tumor_var_feat.support);
  const auto normal_support = static_cast<s32>(bam_feat.normal_var_feat.support);
  const auto tumor_ref_support = static_cast<s32>(bam_feat.tumor_ref_feat.support);
  const auto normal_ref_support = static_cast<s32>(bam_feat.normal_ref_feat.support);
  record->SetFormatField(kFieldAd,
                         tumor_first ? vec<s32>{tumor_ref_support, tumor_support, normal_ref_support, normal_support}
                                     : vec<s32>{normal_ref_support, normal_support, tumor_ref_support, tumor_support});
  const auto tumor_af = static_cast<f32>(bam_feat.tumor_var_feat.duplex_af);
  const auto normal_af = static_cast<f32>(bam_feat.normal_var_feat.duplex_af);
  record->SetFormatField(kFieldAf, tumor_first ? vec<f32>{tumor_af, normal_af} : vec<f32>{normal_af, tumor_af});
  record->SetFormatField(kFieldDp,
                         tumor_first
                             ? vec<s32>{tumor_support + tumor_ref_support, normal_support + normal_ref_support}
                             : vec<s32>{normal_support + normal_ref_support, tumor_support + tumor_ref_support});
  const auto tumor_bq = static_cast<f32>(bam_feat.tumor_var_feat.baseq_mean);
  const auto normal_bq = static_cast<f32>(bam_feat.normal_var_feat.baseq_mean);
  record->SetFormatField(kBaseqQualId, tumor_first ? vec<f32>{tumor_bq, normal_bq} : vec<f32>{normal_bq, tumor_bq});
  const auto tumor_mq = static_cast<f32>(bam_feat.tumor_var_feat.mapq_mean);
  const auto normal_mq = static_cast<f32>(bam_feat.normal_var_feat.mapq_mean);
  record->SetFormatField(kMapQualId, tumor_first ? vec<f32>{tumor_mq, normal_mq} : vec<f32>{normal_mq, tumor_mq});
  const auto tumor_distance = static_cast<f32>(bam_feat.tumor_var_feat.distance_mean);
  const auto normal_distance = static_cast<f32>(bam_feat.normal_var_feat.distance_mean);
  record->SetFormatField(
      kDistanceId, tumor_first ? vec<f32>{tumor_distance, normal_distance} : vec<f32>{normal_distance, tumor_distance});
  record->SetFormatField(kFieldGq, vec<s32>{kTumorNormalGq, kTumorNormalGq});
  record->SetFormatField(kMachineLearningId,
                         tumor_first ? vec<f32>{static_cast<f32>(pred_score.probability), 0}
                                     : vec<f32>{0, static_cast<f32>(pred_score.probability)});
  // Set variant quality based on the ML score alone because this is a binary classification.
  record->SetQuality(pred_score.variant_quality);
  record->SetInfoField(kRefBQId, vec<f32>{static_cast<f32>(bam_feat.ref_feat.baseq_mean)});
  record->SetInfoField(kRefMQId, vec<f32>{static_cast<f32>(bam_feat.ref_feat.mapq_mean)});
  record->SetInfoField(kAltBQId, vec<f32>{static_cast<f32>(bam_feat.var_feat.baseq_mean)});
  record->SetInfoField(kAltMQId, vec<f32>{static_cast<f32>(bam_feat.var_feat.mapq_mean)});
  record->SetInfoField(kFieldNalod, vec<f32>{vcf_feat.nalod});
  record->SetInfoField(kFieldNlod, vec<f32>{vcf_feat.nlod});
  record->SetInfoField(kFieldTlod, vec<f32>{vcf_feat.tlod});
  record->SetInfoField(kFieldMpos, vec<s32>{static_cast<s32>(vcf_feat.mpos)});
  record->SetInfoField(kSubtypeId, vec<s32>{SubstIndex(variant_id.ref, variant_id.alt)});
  record->SetInfoField(kContextId, vec<std::string>{bam_feat.var_feat.context});
  record->SetInfoField(kFieldPopaf, vec<f32>{vcf_feat.popaf});
  if (is_duplex_protocol) {
    const auto& tumor_var = bam_feat.tumor_var_feat;
    const auto& normal_var = bam_feat.normal_var_feat;
    const auto& tumor_ref = bam_feat.tumor_ref_feat;
    const auto& normal_ref = bam_feat.normal_ref_feat;
    const auto t_adc_ref = static_cast<s32>(tumor_ref.duplex_concordant);
    const auto t_adc = static_cast<s32>(tumor_var.duplex_concordant);
    const auto n_adc_ref = static_cast<s32>(normal_ref.duplex_concordant);
    const auto n_adc = static_cast<s32>(normal_var.duplex_concordant);
    const auto t_ads_ref = static_cast<s32>(tumor_ref.duplex_simplex);
    const auto t_ads = static_cast<s32>(tumor_var.duplex_simplex);
    const auto n_ads_ref = static_cast<s32>(normal_ref.duplex_simplex);
    const auto n_ads = static_cast<s32>(normal_var.duplex_simplex);
    const auto t_add_ref = static_cast<s32>(tumor_ref.duplex_discordant);
    const auto t_add = static_cast<s32>(tumor_var.duplex_discordant);
    const auto n_add_ref = static_cast<s32>(normal_ref.duplex_discordant);
    const auto n_add = static_cast<s32>(normal_var.duplex_discordant);
    const auto t_adl_ref = static_cast<s32>(std::round(tumor_ref.duplex_lowbq * kDuplexLowbqToCountFactor));
    const auto t_adl = static_cast<s32>(std::round(tumor_var.duplex_lowbq * kDuplexLowbqToCountFactor));
    const auto n_adl_ref = static_cast<s32>(std::round(normal_ref.duplex_lowbq * kDuplexLowbqToCountFactor));
    const auto n_adl = static_cast<s32>(std::round(normal_var.duplex_lowbq * kDuplexLowbqToCountFactor));
    // Number=R with 2 samples: {sample1_ref, sample1_alt, sample2_ref, sample2_alt}
    record->SetFormatField(
        kDuplexConcordantCountsId,
        tumor_first ? vec<s32>{t_adc_ref, t_adc, n_adc_ref, n_adc} : vec<s32>{n_adc_ref, n_adc, t_adc_ref, t_adc});
    record->SetFormatField(
        kDuplexSimplexCountsId,
        tumor_first ? vec<s32>{t_ads_ref, t_ads, n_ads_ref, n_ads} : vec<s32>{n_ads_ref, n_ads, t_ads_ref, t_ads});
    record->SetFormatField(
        kDuplexDiscordantCountsId,
        tumor_first ? vec<s32>{t_add_ref, t_add, n_add_ref, n_add} : vec<s32>{n_add_ref, n_add, t_add_ref, t_add});
    record->SetFormatField(
        kDuplexLowbqCountsId,
        tumor_first ? vec<s32>{t_adl_ref, t_adl, n_adl_ref, n_adl} : vec<s32>{n_adl_ref, n_adl, t_adl_ref, t_adl});
  }
}

/**
 * @brief Copy a VCF record for a given allele.
 * @param original_record VCF record to be copied
 * @param new_header VCF header for the new VCF record
 * @return The copied VCF record
 */
io::VcfRecordPtr CopyRecord(const io::VcfRecordPtr& original_record,
                            const io::VcfHeaderPtr& new_header,
                            const VcfHeaderInfo& header_info) {
  const auto& record_copy = original_record->Clone(new_header);
  const std::string normal_gt = GenotypeToString(Genotype::kGT00);
  const std::string tumor_gt = GenotypeToString(Genotype::kGT01);
  record_copy->SetAlleles({original_record->Allele(0), original_record->Allele(1)});
  record_copy->SetGTField(normal_gt, header_info.normal_index);
  record_copy->SetGTField(tumor_gt, header_info.tumor_index);
  return record_copy;
}

/**
 * @brief Fail a VCF record and set all FORMAT fields to zero defaults.
 *
 * Without this, the cloned record retains raw Mutect2 FORMAT values which may have incorrect
 * cardinality (e.g. AD with a single value instead of one per allele as required by Number=R).
 * @param record VCF record
 * @param fail_id Filter ID to apply
 */
void FailSomaticTNRecord(const io::VcfRecordPtr& record,
                         const std::string_view fail_id,
                         const bool is_duplex_protocol) {
  record->AddFilter(fail_id);
  record->SetQuality(0);
  // Set all per-sample FORMAT fields to zero defaults (2 samples).
  // AD is Number=R (one value per allele per sample), so 4 values for a biallelic two-sample record.
  record->SetFormatField(kFieldAd, vec<s32>{0, 0, 0, 0});
  record->SetFormatField(kFieldAf, vec<f32>{0, 0});
  record->SetFormatField(kFieldDp, vec<s32>{0, 0});
  record->SetFormatField(kFieldGq, vec<s32>{0, 0});
  record->SetFormatField(kBaseqQualId, vec<f32>{0, 0});
  record->SetFormatField(kMapQualId, vec<f32>{0, 0});
  record->SetFormatField(kDistanceId, vec<f32>{0, 0});
  if (is_duplex_protocol) {
    record->SetFormatField(kDuplexConcordantCountsId, vec<s32>{0, 0, 0, 0});
    record->SetFormatField(kDuplexSimplexCountsId, vec<s32>{0, 0, 0, 0});
    record->SetFormatField(kDuplexDiscordantCountsId, vec<s32>{0, 0, 0, 0});
    record->SetFormatField(kDuplexLowbqCountsId, vec<s32>{0, 0, 0, 0});
  }
}
}  // namespace

/**
 * @brief Check pre-ML hard filters for a tumor-normal variant.
 *
 * Evaluates thresholds in priority order and returns the ID of the first failing filter,
 * or `kFilteringPassId` if all thresholds are satisfied.
 *
 * @param bam_feat         BAM feature tuple for the variant.
 * @param normalize_target Chromosome median depth targets used for DP-ratio filtering.
 * @return Filter ID of the first failing threshold, or `kFilteringPassId`.
 */
std::string_view TumorNormalProcessing::CheckHardFilters(const TumorNormalBamFeatureTuple& bam_feat,
                                                         const VcfFeature& vcf_feat,
                                                         const DepthTuple& normalize_target) const {
  if (bam_feat.tumor_var_feat.support < _min_tumor_support) {
    return kFilteringFailMinTumorSupportId;
  }
  if (bam_feat.normal_var_feat.support > _max_normal_support) {
    return kFilteringFailMaxNormalSupportId;
  }
  // Relax the tumor AF threshold for variants present in the catalog VCF (e.g. TCGA variants), as these
  // are more likely to be true positives. Falls back to the standard threshold when unset.
  const f32 min_tumor_af = vcf_feat.in_catalog ? _catalog_min_tumor_af.value_or(_min_tumor_af) : _min_tumor_af;
  if (bam_feat.tumor_var_feat.duplex_af < min_tumor_af) {
    return kFilteringFailMinTumorAfId;
  }
  if (_min_dp_ratio > 0 && normalize_target.tumor > 0 &&
      static_cast<f32>(vcf_feat.tumor_dp) / static_cast<f32>(normalize_target.tumor) < _min_dp_ratio) {
    return kFilteringFailMinDpRatioId;
  }
  return kFilteringPassId;
}

/**
 * @brief Run the somatic ML model, record SHAP values, and return the resulting filter and score.
 *
 * @param vid              Variant identifier.
 * @param bam_feat         BAM feature tuple for the variant.
 * @param vcf_feat         VCF feature values for the variant.
 * @param normalize_target Chromosome median depth targets.
 * @param calculator       Score calculator wrapping the ML model.
 * @param shap_value_rows  Output accumulator for SHAP value rows.
 * @return Pair of {filter ID, prediction score}.
 */
std::pair<std::string_view, PredictionScore> TumorNormalProcessing::ScoreVariantML(
    const VariantId& vid,
    const TumorNormalBamFeatureTuple& bam_feat,
    const VcfFeature& vcf_feat,
    const DepthTuple& normalize_target,
    const ScoreCalculator& calculator,
    vec<vec<std::string>>& shap_value_rows) const {
  const auto& scoring_target = _normalize_scoring_features ? normalize_target : kZeroDepthTuple;
  const auto& feature_vec = GetFeatureVec(_scoring_cols, vid, bam_feat, vcf_feat, scoring_target);
  const auto score = calculator.CalculateScore(feature_vec);
  if (!score.shap_values.empty()) {
    shap_value_rows.emplace_back(AssembleShapValueTsvRow(vid, score));
  }
  const bool is_snv = vid.type == VariantType::kSNV;
  auto ml_threshold = is_snv ? _snv_min_ml_score : _indel_min_ml_score;
  // Relax the ML score threshold for variants present in the catalog VCF (e.g. TCGA variants), as these
  // are more likely to be true positives. Falls back to the standard threshold when unset.
  if (vcf_feat.in_catalog) {
    const auto& catalog_threshold = is_snv ? _catalog_snv_min_ml_score : _catalog_indel_min_ml_score;
    ml_threshold = catalog_threshold.value_or(ml_threshold);
  }
  const auto* const filter = score.probability < ml_threshold ? kFilteringMLScoreId : kFilteringPassId;
  return {filter, score};
}

void TumorNormalProcessing::ProcessRecord(const io::VcfRecordPtr& record,
                                          vec<io::VcfRecordPtr>& out_records,
                                          const ScoreCalculator& somatic_calculator,
                                          const DepthTuple& normalize_target,
                                          const VcfHeaderInfo& header_info,
                                          vec<vec<std::string>>& shap_value_rows) const {
  using enum VariantType;
  const auto& ref = record->Allele(0);
  const auto& alt = record->Allele(1);
  auto record_copy = CopyRecord(record, _hdr, header_info);

  if (!ContainsOnlyACTG(ref) || !ContainsOnlyACTG(alt)) {
    FailSomaticTNRecord(record_copy, kFilteringNonAcgtRefAltId, _is_duplex_protocol);
    out_records.emplace_back(record_copy);
    return;
  }

  const auto& chrom = record->Chromosome();
  const auto pos = record->Position();
  const auto& [ref_trimmed, alt_trimmed] = TrimVariant(ref, alt);
  const VariantId vid(chrom, ToUnsigned(pos), ref_trimmed, alt_trimmed);

  std::optional<VcfFeature> vcf_feat;
  const auto vcf_features_it = _vcf_features.find(vid);
  if (vcf_features_it != _vcf_features.end()) {
    vcf_feat = vcf_features_it->second;
  }

  const std::optional<TumorNormalBamFeatureTuple> bam_feat =
      vcf_feat.has_value() ? _bam_features.GetTumorNormalBamFeatureTuple(vid) : std::nullopt;

  record_copy->SetAlleles({ref_trimmed, alt_trimmed});

  if (vid.type != kSNV && vid.type != kUnknown) {
    const auto ref_len = ref_trimmed.size();
    const auto alt_len = alt_trimmed.size();
    const auto indel_size = static_cast<u32>(std::max(ref_len, alt_len) - std::min(ref_len, alt_len));
    if (indel_size > _max_indel_size) {
      FailSomaticTNRecord(record_copy, kFilteringFailMaxIndelSizeId, _is_duplex_protocol);
      out_records.emplace_back(record_copy);
      return;
    }
  }

  if (!vcf_feat.has_value() || !bam_feat.has_value()) {
    FailSomaticTNRecord(record_copy, kFilteringMissingFeatureId, _is_duplex_protocol);
    out_records.emplace_back(record_copy);
    return;
  }

  PredictionScore score;
  auto filter = CheckHardFilters(*bam_feat, *vcf_feat, normalize_target);
  if (filter == kFilteringPassId) {
    std::tie(filter, score) =
        ScoreVariantML(vid, *bam_feat, *vcf_feat, normalize_target, somatic_calculator, shap_value_rows);
  }
  SetSomaticTNValues(record_copy,
                     vid,
                     *bam_feat,
                     *vcf_feat,
                     header_info.tumor_index < header_info.normal_index,
                     score,
                     _is_duplex_protocol);
  record_copy->SetFilter(filter);
  out_records.emplace_back(record_copy);
}

}  // namespace xoos::svc
