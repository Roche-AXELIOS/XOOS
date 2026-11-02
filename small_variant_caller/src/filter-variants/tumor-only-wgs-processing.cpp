#include "tumor-only-wgs-processing.h"

#include <cmath>

#include "core/feature-normalization.h"
#include "core/filtering.h"
#include "core/sequencing-protocol.h"
#include "core/variant-feature-extraction.h"
#include "core/vcf-fields.h"
#include "shap-value-tsv.h"
#include "util/seq-util.h"

namespace xoos::svc {

// Genotype quality assigned to passing tumor-only-wgs calls. The single-sample somatic model emits a binary
// classification, so a fixed high GQ is used rather than a per-genotype likelihood.
constexpr s32 kTumorOnlyWgsGq{100};

TumorOnlyWgsProcessing::TumorOnlyWgsProcessing(const GlobalContext& global_ctx,
                                               VarIdToVcfFeatures vcf_features,
                                               BamRegionFeatureCollection bam_features)
    : _vcf_features(std::move(vcf_features)),
      _bam_features(std::move(bam_features)),
      _hdr(global_ctx.hdr),
      _snv_min_ml_score(global_ctx.snv_min_ml_score),
      _indel_min_ml_score(global_ctx.indel_min_ml_score),
      _min_support(global_ctx.min_support),
      _min_af(global_ctx.min_af),
      _max_indel_size(global_ctx.max_indel_size),
      _normalize_scoring_features(global_ctx.normalize_scoring_features),
      _scoring_cols(global_ctx.model_config.scoring_cols),
      _is_duplex_protocol(IsDuplexProtocol(global_ctx.bam_feat_params.sequencing_protocol)) {
}

namespace {

/**
 * @brief Set INFO and FORMAT fields on a tumor-only-wgs VCF record before output.
 *
 * Sets AD, AF, DP, BQ, MQ, DT, GQ, PRED_ML, REF/ALT BQ/MQ, TLOD, MPOS, SUBTYPE, CONTEXT, POPAF, and variant quality.
 * When the sequencing protocol is duplex, also sets ADC, ADS, ADD, and ADL.
 */
void SetSomaticTumorOnlyWgsValues(const io::VcfRecordPtr& record,
                                  const VariantId& variant_id,
                                  const BamFeatureTuple& bam_feat,
                                  const VcfFeature& vcf_feat,
                                  const PredictionScore& pred_score,
                                  const bool is_duplex_protocol) {
  const auto support = static_cast<s32>(bam_feat.var_feat.support);
  const auto ref_support = static_cast<s32>(bam_feat.ref_feat.support);
  record->SetFormatField(kFieldAd, vec<s32>{ref_support, support});
  record->SetFormatField(kFieldAf, vec<f32>{static_cast<f32>(bam_feat.var_feat.duplex_af)});
  record->SetFormatField(kFieldDp, vec<s32>{support + ref_support});
  record->SetFormatField(kBaseqQualId, vec<f32>{static_cast<f32>(bam_feat.var_feat.baseq_mean)});
  record->SetFormatField(kMapQualId, vec<f32>{static_cast<f32>(bam_feat.var_feat.mapq_mean)});
  record->SetFormatField(kDistanceId, vec<f32>{static_cast<f32>(bam_feat.var_feat.distance_mean)});
  record->SetFormatField(kFieldGq, vec<s32>{kTumorOnlyWgsGq});
  record->SetFormatField(kMachineLearningId, vec<f32>{static_cast<f32>(pred_score.probability)});
  if (is_duplex_protocol) {
    record->SetFormatField(kDuplexConcordantCountsId,
                           vec<s32>{static_cast<s32>(bam_feat.ref_feat.duplex_concordant),
                                    static_cast<s32>(bam_feat.var_feat.duplex_concordant)});
    record->SetFormatField(kDuplexSimplexCountsId,
                           vec<s32>{static_cast<s32>(bam_feat.ref_feat.duplex_simplex),
                                    static_cast<s32>(bam_feat.var_feat.duplex_simplex)});
    record->SetFormatField(kDuplexDiscordantCountsId,
                           vec<s32>{static_cast<s32>(bam_feat.ref_feat.duplex_discordant),
                                    static_cast<s32>(bam_feat.var_feat.duplex_discordant)});
    record->SetFormatField(
        kDuplexLowbqCountsId,
        vec<s32>{static_cast<s32>(std::round(bam_feat.ref_feat.duplex_lowbq * kDuplexLowbqToCountFactor)),
                 static_cast<s32>(std::round(bam_feat.var_feat.duplex_lowbq * kDuplexLowbqToCountFactor))});
  }
  record->SetQuality(pred_score.variant_quality);
  record->SetInfoField(kRefBQId, vec<f32>{static_cast<f32>(bam_feat.ref_feat.baseq_mean)});
  record->SetInfoField(kRefMQId, vec<f32>{static_cast<f32>(bam_feat.ref_feat.mapq_mean)});
  record->SetInfoField(kAltBQId, vec<f32>{static_cast<f32>(bam_feat.var_feat.baseq_mean)});
  record->SetInfoField(kAltMQId, vec<f32>{static_cast<f32>(bam_feat.var_feat.mapq_mean)});
  record->SetInfoField(kFieldTlod, vec<f32>{vcf_feat.tlod});
  record->SetInfoField(kFieldMpos, vec<s32>{static_cast<s32>(vcf_feat.mpos)});
  record->SetInfoField(kSubtypeId, vec<s32>{SubstIndex(variant_id.ref, variant_id.alt)});
  record->SetInfoField(kContextId, vec<std::string>{bam_feat.var_feat.context});
  record->SetInfoField(kFieldPopaf, vec<f32>{vcf_feat.popaf});
}

/**
 * @brief Fail a tumor-only-wgs VCF record and set FORMAT fields to zero defaults.
 */
void FailTumorOnlyWgsRecord(const io::VcfRecordPtr& record,
                            const std::string_view fail_id,
                            const bool is_duplex_protocol) {
  record->AddFilter(fail_id);
  record->SetQuality(0);
  record->SetFormatField(kFieldAd, vec<s32>{0, 0});
  record->SetFormatField(kFieldAf, vec<f32>{0});
  record->SetFormatField(kFieldDp, vec<s32>{0});
  record->SetFormatField(kFieldGq, vec<s32>{0});
  record->SetFormatField(kBaseqQualId, vec<f32>{0});
  record->SetFormatField(kMapQualId, vec<f32>{0});
  record->SetFormatField(kDistanceId, vec<f32>{0});
  if (is_duplex_protocol) {
    record->SetFormatField(kDuplexConcordantCountsId, vec<s32>{0, 0});
    record->SetFormatField(kDuplexSimplexCountsId, vec<s32>{0, 0});
    record->SetFormatField(kDuplexDiscordantCountsId, vec<s32>{0, 0});
    record->SetFormatField(kDuplexLowbqCountsId, vec<s32>{0, 0});
  }
}

}  // namespace

std::string_view TumorOnlyWgsProcessing::CheckHardFilters(const BamFeatureTuple& bam_feat) const {
  if (bam_feat.var_feat.support < _min_support) {
    return kFilteringFailMinSupportId;
  }
  if (bam_feat.var_feat.duplex_af < _min_af) {
    return kFilteringFailMinAfId;
  }
  return kFilteringPassId;
}

std::pair<std::string_view, PredictionScore> TumorOnlyWgsProcessing::ScoreVariantML(
    const VariantId& vid,
    const BamFeatureTuple& bam_feat,
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
  const auto ml_threshold = vid.type == VariantType::kSNV ? _snv_min_ml_score : _indel_min_ml_score;
  const auto* const filter = score.probability < ml_threshold ? kFilteringMLScoreId : kFilteringPassId;
  return {filter, score};
}

void TumorOnlyWgsProcessing::ProcessRecord(const io::VcfRecordPtr& record,
                                           vec<io::VcfRecordPtr>& out_records,
                                           const ScoreCalculator& somatic_calculator,
                                           const DepthTuple& normalize_target,
                                           vec<vec<std::string>>& shap_value_rows) const {
  using enum VariantType;
  const auto& ref = record->Allele(0);
  const auto& alt = record->Allele(1);
  auto record_copy = record->Clone(_hdr);
  record_copy->SetAlleles({ref, alt});

  if (!ContainsOnlyACTG(ref) || !ContainsOnlyACTG(alt)) {
    FailTumorOnlyWgsRecord(record_copy, kFilteringNonAcgtRefAltId, _is_duplex_protocol);
    out_records.emplace_back(record_copy);
    return;
  }

  const auto& chrom = record->Chromosome();
  const auto pos = static_cast<u64>(record->Position());
  const auto& [ref_trimmed, alt_trimmed] = TrimVariant(ref, alt);
  const VariantId vid(chrom, pos, ref_trimmed, alt_trimmed);

  const auto vcf_features_it = _vcf_features.find(vid);
  const bool has_vcf_feat = vcf_features_it != _vcf_features.end();

  const std::optional<BamFeatureTuple> bam_feat = has_vcf_feat ? _bam_features.GetBamFeatureTuple(vid) : std::nullopt;

  record_copy->SetAlleles({ref_trimmed, alt_trimmed});

  if (vid.type != kSNV && vid.type != kUnknown) {
    const auto ref_len = ref_trimmed.size();
    const auto alt_len = alt_trimmed.size();
    const auto indel_size = static_cast<u32>(std::max(ref_len, alt_len) - std::min(ref_len, alt_len));
    if (indel_size > _max_indel_size) {
      FailTumorOnlyWgsRecord(record_copy, kFilteringFailMaxIndelSizeId, _is_duplex_protocol);
      out_records.emplace_back(record_copy);
      return;
    }
  }

  if (!has_vcf_feat || !bam_feat.has_value()) {
    FailTumorOnlyWgsRecord(record_copy, kFilteringMissingFeatureId, _is_duplex_protocol);
    out_records.emplace_back(record_copy);
    return;
  }

  const VcfFeature& vcf_feat = vcf_features_it->second;

  PredictionScore score;
  auto filter = CheckHardFilters(*bam_feat);
  if (filter == kFilteringPassId) {
    std::tie(filter, score) =
        ScoreVariantML(vid, *bam_feat, vcf_feat, normalize_target, somatic_calculator, shap_value_rows);
  }
  SetSomaticTumorOnlyWgsValues(record_copy, vid, *bam_feat, vcf_feat, score, _is_duplex_protocol);
  record_copy->SetFilter(filter);
  out_records.emplace_back(record_copy);
}

}  // namespace xoos::svc
