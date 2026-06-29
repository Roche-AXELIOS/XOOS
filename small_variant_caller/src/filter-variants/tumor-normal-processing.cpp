#include "tumor-normal-processing.h"

#include "core/filtering.h"
#include "core/variant-feature-extraction.h"
#include "core/vcf-fields.h"
#include "util/seq-util.h"

namespace xoos::svc {

TumorNormalProcessing::TumorNormalProcessing(const GlobalContext& global_ctx,
                                             const u64 prev_pos,
                                             VarIdToVcfFeatures vcf_features,
                                             BamRegionFeatureCollection bam_features)
    : _prev_pos(prev_pos),
      _vcf_features(std::move(vcf_features)),
      _bam_features(std::move(bam_features)),
      _hdr(global_ctx.hdr),
      _somatic_ml_snv_threshold(global_ctx.somatic_tn_snv_ml_threshold),
      _somatic_ml_indel_threshold(global_ctx.somatic_tn_indel_ml_threshold),
      _tumor_support_threshold(global_ctx.tumor_support_threshold),
      _scoring_cols(global_ctx.model_config.scoring_cols) {
}

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
static void SetSomaticTNValues(const io::VcfRecordPtr& record,
                               const VariantId& variant_id,
                               const TumorNormalBamFeatureTuple& bam_feat,
                               const VcfFeature& vcf_feat,
                               const bool tumor_first,
                               const f64 pred_score) {
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
  record->SetFormatField(
      kMachineLearningId,
      tumor_first ? vec<f32>{static_cast<f32>(pred_score), 0} : vec<f32>{0, static_cast<f32>(pred_score)});
  // Set variant quality based on the ML score alone because this is a binary classification.
  record->SetQuality(static_cast<f32>(CalculatePhredQualityScore(pred_score)));
  record->SetInfoField<f32>(kRefBQId, {static_cast<f32>(bam_feat.ref_feat.baseq_mean)});
  record->SetInfoField<f32>(kRefMQId, {static_cast<f32>(bam_feat.ref_feat.mapq_mean)});
  record->SetInfoField<f32>(kAltBQId, {static_cast<f32>(bam_feat.var_feat.baseq_mean)});
  record->SetInfoField<f32>(kAltMQId, {static_cast<f32>(bam_feat.var_feat.mapq_mean)});
  record->SetInfoField<f32>(kFieldNalod, {vcf_feat.nalod});
  record->SetInfoField<f32>(kFieldNlod, {vcf_feat.nlod});
  record->SetInfoField<f32>(kFieldTlod, {vcf_feat.tlod});
  record->SetInfoField<s32>(kFieldMpos, {static_cast<s32>(vcf_feat.mpos)});
  record->SetInfoField<s32>(kSubtypeId, {SubstIndex(variant_id.ref, variant_id.alt)});
  record->SetInfoField<std::string>(kContextId, {bam_feat.var_feat.context});
  record->SetInfoField<f32>(kFieldPopaf, {vcf_feat.popaf});
}

/**
 * @brief Copy a VCF record for a given allele.
 * @param original_record VCF record to be copied
 * @param new_header VCF header for the new VCF record
 * @return The copied VCF record
 */
static io::VcfRecordPtr CopyRecord(const io::VcfRecordPtr& original_record,
                                   const io::VcfHeaderPtr& new_header,
                                   const VcfHeaderInfo& header_info) {
  const auto& record_copy = original_record->Clone(new_header);
  const std::string normal_gt = GenotypeToString(Genotype::kGT00);
  const std::string tumor_gt = GenotypeToString(Genotype::kGT01);
  record_copy->SetGTField(normal_gt, header_info.normal_index);
  record_copy->SetGTField(tumor_gt, header_info.tumor_index);
  record_copy->SetAlleles({original_record->Allele(0), original_record->Allele(0)});
  return record_copy;
}

/**
 * @brief Fail a VCF record and update relevant fields.
 * @param record VCF record
 */
static void FailSomaticTNRecord(const io::VcfRecordPtr& record, const std::string& fail_id) {
  record->AddFilter(fail_id);
  // Set VCF GQ and QUAL to 0
  record->SetFormatField(kFieldGq, vec<s32>{0, 0});
  record->SetQuality(0);
}

void TumorNormalProcessing::ProcessRecord(const io::VcfRecordPtr& record,
                                          vec<io::VcfRecordPtr>& out_records,
                                          const ScoreCalculator& somatic_calculator,
                                          const DepthTuple& normalize_target,
                                          const VcfHeaderInfo& header_info) const {
  using enum VariantType;
  const auto& ref = record->Allele(0);
  const bool ref_acgt_only = ContainsOnlyACTG(ref);
  const auto& alt = record->Allele(1);
  auto record_copy = CopyRecord(record, _hdr, header_info);
  if (!ref_acgt_only) {
    // REF/ALT contains non-ACGT; copy the record and set FILTER to FAIL
    FailSomaticTNRecord(record_copy, kFilteringFailId);
    out_records.emplace_back(record_copy);
    return;
  }
  if (alt == "*") {
    out_records.emplace_back(record_copy);
  } else {
    const auto& chrom = record->Chromosome();
    const auto pos = record->Position();
    const auto& [ref_trimmed, alt_trimmed] = TrimVariant(ref, alt);
    const VariantId vid(chrom, pos, ref_trimmed, alt_trimmed);

    std::optional<VcfFeature> vcf_feat;
    const auto vcf_features_it = _vcf_features.find(vid);
    if (vcf_features_it != _vcf_features.end()) {
      vcf_feat = vcf_features_it->second;
    }

    std::optional<TumorNormalBamFeatureTuple> bam_feat;
    if (vcf_feat.has_value()) {
      bam_feat = _bam_features.GetTumorNormalBamFeatureTuple(vid);
    }

    record_copy->SetAlleles({ref_trimmed, alt_trimmed});
    if (vcf_feat.has_value() && bam_feat.has_value()) {
      // Check feature specific values first, then apply ML only if needed.
      f64 score = 0;
      auto filter_to_set = kFilteringPassId;
      if (bam_feat->tumor_var_feat.support < _tumor_support_threshold) {
        filter_to_set = kFilteringCountsId;
      } else if (bam_feat->normal_var_feat.support > 2) {
        filter_to_set = kFilteringFailSomaticTNNormalADId;
      } else if (bam_feat->tumor_var_feat.duplex_af < 0.05) {
        filter_to_set = kFilteringAFId;
      } else {
        const auto& feature_vec =
            GetFeatureVec(_scoring_cols, vid, bam_feat.value(), vcf_feat.value(), normalize_target);
        score = somatic_calculator.CalculateScore(feature_vec);
        const auto ml_threshold_to_use = vid.type == kSNV ? _somatic_ml_snv_threshold : _somatic_ml_indel_threshold;
        if (score < ml_threshold_to_use) {
          filter_to_set = kFilteringFailSomaticTNMLId;
        }
      }
      SetSomaticTNValues(record_copy,
                         vid,
                         bam_feat.value(),
                         vcf_feat.value(),
                         header_info.tumor_index < header_info.normal_index,
                         score);
      record_copy->SetFilter(filter_to_set);
      out_records.emplace_back(record_copy);
    } else {
      // Either BAM or VCF features not found; copy the record and set FILTER to FAIL
      FailSomaticTNRecord(record_copy, kFilteringMissingFeatureId);
      out_records.emplace_back(record_copy);
    }
  }
}
}  // namespace xoos::svc
