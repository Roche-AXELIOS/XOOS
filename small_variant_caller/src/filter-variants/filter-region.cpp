#include "filter-region.h"

#include <algorithm>
#include <memory>

#include <taskflow/taskflow.hpp>

#include <xoos/error/error.h>
#include <xoos/io/fasta-reader.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>
#include <xoos/util/container-functions.h>
#include <xoos/util/string-functions.h>

#include "compute-vcf-features/compute-vcf-features.h"
#include "core/filtering.h"
#include "core/variant-feature-extraction.h"
#include "core/vcf-fields.h"
#include "germline-tagging.h"
#include "reconcile-germline.h"
#include "tumor-normal-processing.h"
#include "tumor-only-te-processing.h"
#include "util/log-util.h"
#include "util/parallel-compute-utils.h"
#include "util/seq-util.h"
#include "util/vcf-util.h"

namespace xoos::svc {

/**
 * @brief Update a VCF record for the given features within the germline workflow.
 * @param record VCF record to be updated
 * @param bam_feat BAM features
 * @param vcf_feat VCF features
 * @param ref_feat Reference features
 */
static void UpdateGermlineRecordWithFeatures(const io::VcfRecordPtr& record,
                                             const BamFeatureTuple& bam_feat,
                                             const VcfFeature& vcf_feat) {
  const auto& bam_ref_feat = bam_feat.ref_feat;
  const auto& bam_var_feat = bam_feat.var_feat;
  record->SetFormatField<s32>(kGermlineMLId, {1});
  // "lowbq" features are incremented by 0.5 per read, so they are not always an integer
  // To be consistent for both AD and DP fields, use `std::round()` to convert to nearest integer
  record->SetFormatField<s32>(kAlleleDepthId,
                              {static_cast<s32>(std::round(bam_ref_feat.support + bam_ref_feat.duplex_lowbq)),
                               static_cast<s32>(std::round(bam_var_feat.support + bam_var_feat.duplex_lowbq))});
  record->SetFormatField<s32>(kFieldDp, {static_cast<s32>(std::round(bam_ref_feat.duplex_dp))});
  record->SetFormatField<f32>(kGermlineGnomadAFId, {vcf_feat.popaf});
  record->SetFormatField<f32>(kGermlineRefAvgMapqId, {static_cast<f32>(bam_ref_feat.mapq_mean)});
  record->SetFormatField<f32>(kGermlineAltAvgMapqId, {static_cast<f32>(bam_var_feat.mapq_mean)});
  record->SetFormatField<f32>(kGermlineRefAvgDistId, {static_cast<f32>(bam_ref_feat.distance_mean)});
  record->SetFormatField<f32>(kGermlineAltAvgDistId, {static_cast<f32>(bam_var_feat.distance_mean)});
  record->SetFormatField<f32>(kGermlineDensity100BPId, {static_cast<f32>(vcf_feat.variant_density)});
}

/**
 * @brief Update a copied VCF record for the germline workflow. Copies the original DP, AD, GT, and ALT fields from the
 * original VCF record to the copied record using GATK-specific field IDs.
 * @param original_record Original VCF record
 * @param record_copy Copied VCF record
 */
static void UpdateGermlineRecordCopy(const io::VcfRecordPtr& original_record, const io::VcfRecordPtr& record_copy) {
  const auto& dp_values = original_record->GetFormatFieldNoCheck<s32>(kFieldDp);
  const auto& ad_values = original_record->GetFormatFieldNoCheck<s32>(kFieldAd);
  // note that AD and ALT values from multi-allelic records are not split here
  const std::string& gt_value = original_record->GetGTField();
  vec<std::string> alt_vec;
  for (s32 i = 1; i < original_record->NumAlleles(); ++i) {
    alt_vec.emplace_back(original_record->Allele(i));
  }
  const std::string alt_string = string::Join(alt_vec, ",");
  record_copy->SetFormatField<s32>(kGermlineGATKDPId, dp_values);
  record_copy->SetFormatField<s32>(kGermlineGATKADId, ad_values);
  record_copy->SetFormatField<std::string>(kGermlineGATKGTId, {gt_value});
  record_copy->SetFormatField<std::string>(kGermlineGATKAltId, {alt_string});
}

void FilterRegionClass::FilterGermlineRecord(const io::VcfRecordPtr& record,
                                             const VarIdToVcfFeatures& vcf_features,
                                             const BamRegionFeatureCollection& bam_features,
                                             const DepthTuple& normalize_target,
                                             vec<io::VcfRecordPtr>& out_records) {
  using enum VariantType;
  using enum Genotype;
  const auto& chrom = record->Chromosome();
  const auto pos = static_cast<u64>(record->Position());
  const auto snv_calculator = [&snv_calc = _worker_ctx->calculators.at(0)](const auto& feature_vec) {
    return snv_calc.CalculateScoreGermline(feature_vec);
  };
  const auto& indel_calculator = [&indel_calc = _worker_ctx->calculators.at(1)](const auto& feature_vec) {
    return indel_calc.CalculateScoreGermline(feature_vec);
  };
  const auto& snv_scoring_cols = _global_ctx.model_config.snv_scoring_cols;
  const auto& indel_scoring_cols = _global_ctx.model_config.indel_scoring_cols;
  const bool ref_acgt_only = ContainsOnlyACTG(record->Allele(0));

  // a VCF record can have >=1 ALT allele
  for (auto& record_copy :
       SplitMultiAllelicRecord(record, _hdr, _global_ctx.vcf_info_metadata, _global_ctx.vcf_fmt_metadata, true)) {
    UpdateGermlineRecordCopy(record, record_copy);
    if (!ref_acgt_only) {
      // REF/ALT contains non-ACGT; copy the record and set FILTER to FAIL
      FailRecord(record_copy, kFilteringNonAcgtRefAltId, kGT00);
      out_records.emplace_back(record_copy);
      continue;
    }
    const auto& ref = record_copy->Allele(0);
    const auto& alt = record_copy->Allele(1);
    if (alt == "*") {
      _tmp_alt_wildcard_records.emplace_back(record_copy);
      continue;
    }
    const VariantId vid(chrom, pos, ref, alt);

    std::optional<VcfFeature> vcf_feat;
    auto vcf_features_it = vcf_features.find(vid);
    if (vcf_features_it != vcf_features.end() && vcf_features_it->second.genotype != kGTNA) {
      vcf_feat = vcf_features_it->second;
    }

    std::optional<BamFeatureTuple> bam_feat;
    if (vcf_feat.has_value()) {
      bam_feat = bam_features.GetBamFeatureTuple(vid);
    }

    if (vcf_feat.has_value() && bam_feat.has_value()) {
      _tmp_vids.emplace_back(vid);
      _tmp_records.emplace_back(record_copy);

      // update the VCF record with ML features
      UpdateGermlineRecordWithFeatures(record_copy, bam_feat.value(), vcf_feat.value());

      // classify the preliminary genotype and append it to the temporary vector
      if (vid.type == kSNV) {
        const auto& feature_vec =
            GetFeatureVec(snv_scoring_cols, vid, bam_feat.value(), vcf_feat.value(), normalize_target);
        _tmp_genotypes.emplace_back(snv_calculator(feature_vec));
      } else {
        const auto& feature_vec =
            GetFeatureVec(indel_scoring_cols, vid, bam_feat.value(), vcf_feat.value(), normalize_target);
        _tmp_genotypes.emplace_back(indel_calculator(feature_vec));
      }
    } else {
      // features not found; copy the record and set FILTER to FAIL
      FailRecord(record_copy, kFilteringMissingFeatureId, kGT00);
      out_records.emplace_back(record_copy);
    }
  }
}

void FilterRegionClass::FilterGermlineTaggingRecord(const io::VcfRecordPtr& record,
                                                    const VarIdToVcfFeatures& vcf_features,
                                                    const BamRegionFeatureCollection& bam_features,
                                                    const DepthTuple& normalize_target) {
  using enum VariantType;
  using enum Genotype;
  const auto& chrom = record->Chromosome();
  const auto pos = static_cast<u64>(record->Position());
  const auto snv_min_score = _global_ctx.somatic_tn_snv_ml_threshold;
  const auto indel_min_score = _global_ctx.somatic_tn_indel_ml_threshold;
  const auto calculator = [&calc = _worker_ctx->calculators.front()](const auto& feature_vec, const f64 min_score) {
    return calc.CalculateScoreGermline(feature_vec, min_score);
  };
  const auto& scoring_cols = _global_ctx.model_config.scoring_cols;
  const bool ref_acgt_only = ContainsOnlyACTG(record->Allele(0));
  // a VCF record can have >=1 ALT allele
  for (auto& record_copy :
       SplitMultiAllelicRecord(record, _hdr, _global_ctx.vcf_info_metadata, _global_ctx.vcf_fmt_metadata, true)) {
    const auto& ref = record_copy->Allele(0);
    const auto& alt = record_copy->Allele(1);
    if (!ref_acgt_only) {
      _tmp_vids.emplace_back(chrom, pos, ref, alt);
      _tmp_records.emplace_back(record_copy);
      _tmp_genotypes.emplace_back(kGT00);
      continue;
    }
    if (alt == "*") {
      continue;
    }
    const VariantId vid(chrom, pos, ref, alt);

    std::optional<VcfFeature> vcf_feat;
    const auto vcf_features_it = vcf_features.find(vid);
    if (vcf_features_it != vcf_features.end() && vcf_features_it->second.genotype != kGTNA) {
      vcf_feat = vcf_features_it->second;
    }

    std::optional<TumorNormalBamFeatureTuple> bam_feat;
    if (vcf_feat.has_value()) {
      bam_feat = bam_features.GetTumorNormalBamFeatureTuple(vid);
    }

    if (vcf_feat.has_value() && bam_feat.has_value()) {
      _tmp_vids.emplace_back(vid);
      _tmp_records.emplace_back(record_copy);
      // classify the preliminary genotype and append it to the temporary vector
      const auto& feature_vec = GetFeatureVec(scoring_cols, vid, bam_feat.value(), vcf_feat.value(), normalize_target);
      const auto min_score = vid.type == kSNV ? snv_min_score : indel_min_score;
      _tmp_genotypes.emplace_back(calculator(feature_vec, min_score));
    } else {
      _tmp_vids.emplace_back(vid);
      // features not found; copy the record and set FILTER to FAIL
      _tmp_records.emplace_back(record_copy);
      _tmp_genotypes.emplace_back(kGT00);
    }
  }
}

void FilterRegionClass::ReconcileGermlineFeatures(const TargetRegion& region,
                                                  const u64 pos,
                                                  std::optional<u64>& prev_pos,
                                                  std::optional<u64>& prev_pass_ref_max_pos,
                                                  vec<io::VcfRecordPtr>& out_records) {
  if (!_tmp_records.empty() && (prev_pos.has_value() && prev_pos.value() != pos)) {
    // Reconcile predicted genotypes at the previous position and add updated VCF records to the output vector
    if (!prev_pass_ref_max_pos.has_value() || prev_pass_ref_max_pos.value() < _tmp_vids.front().pos) {
      // Wildcard ALT alleles are valid only if an upstream variant overlaps this position
      _tmp_alt_wildcard_records = {};
    }

    std::optional<u64> pass_ref_max_pos;
    const u64 prev_pos_end = prev_pos.value() + _tmp_vids.front().ref.size();
    const bool is_haploid = IsHaploid(region, prev_pos.value(), prev_pos_end);
    if (is_haploid) {
      // haploid regions:
      // - biological male's chr X/Y non-PARs
      ReconcileGermlineHaploidRecords(_tmp_vids, _tmp_records, _tmp_genotypes, out_records);
    } else {
      // diploid regions:
      // - all autosomes
      // - biological male's chr X/Y PARs
      // - biological female's chr X
      pass_ref_max_pos = ReconcileGermlineDiploidRecords(_tmp_vids,
                                                         _tmp_records,
                                                         _tmp_alt_wildcard_records,
                                                         _tmp_genotypes,
                                                         out_records,
                                                         _global_ctx.vcf_info_metadata,
                                                         _global_ctx.vcf_fmt_metadata);
    }
    if (pass_ref_max_pos.has_value()) {
      if (prev_pass_ref_max_pos.has_value()) {
        prev_pass_ref_max_pos = std::max(prev_pass_ref_max_pos.value(), pass_ref_max_pos.value());
      } else {
        prev_pass_ref_max_pos = pass_ref_max_pos.value();
      }
    }

    // reset temp storage of variant records and features
    _tmp_vids = {};
    _tmp_records = {};
    _tmp_genotypes = {};
    _tmp_alt_wildcard_records = {};
  }
  prev_pos = pos;
}

void FilterRegionClass::ReconcileGermlineTaggingFeatures(const TargetRegion& region,
                                                         const u64 pos,
                                                         std::optional<u64>& prev_pos,
                                                         vec<io::VcfRecordPtr>& out_records) {
  if (!_tmp_records.empty() && (prev_pos.has_value() && prev_pos.value() != pos)) {
    // Reconcile predicted genotypes at the previous position and add updated VCF records to the output vector
    const u64 prev_pos_end = prev_pos.value() + _tmp_vids.front().ref.size();
    const bool is_haploid = IsHaploid(region, prev_pos.value(), prev_pos_end);
    ReconcileGermlineTaggingRecords(
        _tmp_records, _tmp_genotypes, is_haploid, _global_ctx.vcf_normal_index.value(), out_records);

    // reset temp storage of variant records and features
    _tmp_vids = {};
    _tmp_records = {};
    _tmp_genotypes = {};
    _tmp_alt_wildcard_records = {};
  }
  prev_pos = pos;
}

bool FilterRegionClass::IsHaploid(const TargetRegion& region, const u64 prev_pos, const u64 prev_pos_end) const {
  const bool is_male = _global_ctx.sex == sex_predict::Sex::kMale;
  if (!is_male) {
    // female has no haploid regions
    return false;
  }
  const bool is_chr_x = region.chrom == _global_ctx.chr_x_name;
  const bool is_chr_y = region.chrom == _global_ctx.chr_y_name;
  if (!is_chr_x && !is_chr_y) {
    // only chr X and chr Y in male can have haploid regions
    return false;
  }
  const bool at_chr_x_par = is_chr_x && IntervalOverlap(_global_ctx.chr_x_par, region.start, region.end);
  const bool at_chr_y_par = is_chr_y && IntervalOverlap(_global_ctx.chr_y_par, region.start, region.end);
  const bool x_overlap =
      (is_chr_x && (!at_chr_x_par || !IntervalOverlap(_global_ctx.chr_x_par, prev_pos, prev_pos_end)));
  const bool y_overlap =
      (is_chr_y && (!at_chr_y_par || !IntervalOverlap(_global_ctx.chr_y_par, prev_pos, prev_pos_end)));
  return x_overlap || y_overlap;
}

void FilterRegionClass::UpdateWorkerCtx(std::unique_ptr<WorkerContext>& worker_ctx) {
  _worker_ctx = std::move(worker_ctx);
}

void FilterRegionClass::FilterTumorNormalRegion(const xoos::svc::TargetRegion& region,
                                                xoos::svc::FlowContext& flow_ctx) {
  // This function is designed to filter a given region of a GATK Mutect2 output VCF in a paired somatic tumor normal
  // analysis. It is designed to use a single thread and be run in parallel for multiple regions simultaneously as part
  // of a region based end to end filtering workflow. If there are variants that need to be filtered in the specified
  // region of the genome, BAM and VCF features for the region are first computed before variants are read in and
  // filtered using a LightGBM ML model to determine if they are indeed somatic variants. Variants that are deemed to
  // not be somatic from this model are then passed to a second LightGBM ML model that attempts to determine if they
  // represent germline variation between the sample and the reference genome or noise. Finally, all variants are
  // written to the output buffer.
  using enum VariantType;
  // Load values for worker and global contexts
  auto& reader = _worker_ctx->filter_variants_reader;

  const auto& somatic_calculator = _worker_ctx->calculators.front();

  auto& out_records = flow_ctx.out_records;
  vec<io::VcfRecordPtr> fail_records;

  const auto& hdr = _global_ctx.hdr;
  const auto& bed_regions = _global_ctx.bed_regions;
  const auto& interest_regions = _global_ctx.interest_regions;
  const auto& normalize_targets = _global_ctx.normalize_targets;

  if (!reader.SetRegion(region.chrom, region.start, region.end)) {
    return;
  }

  auto header_info = GetVcfHeaderInfo(hdr);
  if (!header_info.has_tumor_normal) {
    // Somatic TN Workflow but VCF doesn't have both tumor and normal
    return;
  }
  // Compute BAM and VCF Features here
  auto [vcf_features, bam_features] =
      ComputeBamAndVcfFeaturesForRegion(_global_ctx, *_worker_ctx, region, bed_regions, interest_regions);
  TumorNormalProcessing tumor_normal_processing(_global_ctx, region.start, vcf_features, bam_features);

  while (const auto& record = reader.GetNextRegionRecord(BCF_UN_ALL)) {
    for (auto& record_copy :
         SplitMultiAllelicRecord(record, hdr, _global_ctx.vcf_info_metadata, _global_ctx.vcf_fmt_metadata, true)) {
      const auto& chrom = record_copy->Chromosome();
      if (record_copy->Position() < 0) {
        WarnAsErrorIfSet("Found VCF record with position less than 0");
        continue;
      }
      const auto pos = record_copy->Position();

      if (chrom != region.chrom || std::cmp_less(pos, region.start) || std::cmp_greater_equal(pos, region.end)) {
        continue;
      }

      const DepthTuple& normalize_target = normalize_targets.GetValue(region.chrom);

      // Filter Somatic for region
      // TODO: Filter per sample
      tumor_normal_processing.ProcessRecord(
          record_copy, out_records, somatic_calculator, normalize_target, header_info);
    }
  }
}

void FilterRegionClass::FilterTumorOnlyTeRegion(const TargetRegion& region, FlowContext& flow_ctx) {
  // This function is intended to be used within a single thread (e.g. one TaskFlow task), but multiple regions can be
  // processed in parallel in the tumor-only te workflow.
  // For the specified target region, BAM and VCF features are computed, then variants are filtered by
  // chromosomal position, with a failure reason assigned to a variant if it is not considered to be somatic variation
  // based on the ML filtering or other filtering criteria, and written to the output buffer.

  // Load values for worker and global contexts
  auto& reader = _worker_ctx->filter_variants_reader;

  const auto& somatic_calculator = _worker_ctx->calculators.front();

  auto& out_records = flow_ctx.out_records;

  const auto& scoring_cols = _global_ctx.model_config.scoring_cols;
  const auto& bed_regions = _global_ctx.bed_regions;
  const auto& interest_regions = _global_ctx.interest_regions;
  const auto& normalize_targets = _global_ctx.normalize_targets;
  const bool make_forcecalls = _global_ctx.force_calls.has_value();
  const auto& block_list = _global_ctx.block_list.value_or(StrUnorderedSet{});
  const auto& hotspots = _global_ctx.hotspots.value_or(StrUnorderedSet{});

  if (!reader.SetRegion(region.chrom, region.start, region.end)) {
    return;
  }

  // Compute BAM and VCF Features here
  auto [vcf_features, bam_features] =
      ComputeBamAndVcfFeaturesForRegion(_global_ctx, *_worker_ctx, region, bed_regions, interest_regions);

  // Compute thresholds and setup filtering settings for somatic filtering
  auto weight_counts_thresholds =
      CalculateWeightedCountThresholdsPerSubstitutionType(bam_features, 0, _global_ctx.weighted_counts_threshold);
  const auto filter_settings = FilterSettings(_global_ctx.bam_feat_params.min_mapq,
                                              _global_ctx.bam_feat_params.min_bq,
                                              weight_counts_thresholds,
                                              _global_ctx.ml_threshold,
                                              block_list,
                                              _global_ctx.hotspot_weighted_counts_threshold,
                                              _global_ctx.hotspot_ml_threshold,
                                              _global_ctx.min_allele_freq_threshold,
                                              hotspots);
  const auto phased_filter_settings = PhasedFilterSettings(_global_ctx.bam_feat_params.min_mapq,
                                                           _global_ctx.bam_feat_params.min_bq,
                                                           block_list,
                                                           _global_ctx.min_phased_allele_freq,
                                                           _global_ctx.max_phased_allele_freq,
                                                           _global_ctx.min_alt_counts);

  TumorOnlyTeProcessing tumor_only_te_processing(
      _global_ctx, filter_settings, phased_filter_settings, region.start, vcf_features, bam_features);

  while (const auto& record = reader.GetNextRegionRecord(BCF_UN_ALL)) {
    const auto chrom = record->Chromosome();
    if (record->Position() < 0) {
      WarnAsErrorIfSet("Found VCF record with position less than 0");
      continue;
    }
    const auto pos = static_cast<u64>(record->Position());
    if (chrom != region.chrom || pos < region.start || pos >= region.end) {
      continue;
    }

    if (IsMultiAllelicRecord(record)) {
      throw error::Error(
          fmt::format("VCF file contains a variant with more than 2 alleles: {}:{}. "
                      "Run 'gatk LeftAlignAndTrimVariants -R <genome> -V <this VCF> -O <new VCF> "
                      "--split-multi-allelics' to normalize the VCF file.",
                      chrom,
                      pos + 1));
    }

    const auto& [ref_trimmed, alt_trimmed] = TrimVariant(record->Allele(0), record->Allele(1));
    auto key = GetVariantCorrelationKey(chrom, pos, ref_trimmed, alt_trimmed);
    auto key_unpadded = GetVariantCorrelationKey(chrom, pos, ref_trimmed, alt_trimmed, false);
    const VariantId vid(chrom, pos, ref_trimmed, alt_trimmed);
    auto variant_feature =
        tumor_only_te_processing.ScoreVariant(vid, normalize_targets, scoring_cols, somatic_calculator);

    if (make_forcecalls) {
      tumor_only_te_processing.CreateNewForceRecordsForInbetweenPositions(pos, chrom, out_records);
    }

    tumor_only_te_processing.ProcessRecord(
        vid, record, variant_feature, out_records, make_forcecalls, key, key_unpadded);
  }

  // Final force calls

  if (make_forcecalls) {
    tumor_only_te_processing.CreateNewForceRecordsForInbetweenPositions(region.end, region.chrom, out_records);
  }
}

void FilterRegionClass::FilterGermlineRegion(const TargetRegion& region, FlowContext& flow_ctx) {
  // This function is intended to be used within a single thread (e.g. one TaskFlow task), but multiple regions can be
  // processed in parallel in the germline filtering workflow.
  // For the specified target region, BAM and VCF features are computed, then variants are filtered by
  // chromosomal position, assigned a genotype, and written to the output buffer.

  // Load values for worker and global contexts
  auto& reader = _worker_ctx->filter_variants_reader;

  auto& out_records = flow_ctx.out_records;

  const auto& bed_regions = _global_ctx.bed_regions;
  const auto& interest_regions = _global_ctx.interest_regions;
  const auto& normalize_targets = _global_ctx.normalize_targets;

  if (!reader.SetRegion(region.chrom, region.start, region.end)) {
    return;
  }

  // Flush the temporary storage of variant records and features from any prior region
  _tmp_vids = {};
  _tmp_records = {};
  _tmp_genotypes = {};
  _tmp_alt_wildcard_records = {};

  // Compute BAM and VCF Features here
  auto [vcf_features, bam_features] =
      ComputeBamAndVcfFeaturesForRegion(_global_ctx, *_worker_ctx, region, bed_regions, interest_regions);

  std::optional<u64> prev_pos = std::nullopt;
  std::optional<u64> prev_pass_ref_max_pos = std::nullopt;  // max position of previous variant's ref allele
  // @TODO: initialize this variable using the last variant from the previous region

  const DepthTuple& normalize_target = normalize_targets.GetValue(region.chrom);
  while (const auto& record = reader.GetNextRegionRecord(BCF_UN_ALL)) {
    const auto& chrom = record->Chromosome();
    if (record->Position() < 0) {
      WarnAsErrorIfSet("Found VCF record with position less than 0");
      continue;
    }
    const auto& pos = record->Position();

    // reset the previous pass position if the current record is not in the target region
    if (chrom != region.chrom || std::cmp_less(pos, region.start) || std::cmp_greater_equal(pos, region.end)) {
      prev_pass_ref_max_pos = std::nullopt;
      continue;
    }
    // First reconcile and output the records for any upstream positions, then process the current record.
    ReconcileGermlineFeatures(region, pos, prev_pos, prev_pass_ref_max_pos, out_records);
    FilterGermlineRecord(record, vcf_features, bam_features, normalize_target, out_records);
  }
  if (!_tmp_records.empty() && prev_pos.has_value()) {
    // process variants at the very last position in the target region
    if (IsHaploid(region, prev_pos.value(), prev_pos.value() + 1)) {
      ReconcileGermlineHaploidRecords(_tmp_vids, _tmp_records, _tmp_genotypes, out_records);
    } else {
      ReconcileGermlineDiploidRecords(_tmp_vids,
                                      _tmp_records,
                                      _tmp_alt_wildcard_records,
                                      _tmp_genotypes,
                                      out_records,
                                      _global_ctx.vcf_info_metadata,
                                      _global_ctx.vcf_fmt_metadata);
    }
  }
}

void FilterRegionClass::FilterGermlineTaggingRegion(const TargetRegion& region, FlowContext& flow_ctx) {
  // Load values for worker and global contexts
  auto& reader = _worker_ctx->filter_variants_reader;

  auto& out_records = flow_ctx.out_records;

  const auto& bed_regions = _global_ctx.bed_regions;
  const auto& interest_regions = _global_ctx.interest_regions;
  const auto& normalize_targets = _global_ctx.normalize_targets;

  if (!reader.SetRegion(region.chrom, region.start, region.end)) {
    return;
  }

  // Flush the temporary storage of variant records and features from any prior region
  _tmp_vids = {};
  _tmp_records = {};
  _tmp_genotypes = {};
  _tmp_alt_wildcard_records = {};

  // Compute BAM and VCF Features here
  auto [vcf_features, bam_features] =
      ComputeBamAndVcfFeaturesForRegion(_global_ctx, *_worker_ctx, region, bed_regions, interest_regions);

  std::optional<u64> prev_pos = std::nullopt;

  const DepthTuple& normalize_target = normalize_targets.GetValue(region.chrom);
  while (const auto& record = reader.GetNextRegionRecord(BCF_UN_ALL)) {
    const auto& chrom = record->Chromosome();
    if (record->Position() < 0) {
      WarnAsErrorIfSet("Found VCF record with position less than 0");
      continue;
    }
    const auto& pos = record->Position();

    if (chrom != region.chrom || std::cmp_less(pos, region.start) || std::cmp_greater_equal(pos, region.end)) {
      continue;
    }

    // Reconcile and output the records for any upstream positions
    ReconcileGermlineTaggingFeatures(region, pos, prev_pos, out_records);

    if (HasPassFilter(record)) {
      // PASS filter found; this is not a candidate for germline tagging
      auto record_copy = record->Clone(_hdr);
      // add the SOMATIC INFO field flag
      record_copy->AddInfoFieldFlag(kGermlineTaggingInfoSomaticId);
      out_records.emplace_back(record_copy);
    } else {
      // PASS filter not found; this is a candidate for germline tagging
      FilterGermlineTaggingRecord(record, vcf_features, bam_features, normalize_target);
    }
  }
  if (!_tmp_records.empty() && prev_pos.has_value()) {
    // process variants at the very last position in the target region
    const u64 prev_pos_end = prev_pos.value() + _tmp_vids.front().ref.size();
    const bool is_haploid = IsHaploid(region, prev_pos.value(), prev_pos_end);
    ReconcileGermlineTaggingRecords(
        _tmp_records, _tmp_genotypes, is_haploid, _global_ctx.vcf_normal_index.value(), out_records);
  }
}

}  // namespace xoos::svc
