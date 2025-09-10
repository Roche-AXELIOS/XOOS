#include "filter-region.h"

#include <algorithm>
#include <memory>
#include <unordered_map>

#include <taskflow/taskflow.hpp>

#include <xoos/error/error.h>
#include <xoos/io/bed-region.h>
#include <xoos/io/fasta-reader.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/io/vcf/vcf-writer.h>
#include <xoos/log/logging.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>
#include <xoos/util/container-functions.h>
#include <xoos/util/string-functions.h>

#include "core/filtering.h"
#include "core/variant-feature-extraction.h"
#include "core/vcf-fields.h"
#include "reconcile-germline.h"
#include "util/log-util.h"
#include "util/parallel-compute-utils.h"
#include "util/seq-util.h"

namespace xoos::svc {

/**
 * @brief Update a VCF record for the given features within the germline workflow.
 * @param record VCF record to be updated
 * @param bam_feat BAM features
 * @param vcf_feat VCF features
 * @param ref_feat Reference features
 */
static void UpdateGermlineRecordWithFeatures(const io::VcfRecordPtr& record,
                                             const UnifiedVariantFeature& bam_feat,
                                             const VcfFeature& vcf_feat,
                                             const UnifiedReferenceFeature& ref_feat) {
  record->SetFormatField<s32>(kGermlineMLId, {1});
  // "lowbq" features are incremented by 0.5 per read, so they are not always an integer
  // To be consistent for both AD and DP fields, use `std::round()` to convert to nearest integer
  record->SetFormatField<s32>(kAlleleDepthId,
                              {static_cast<s32>(std::round(ref_feat.support + ref_feat.duplex_lowbq)),
                               static_cast<s32>(std::round(bam_feat.support + bam_feat.duplex_lowbq))});
  record->SetFormatField<s32>(kFieldDp, {static_cast<s32>(std::round(ref_feat.duplex_dp))});
  record->SetFormatField<f32>(kGermlineGnomadAFId, {vcf_feat.popaf});
  record->SetFormatField<f32>(kGermlineRefAvgMapqId, {static_cast<f32>(ref_feat.mapq_mean)});
  record->SetFormatField<f32>(kGermlineAltAvgMapqId, {static_cast<f32>(bam_feat.mapq_mean)});
  record->SetFormatField<f32>(kGermlineRefAvgDistId, {static_cast<f32>(ref_feat.distance_mean)});
  record->SetFormatField<f32>(kGermlineAltAvgDistId, {static_cast<f32>(bam_feat.distance_mean)});
  record->SetFormatField<f32>(kGermlineDensity100BPId, {static_cast<f32>(vcf_feat.variant_density)});
}

/**
 * @brief Copy a VCF record for a given allele within the germline workflow.
 * @param original_record VCF record to be copied
 * @param new_header VCF header for the new VCF record
 * @param allele_idx Allele index in the VCF record to be copied
 * @return The copied VCF record
 */
static io::VcfRecordPtr CopyGermlineRecord(const io::VcfRecordPtr& original_record,
                                           const io::VcfHeaderPtr& new_header,
                                           const s32 allele_idx) {
  const auto& record_copy = original_record->Clone(new_header);
  record_copy->SetAlleles({original_record->Allele(0), original_record->Allele(allele_idx)});

  // TODO: check header for IDs with `Number=A` or `Number=R` and copy the record values systematically

  // INFO fields with ALT alleles
  const auto alt_val_idx = static_cast<size_t>(allele_idx) - 1;
  const auto allele_idx_u = static_cast<size_t>(allele_idx);
  static const vec<std::string> kIntInfoFields{kFieldAc, kFieldHapcomp, kFieldMleac};
  for (const auto& field : kIntInfoFields) {
    const auto& values = original_record->GetInfoFieldNoCheck<s32>(field);
    if (!values.empty()) {
      record_copy->SetInfoField<s32>(field, {values[alt_val_idx]});
    }
  }
  static const vec<std::string> kFloatInfoFields{kFieldAf, kFieldHapdom, kFieldMleaf};
  for (const auto& field : kFloatInfoFields) {
    const auto& values = original_record->GetInfoFieldNoCheck<f32>(field);
    if (!values.empty()) {
      record_copy->SetInfoField<f32>(field, {values[alt_val_idx]});
    }
  }

  // INFO fields with REF and ALT alleles
  const auto& rpa_values = original_record->GetInfoFieldNoCheck<s32>(kFieldRpa);
  if (!rpa_values.empty()) {
    record_copy->SetInfoField<s32>(kFieldRpa, {rpa_values[0], rpa_values[allele_idx_u]});
  }
  // FORMAT fields
  const auto& dp_values = original_record->GetFormatFieldNoCheck<s32>(kFieldDp);
  // FORMAT fields with REF and ALT alleles
  const auto& ad_values = original_record->GetFormatFieldNoCheck<s32>(kFieldAd);
  if (!ad_values.empty()) {
    record_copy->SetFormatField<s32>(kFieldAd, {ad_values[0], ad_values[allele_idx_u]});
  }
  // add GATK's original DP, AD, GT, ALT
  // note that AD and ALT values from multi-allelic records are not split here
  // extract original field values and copy them to output records
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

  return record_copy;
}

/**
 * @brief Find a reference feature at a specific chromosome and position.
 * @param ref_features Map of reference features indexed by chromosome and position
 * @param chrom Chromosome name
 * @param pos Position in the chromosome
 * @return UnifiedReferenceFeature at the specified position or a zeroed feature if not found
 */
static const UnifiedReferenceFeature& FindRefFeature(const RefInfoMap& ref_features,
                                                     const std::string& chrom,
                                                     u64 pos) {
  // Find the reference feature at the specified chromosomal position.
  // Either return the non-empty UnifiedReferenceFeature struct for that position OR an empty struct with zero values.
  auto chrom_it = ref_features.find(chrom);
  if (chrom_it == ref_features.end()) {
    return kZeroUnifiedReferenceFeature;
  }
  auto pos_it = chrom_it->second.find(pos);
  if (pos_it == chrom_it->second.end()) {
    return kZeroUnifiedReferenceFeature;
  }
  return pos_it->second;
}

/**
 * @brief Helper function to get a VCF feature for a specific chromosome, position, and variant ID.
 * @param vcf_features Map of VCF features indexed by chromosome, position, and variant ID
 * @param chrom Chromosome name
 * @param pos Position in the chromosome
 * @param vid Variant ID
 * @return An optional VcfFeature if found, otherwise std::nullopt
 */
static std::optional<VcfFeature> GetVcfFeature(const ChromToVcfFeaturesMap& vcf_features,
                                               const std::string& chrom,
                                               const u64 pos,
                                               const VariantId& vid) {
  auto chrom_itr = vcf_features.find(chrom);
  if (chrom_itr != vcf_features.end()) {
    auto pos_itr = chrom_itr->second.find(pos);
    if (pos_itr != chrom_itr->second.end()) {
      auto vid_itr = pos_itr->second.find(vid);
      if (vid_itr != pos_itr->second.end()) {
        return vid_itr->second;
      }
    }
  }
  return std::nullopt;
}

/**
 * @brief Helper function to get a BAM feature for a specific chromosome, position, and variant ID.
 * @param bam_features Map of BAM features indexed by chromosome, position, and variant ID
 * @param chrom Chromosome name
 * @param pos Position in the chromosome
 * @param vid Variant ID
 * @return An optional UnifiedVariantFeature if found, otherwise std::nullopt
 */
static std::optional<UnifiedVariantFeature> GetBamFeature(const ChromToVariantInfoMap& bam_features,
                                                          const std::string& chrom,
                                                          const u64 pos,
                                                          const VariantId& vid) {
  auto chrom_itr = bam_features.find(chrom);
  if (chrom_itr != bam_features.end()) {
    auto pos_itr = chrom_itr->second.find(pos);
    if (pos_itr != chrom_itr->second.end()) {
      auto vid_itr = pos_itr->second.find(vid);
      if (vid_itr != pos_itr->second.end()) {
        return vid_itr->second;
      }
    }
  }
  return std::nullopt;
}

void FilterRegionClass::FilterMultipleAlleles(const io::VcfRecordPtr& record,
                                              const RefInfoMap& ref_features,
                                              const ChromToVcfFeaturesMap& vcf_features,
                                              const ChromToVariantInfoMap& bam_features,
                                              const std::optional<u32> normalize_target,
                                              vec<io::VcfRecordPtr>& out_records) {
  using enum VariantType;
  using enum Genotype;
  const auto& chrom = record->Chromosome();
  const auto pos = static_cast<u64>(record->Position());
  const auto snv_calculator = [&snv_calc = _worker_ctx->calculators[0]](const auto& feature_vec) {
    return snv_calc.CalculateScoreGermline(feature_vec);
  };
  const auto& indel_calculator = [&indel_calc = _worker_ctx->calculators[1]](const auto& feature_vec) {
    return indel_calc.CalculateScoreGermline(feature_vec);
  };
  const auto& snv_scoring_cols = _global_ctx.model_config.snv_scoring_cols;
  const auto& indel_scoring_cols = _global_ctx.model_config.indel_scoring_cols;
  const auto& ref = record->Allele(0);
  const bool ref_acgt_only = ContainsOnlyACTG(ref);
  const UnifiedReferenceFeature& snv_ins_ref_feat = FindRefFeature(ref_features, chrom, pos);
  // deletion reference feature position needs to be offset by 1
  const UnifiedReferenceFeature& del_ref_feat = FindRefFeature(ref_features, chrom, pos + 1);
  // a VCF record can have >=1 ALT allele
  for (s32 allele_idx = 1; allele_idx < record->NumAlleles(); ++allele_idx) {
    const auto& alt = record->Allele(allele_idx);
    const auto& record_copy = CopyGermlineRecord(record, _hdr, allele_idx);
    if (!ref_acgt_only) {
      // REF/ALT contains non-ACGT; copy the record and set FILTER to FAIL
      FailRecord(record_copy, kFilteringNonAcgtRefAltId, kGT00);
      out_records.emplace_back(record_copy);
      continue;
    }
    if (alt == "*") {
      _tmp_alt_wildcard_records.emplace_back(record_copy);
      continue;
    }
    const auto& [ref_trimmed, alt_trimmed] = TrimVariant(ref, alt);
    const VariantId vid(chrom, pos, ref_trimmed, alt_trimmed);
    std::optional<VcfFeature> found_vcf_feat = GetVcfFeature(vcf_features, chrom, pos, vid);
    std::optional<UnifiedVariantFeature> found_var_feat;
    bool feat_found = false;
    if (found_vcf_feat.has_value()) {
      feat_found = found_vcf_feat->genotype != Genotype::kGTNA;
    }
    if (feat_found) {
      found_var_feat = GetBamFeature(bam_features, chrom, pos, vid);
    }

    if (feat_found && found_var_feat.has_value()) {
      auto& var_feat = found_var_feat.value();
      auto& vcf_feat = found_vcf_feat.value();
      record_copy->SetAlleles({ref_trimmed, alt_trimmed});
      _tmp_vids.emplace_back(vid);
      _tmp_records.emplace_back(record_copy);

      // update the VCF record with ML features
      const auto& ref_feat = vid.type == kDeletion ? del_ref_feat : snv_ins_ref_feat;
      UpdateGermlineRecordWithFeatures(record_copy, var_feat, vcf_feat, ref_feat);

      // classify the preliminary genotype and append it to the temporary vector
      if (vid.type == kSNV) {
        const auto& feature_vec = GetFeatureVec(snv_scoring_cols, vid, var_feat, ref_feat, vcf_feat, normalize_target);
        _tmp_genotypes.emplace_back(snv_calculator(feature_vec));
      } else {
        const auto& feature_vec =
            GetFeatureVec(indel_scoring_cols, vid, var_feat, ref_feat, vcf_feat, normalize_target);
        _tmp_genotypes.emplace_back(indel_calculator(feature_vec));
      }
    } else {
      // features not found; copy the record and set FILTER to FAIL
      record_copy->SetAlleles({ref_trimmed, alt_trimmed});
      FailRecord(record_copy, kFilteringMissingFeatureId, kGT00);
      out_records.emplace_back(record_copy);
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
    if (!prev_pass_ref_max_pos.has_value() || prev_pass_ref_max_pos.value() < _tmp_vids[0].pos) {
      // Wildcard ALT alleles are valid only if an upstream variant overlaps this position
      _tmp_alt_wildcard_records = {};
    }

    std::optional<u64> pass_ref_max_pos;
    const u64 prev_pos_end = prev_pos.value() + _tmp_vids[0].ref.size();
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
      pass_ref_max_pos = ReconcileGermlineDiploidRecords(
          _tmp_vids, _tmp_records, _tmp_alt_wildcard_records, _tmp_genotypes, out_records);
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

bool FilterRegionClass::IsHaploid(const TargetRegion& region, const u64 prev_pos, const u64 prev_pos_end) const {
  const bool is_male = _global_ctx.sex == Sex::kMale;
  const bool is_chr_x = region.chrom == _global_ctx.chr_x_name;
  const bool is_chr_y = region.chrom == _global_ctx.chr_y_name;
  const bool at_chr_x_par = is_chr_x && IntervalOverlap(_global_ctx.chr_x_par, region.start, region.end);
  const bool at_chr_y_par = is_chr_y && IntervalOverlap(_global_ctx.chr_y_par, region.start, region.end);
  const bool x_overlap =
      (is_chr_x && (!at_chr_x_par || !IntervalOverlap(_global_ctx.chr_x_par, prev_pos, prev_pos_end)));
  const bool y_overlap =
      (is_chr_y && (!at_chr_y_par || !IntervalOverlap(_global_ctx.chr_y_par, prev_pos, prev_pos_end)));
  return is_male && (x_overlap || y_overlap);
}

void FilterRegionClass::UpdateWorkerCtx(std::unique_ptr<WorkerContext>& worker_ctx) {
  _worker_ctx = std::move(worker_ctx);
}

void FilterRegionClass::FilterGermlineRegion(const TargetRegion& region, FlowContext& flow_ctx) {
  // This function is intended to be used within a single thread (e.g. one TaskFlow task), but multiple regions can be
  // processed in parallel in the filtering workflow.
  // For the specified target region, BAM and VCF features are computed, then variants are filtered by
  // chromosomal position, assigned a genotype, and written to the output buffer.

  using enum VariantType;

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
  auto [vcf_features, bam_features, ref_features] =
      ComputeBamAndVcfFeaturesForRegion(_global_ctx, *_worker_ctx, region, bed_regions, interest_regions);

  std::optional<u64> prev_pos = std::nullopt;
  std::optional<u64> prev_pass_ref_max_pos = std::nullopt;  // max position of previous variant's ref allele
  // @TODO: initialize this variable using the last variant from the previous region

  const std::optional<u32> normalize_target = util::container::Find(normalize_targets, region.chrom);
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
    FilterMultipleAlleles(record, ref_features, vcf_features, bam_features, normalize_target, out_records);
  }
  if (!_tmp_records.empty() && prev_pos.has_value()) {
    // process variants at the very last position in the target region
    if (IsHaploid(region, prev_pos.value(), prev_pos.value() + 1)) {
      ReconcileGermlineHaploidRecords(_tmp_vids, _tmp_records, _tmp_genotypes, out_records);
    } else {
      ReconcileGermlineDiploidRecords(_tmp_vids, _tmp_records, _tmp_alt_wildcard_records, _tmp_genotypes, out_records);
    }
  }
}

}  // namespace xoos::svc
