#pragma once

#include <optional>

#include <xoos/io/vcf/vcf-header.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "core/variant-info.h"
#include "util/parallel-compute-utils.h"
#include "util/region-util.h"
#include "util/vcf-util.h"

namespace xoos::svc {

/**
 * @brief Update a VCF record for the given features within the germline workflow.
 * @param record VCF record to be updated
 * @param bam_feat BAM features
 * @param vcf_feat VCF features
 * @param is_duplex_protocol Whether the sequencing protocol is duplex
 */
void UpdateGermlineRecordWithFeatures(const io::VcfRecordPtr& record,
                                      const BamFeatureTuple& bam_feat,
                                      const VcfFeature& vcf_feat,
                                      bool is_duplex_protocol);

class FilterRegionClass {
 public:
  FilterRegionClass(const GlobalContext& global_ctx, std::unique_ptr<WorkerContext>& worker_ctx, io::VcfHeaderPtr hdr)
      : _global_ctx(&global_ctx), _worker_ctx(std::move(worker_ctx)), _hdr(std::move(hdr)) {
  }

  /**
   * @brief Update the worker context with a new WorkerContext instance.
   * @param worker_ctx A unique pointer to the new WorkerContext instance.
   */
  void UpdateWorkerCtx(std::unique_ptr<WorkerContext>& worker_ctx);

  /**
   * @brief Filter variants in a specified region for the `germline-wgs` workflow.
   *
   * This function is intended to be used within a single thread, but can be run in parallel for different regions. For
   * the given region BAM and VCF features are computed, then variants are filtered by chromosomal position, assigned a
   * genotype, and written to the output buffer.
   *
   * @param result RegionResult for output records and SHAP value rows.
   * @param region Target region to filter variants.
   *
   * @note This function assumes that the `_worker_ctx` and `_global_ctx` are properly initialized and contain the
   * necessary information for filtering.
   */
  void FilterGermlineRegion(const TargetRegion& region, RegionResult& result);

  /**
   * @brief Filter variants in a specified region for the `tumor-only-te` workflow.
   *
   * This function is intended to be used within a single thread, but can be run in parallel for different regions. For
   * the given region BAM and VCF features are computed, then variants are filtered by chromosomal position, assigned a
   * failure reason if not considered passing somatic variation based on the ML filtering or other filtering criteria,
   * and written to the output buffer.
   *
   * @param result RegionResult for output records and SHAP value rows.
   * @param region Target region to filter variants.
   *
   * @note This function assumes that the `_worker_ctx` and `_global_ctx` are properly initialized and contain the
   * necessary information for filtering.
   */
  void FilterTumorOnlyTeRegion(const TargetRegion& region, RegionResult& result);

  /**
   * @brief Filter variants in a specified region for the `tumor-normal-wgs` workflow.
   *
   * This function is intended to be used within a single thread, but can be run in parallel for different regions. For
   * the given region BAM and VCF features are computed, then variants are filtered by chromosomal position, assigned a
   * failure reason if not considered passing somatic variation based on the ML filtering or other filtering criteria,
   * and written to the output buffer.
   *
   * @param result RegionResult for output records and SHAP value rows.
   * @param region Target region to filter variants.
   *
   * @note This function assumes that the `_worker_ctx` and `_global_ctx` are properly initialized and contain the
   * necessary information for filtering.
   */
  void FilterTumorNormalRegion(const TargetRegion& region, RegionResult& result);

  /**
   * @brief Filter variants in the given region for the `tumor-only-wgs` workflow.
   *
   * Single-threaded per-region worker. Computes single-sample BAM and VCF features for the region, then scores each
   * variant with the somatic ML model and applies hard filters. There is no normal sample and no germline tagging.
   *
   * @param region Target region to filter variants.
   * @param result RegionResult for output records and SHAP value rows.
   *
   * @note Assumes `_worker_ctx` and `_global_ctx` are properly initialized.
   */
  void FilterTumorOnlyWgsRegion(const TargetRegion& region, RegionResult& result);

  void FilterGermlineTaggingRegion(const TargetRegion& region, RegionResult& result);

  /**
   * @brief Filter variants in a specified region for the `germline-mitochondria` workflow.
   *
   * This function is intended to be used within a single thread, but can be run in parallel for different regions. For
   * the given region BAM features are computed (honoring `duplex_lowbq`), then each variant is assigned a genotype and
   * `low_af` / `low_dp` FILTER tags based on its allele frequency and depth against the mitochondria thresholds, and
   * written to the output buffer. This workflow does not use an ML model.
   *
   * @param result RegionResult for output records.
   * @param region Target region to filter variants.
   *
   * @note This function assumes that the `_worker_ctx` and `_global_ctx` are properly initialized and contain the
   * necessary information for filtering.
   */
  void FilterMitochondriaRegion(const TargetRegion& region, RegionResult& result);

  /**
   * @brief Reconcile germline features for the given region and position, updating the previous position and maximum
   * reference position as needed.
   *
   * This function reconciles predicted genotypes and updates VCF records based on the features at the specified
   * position within the target region. It also manages the state of previous positions to ensure correct processing of
   * variants.
   *
   * @param region Target region being processed.
   * @param pos Current 0-based genomic position in the VCF record.
   * @param prev_pos Reference to the previous 0-based position, updated if necessary.
   * @param prev_pass_ref_max_pos Reference to the maximum 0-based position of the previous variant's reference allele,
   * updated if necessary.
   * @param prev_span_deletion_fields Per-allele FORMAT field snapshot of the most recent spanning deletion, carried
   * across positions alongside prev_pass_ref_max_pos and applied to wildcard (`*`) records so they inherit the
   * deletion's computed counts.
   * @param out_records Output records to which reconciled records will be added.
   *
   * @note This function assumes _temp_vids, _temp_records, _temp_genotypes, and _temp_alt_wildcard_records have been
   * populated for the previous position being processed.
   */
  void ReconcileGermlineFeatures(const TargetRegion& region,
                                 u64 pos,
                                 std::optional<u64>& prev_pos,
                                 std::optional<u64>& prev_pass_ref_max_pos,
                                 PerAlleleFieldSnapshot& prev_span_deletion_fields,
                                 vec<io::VcfRecordPtr>& out_records);

  /**
   * @brief Reconcile germline-tagging features for the given region and position, updating the previous position as
   * needed.
   *
   * This function reconciles predicted genotypes and updates VCF records based on the features at the specified
   * position within the target region. It also manages the state of previous positions to ensure correct processing of
   * variants.
   *
   * @param region Target region being processed.
   * @param pos Current 0-based genomic position in the VCF record.
   * @param prev_pos Reference to the previous 0-based position, updated if necessary.
   * @param out_records Output records to which reconciled records will be added.
   *
   * @note This function assumes _tmp_vids, _temp_records, _temp_genotypes have been populated for the previous position
   * being processed.
   */
  void ReconcileGermlineTaggingFeatures(const TargetRegion& region,
                                        u64 pos,
                                        std::optional<u64>& prev_pos,
                                        vec<io::VcfRecordPtr>& out_records);

  /**
   * @brief Filter multiple alleles in a VCF record for the `germline-wgs` workflow. This function processes a VCF
   * record, pulling the relevant features, scoring them to assign a genotype and updating the records as needed.
   * @param record A VCF record to filter.
   * @param vcf_features A map of VariantId to VCF features.
   * @param bam_features A collection of BAM features.
   * @param normalize_target A normalization target for feature values.
   * @param result RegionResult for output VCF records and SHAP value rows
   *
   * @note This function directly sets _temp_vids, _temp_records, _temp_genotypes, and _temp_alt_wildcard_records as
   * needed for later reconciliation.
   */
  void FilterGermlineRecord(const io::VcfRecordPtr& record,
                            const VarIdToVcfFeatures& vcf_features,
                            const BamRegionFeatureCollection& bam_features,
                            const DepthTuple& normalize_target,
                            RegionResult& result);

  /**
   * @brief Filter multiple alleles in a VCF record for the `germline-tagging` workflow. This function processes a VCF
   * record, pulling the relevant features, scoring them to assign a genotype.
   * @param record A VCF record to filter.
   * @param vcf_features A map of VariantId to VCF features.
   * @param bam_features A collection of BAM features.
   * @param normalize_target An optional normalization target for feature values.
   * @param shap_value_rows Vector for output SHAP value rows
   *
   * @note This function directly sets _temp_vids, _temp_records, _temp_genotypes, and _temp_alt_wildcard_records as
   * needed for later reconciliation.
   */
  void FilterGermlineTaggingRecord(const io::VcfRecordPtr& record,
                                   const VarIdToVcfFeatures& vcf_features,
                                   const BamRegionFeatureCollection& bam_features,
                                   const DepthTuple& normalize_target,
                                   vec<vec<std::string>>& shap_value_rows);

  /**
   * Checks if the region is haploid based on the previous position and end position. The region is compared to the PAR
   * regions of ChrX and ChrY if it is from the X or Y chromosome, respectively and checks done to see if it is a
   * diploid position. If the region is not on these chromosomes, it is considered diploid.
   * @param region A TargetRegion object representing the genomic region.
   * @param prev_pos A previous position to compare
   * @param prev_pos_end A previous position end to compare
   * @return True is region is haploid, false otherwise.
   */
  bool IsHaploid(const TargetRegion& region, u64 prev_pos, u64 prev_pos_end) const;

  /**
   * @brief Reset the temporary storage vectors used for holding variant records, features, and genotypes during
   * processing. This function is used during region processing to clear out data from previously processed positions
   * and ensure that the temporary storage is ready for the next set of variants to be processed.
   */
  void ResetTmpStorage();

 private:
  /**
   * @brief Score a single germline allele that has both VCF and BAM features, updating the record and temporary
   * storage.
   * @param record_copy Single-allele VCF record to update and retain.
   * @param vid VariantId of the allele.
   * @param bam_feat BAM features for the allele.
   * @param vcf_feat VCF features for the allele.
   * @param normalize_target Normalization target for feature values.
   * @param result RegionResult for SHAP value rows.
   */
  void ScoreGermlineAllele(const io::VcfRecordPtr& record_copy,
                           const VariantId& vid,
                           const BamFeatureTuple& bam_feat,
                           const VcfFeature& vcf_feat,
                           const DepthTuple& normalize_target,
                           RegionResult& result);

  const GlobalContext* _global_ctx{nullptr};
  std::unique_ptr<WorkerContext> _worker_ctx;
  io::VcfHeaderPtr _hdr;

  vec<VariantId> _tmp_vids;
  vec<io::VcfRecordPtr> _tmp_records;
  vec<PredictionScore> _tmp_genotypes;
  vec<io::VcfRecordPtr> _tmp_alt_wildcard_records;
};

}  // namespace xoos::svc
