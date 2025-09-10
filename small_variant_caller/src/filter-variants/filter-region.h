#pragma once

#include <optional>

#include <xoos/io/vcf/vcf-header.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/io/vcf/vcf-writer.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/yc-decode/yc-decoder.h>

#include "core/genotype.h"
#include "core/variant-info.h"
#include "util/parallel-compute-utils.h"
#include "util/region-util.h"

namespace xoos::svc {

class FilterRegionClass {
 private:
  const GlobalContext& _global_ctx;
  std::unique_ptr<WorkerContext> _worker_ctx;
  io::VcfHeaderPtr _hdr;

  vec<VariantId> _tmp_vids;
  vec<io::VcfRecordPtr> _tmp_records;
  vec<Genotype> _tmp_genotypes;
  vec<io::VcfRecordPtr> _tmp_alt_wildcard_records;

 public:
  FilterRegionClass(const GlobalContext& global_ctx, std::unique_ptr<WorkerContext>& worker_ctx, io::VcfHeaderPtr hdr)
      : _global_ctx(global_ctx), _worker_ctx(std::move(worker_ctx)), _hdr(std::move(hdr)) {
  }

  /**
   * @brief Update the worker context with a new WorkerContext instance.
   * @param worker_ctx A unique pointer to the new WorkerContext instance.
   */
  void UpdateWorkerCtx(std::unique_ptr<WorkerContext>& worker_ctx);

  /**
   * @brief Filter variants in a specified region for the `germline` workflow.
   *
   * This function is intended to be used within a single thread, but can be run in parallel for different regions. For
   * the given region BAM and VCF features are computed, then variants are filtered by chromosomal position, assigned a
   * genotype, and written to the output buffer.
   *
   * @param flow_ctx Flow context for managing output records.
   * @param region Target region to filter variants.
   *
   * @note This function assumes that the `_worker_ctx` and `_global_ctx` are properly initialized and contain the
   * necessary information for filtering.
   */
  void FilterGermlineRegion(const TargetRegion& region, FlowContext& flow_ctx);

  /**
   * @brief Reconcile germline features for the given region and position, updating the previous position and maximum
   * reference position as needed.
   *
   * This function reconciles predicted genotypes and updates VCF records based on the features at the specified
   * position within the target region. It also manages the state of previous positions to ensure correct processing of
   * variants.
   *
   * @param region Target region being processed.
   * @param pos Current position in the VCF record.
   * @param prev_pos Reference to the previous position, updated if necessary.
   * @param prev_pass_ref_max_pos Reference to the maximum position of the previous variant's reference allele, updated
   * if necessary.
   * @param out_records Output records to which reconciled records will be added.
   *
   * @note This function assumes _temp_vids, _temp_records, _temp_genotypes, and _temp_alt_wildcard_records have been
   * populated for the previous position being processed.
   */
  void ReconcileGermlineFeatures(const TargetRegion& region,
                                 u64 pos,
                                 std::optional<u64>& prev_pos,
                                 std::optional<u64>& prev_pass_ref_max_pos,
                                 vec<io::VcfRecordPtr>& out_records);

  /**
   * @brief Filter multiple alleles in a VCF record for the `germline` workflow. This function processes a VCF record,
   * pulling the relevant features, scoring them to assign a genotype and updating the records as needed.
   * @param record A VCF record to filter.
   * @param ref_features A map of reference features indexed by chromosome.
   * @param vcf_features A map of VCF features indexed by chromosome.
   * @param bam_features A map of BAM features indexed by chromosome.
   * @param normalize_target An optional normalization target for feature values.
   * @param out_records An output vector to store filtered VCF records.
   *
   * @note This function directly sets _temp_vids, _temp_records, _temp_genotypes, and _temp_alt_wildcard_records as
   * needed for later reconciliation.
   */
  void FilterMultipleAlleles(const io::VcfRecordPtr& record,
                             const RefInfoMap& ref_features,
                             const ChromToVcfFeaturesMap& vcf_features,
                             const ChromToVariantInfoMap& bam_features,
                             std::optional<u32> normalize_target,
                             vec<io::VcfRecordPtr>& out_records);

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
};

}  // namespace xoos::svc
