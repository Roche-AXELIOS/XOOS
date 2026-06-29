#pragma once
#include "core/genotype.h"
#include "core/score-calculator.h"
#include "xoos/io/vcf/vcf-record.h"
#include "xoos/types/vec.h"

namespace xoos::svc {

/**
 * @brief Update a VCF record with the specified genotype and adjust FILTER and INFO fields accordingly.
 * @param record VCF record to be updated
 * @param genotype Predicted genotype
 * @param score GenotypeScore struct for updating genotype quality and variant quality
 * @param sample_idx Index of the sample in the VCF record to update (-1 if not applicable)
 */
void UpdateGermlineTaggingRecord(const io::VcfRecordPtr& record,
                                 Genotype genotype,
                                 const GenotypeScore& score,
                                 s32 sample_idx);

/**
 * @brief Reconcile germline tagging records based on ML genotypes and haploid/diploid status.
 * @param in_records Input VCF records
 * @param ml_genotypes Predicted genotypes for the records
 * @param is_haploid Whether the region is haploid
 * @param normal_sample_idx Index of the normal sample in the VCF record (-1 if not applicable)
 * @param out_records Output VCF records after reconciliation
 */
void ReconcileGermlineTaggingRecords(const vec<io::VcfRecordPtr>& in_records,
                                     const vec<GenotypeScore>& ml_genotypes,
                                     bool is_haploid,
                                     s32 normal_sample_idx,
                                     vec<io::VcfRecordPtr>& out_records);

}  // namespace xoos::svc
