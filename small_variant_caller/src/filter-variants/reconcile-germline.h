#pragma once

#include <optional>
#include <string>

#include <xoos/io/vcf/vcf-header.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/io/vcf/vcf-writer.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/yc-decode/yc-decoder.h>

#include "core/genotype.h"
#include "core/variant-id.h"

namespace xoos::svc {

/**
 * @brief Get the indexes of genotypes in a vector of VariantId objects.
 * @param vids Vector of VariantId objects representing variants.
 * @param gt01_12_indexes Vector of indexes for genotypes classified as kGT01 or kGT12.
 * @param is_snv_ins_vec Vector indicating whether each variant is a SNV or insertion.
 * @param gt12_indexes Vector of indexes for genotypes classified as kGT12.
 * @param num_fail Count of genotypes classified as failures (not kGT01, kGT11, or kGT12).
 * @param gt01_11_indexes Vector of indexes for genotypes classified as kGT01 or kGT11.
 * @param ml_genotypes Vector of predicted genotypes for the variants.
 */
void GetGTIndexes(const vec<VariantId>& vids,
                  vec<size_t>& gt01_12_indexes,
                  vec<bool>& is_snv_ins_vec,
                  vec<size_t>& gt12_indexes,
                  u32& num_fail,
                  vec<size_t>& gt01_11_indexes,
                  const vec<Genotype>& ml_genotypes);

/**
 * @brief Creates a multi-allelic variant record by combining two single-allelic VCF records.
 *
 * This function takes two VCF records that represent different alleles at the same genomic position
 * and combines them into a single multi-allelic record. The records must share the same chromosome
 * and position but have different alternative alleles.
 *
 * @param record1 Pointer to the first VCF record containing one allele
 * @param record2 Pointer to the second VCF record containing another allele at the same position
 *
 * @return true if the multi-allelic record was successfully created, false otherwise
 *         (e.g., if records are incompatible, at different positions, or have conflicting data)
 *
 * @note The function assumes both input records are valid and non-null
 * @note Sample-specific data (genotypes, format fields) from both records will be merged appropriately
 * @note INFO fields will be combined according to VCF specification for multi-allelic variants
 */
bool MergeMultiAllelicRecords(const io::VcfRecordPtr& record1, const io::VcfRecordPtr& record2);

/**
 * @brief Update VCF record fields based on genotype classification for germline variants.
 *
 * This function updates VCF INFO and FORMAT fields to reflect the predicted genotype classification
 * from germline variant filtering. It is specifically designed for updating records from GATK
 * Haplotype Caller VCF output files in the germline workflow.
 *
 * The function performs the following updates based on the genotype:
 * - Sets appropriate FILTER status (PASS for valid genotypes, FAIL for invalid ones)
 * - Updates INFO fields: AC (allele count), AF (allele frequency), AN (allele number)
 * - Updates FORMAT field: GT (genotype)
 *
 * For failed genotypes (kGT00, kGT0, and other invalid genotypes), the function normalizes
 * the genotype values to ensure VCF compliance:
 * - Haploid failures are forced to GT=1, AC=1, AF=0, AN=1
 * - Diploid failures are forced to GT=0/1, AC=1, AF=0, AN=2
 *
 * @param record VCF record to be updated. Must be a valid VcfRecordPtr.
 * @param genotype The predicted genotype classification from the germline filtering model.
 *                 Valid values include:
 *                 - kGT01: Heterozygous diploid (0/1)
 *                 - kGT11: Homozygous alternate diploid (1/1)
 *                 - kGT12: Compound heterozygous diploid (1/2)
 *                 - kGT1: Haploid alternate (1)
 *                 - kGT0: Haploid reference (0) - treated as failure
 *                 - kGT00 and others: Treated as failures
 *
 * @throws std::runtime_error If required INFO/FORMAT fields (AC, AF, AN, GT) are not present
 *                           in the VCF record header or cannot be set.
 *
 * @note This function assumes the VCF record contains the standard GATK fields and may not
 *       work correctly with VCF files from other variant callers.
 * @note The function modifies the input record in-place.
 * @note For multi-allelic records that have been split, ensure AC and AF arrays are properly
 *       sized before calling this function.
 */
void UpdateGenotypeRelatedFields(const io::VcfRecordPtr& record, Genotype genotype);

/**
 * @brief Fail a VCF record and update relevant fields.
 * @param record VCF record
 */
void FailRecord(const io::VcfRecordPtr& record, const std::string& fail_id, const Genotype& genotype);

/**
 * @brief Helper function to update and add failed VCF records at a position for the germline workflow.
 * @param in_records Vector of VCF records
 * @param out_records Vector to hold output VCF records
 * @param skip_indexes Vector of indexes to skip
 * @param ml_genotypes Vector of predicted genotype for the records
 */
void FailAndAddGermlineDiploidRecords(const vec<io::VcfRecordPtr>& in_records,
                                      vec<io::VcfRecordPtr>& out_records,
                                      const vec<size_t>& skip_indexes,
                                      const vec<Genotype>& ml_genotypes);

/**
 * @brief Reconcile one or more ML-predicted genotypes of variant records at a single diploid chromosomal position for
 * the germline workflow.
 *
 * This function is intended for variants in diploid regions (autosomes, male XY PARs, female X chromosome) only.
 * There should be at most 2 variants at a position because the organism is expected to be diploid.
 * Since the genotypes are expected to be independently predicted by machine learning models, this function resolves
 * ambiguities and conflicts in genotype assignments for variants at the same position.
 *
 * The function implements sophisticated logic to handle multiple variants at the same position:
 * - Single variants: Direct genotype assignment with wildcard ALT support for GT=1/2 or GT=1/1
 * - Two GT=1/2 variants: Combined into multi-allelic record with GT=1/2
 * - Single passing variants: Assigns appropriate genotype while failing others
 * - Mixed variant types: Combines SNVs/insertions with deletions when both have GT=0/1 or GT=1/2
 * - Ambiguous cases: Fails all variants with appropriate filter reasons
 *
 * @param vids Vector of variant IDs for all variants at this position. Must be non-empty.
 * @param in_records Vector of VCF records corresponding to the variant IDs. Must match `vids` size.
 * @param alt_wildcard_records Vector of VCF records with ALT="*" at this position. These correspond to an upstream
 * deletion overlapping with the GT=1/2 and GT=1/1 variants at this position. Can be empty.
 * @param out_records Output vector to append processed VCF records. Records are modified in-place
 *                   and include updated genotype fields, filter status, and format fields.
 * @param ml_genotypes Vector of genotype predicted for each input record. Must match `vids` size.
 *
 * @return Optional maximum reference position spanned by the REF allele of any variant that
 *         were assigned GT=0/1 or GT=1/2. Used to determine validity of wildcard ALT alleles at downstream positions.
 *         Returns nullopt if no variants were assigned GT=0/1 or GT=1/2.
 *
 * @details Genotype Assignment Logic:
 * - CASE 1: Single variant
 *   - GT=1/2 or GT=1/1 + wildcard: Creates multi-allelic record with GT=1/2
 *   - GT=1/2 alone: Fails with "multiallele_partner" filter
 *   - GT=0/1, GT=1/1: Direct assignment
 *   - GT=0/0: Fails with "false_positive" filter
 *
 * - CASE 2: Two passing variants, both GT=1/2
 *   - Creates multi-allelic record (shorter ALT for allele 1, longer ALT for allele 2)
 *   - All other variants, if they exist: Fail with "false_positive" filter
 *
 * - CASE 3: Single passing variant, either GT=0/1 or GT=1/1
 *   - Assigns genotype to passing variant
 *   - All other variants: Fail with "false_positive" filter
 *
 * - CASE 4: Two passing variants, a deletion and a SNV/insertion, both with either GT=0/1 or GT=1/2
 *   - Creates multi-allelic record with GT=1/2
 *   - All other variants, if they exist: Fail with "false_positive" filter
 *
 * - CASE 5: GT=0/0 or conflicting genotypes
 *   - Fails all variants with appropriate filters
 *
 * @note Multi-allelic Record Creation:
 * - Combines INFO fields (AC, AF, AN, etc.) according to VCF specification
 * - Merges FORMAT fields (AD, DP, etc.) for both alleles
 * - Preserves original GATK fields in custom FORMAT fields
 * - May fail if REF/ALT sequences cannot be properly formatted
 *
 * @note Filter Assignment:
 * - PASS: Valid genotypes (GT=0/1, GT=1/1, GT=1/2)
 * - FAIL: Invalid genotypes with specific failure reasons:
 *   - false_positive: GT=0/0 (classified as noise)
 *   - multiallele_partner: GT=1/2 without valid partner
 *   - multiallele_conflict: Valid genotype but conflicts with position consensus
 *   - multiallele_format: Cannot format multi-allelic REF/ALT sequences
 *
 * @warning Input Validation:
 * - All input vectors must have matching sizes (except alt_wildcard_records)
 * - VCF records must contain required INFO/FORMAT fields for updates
 * - Variant IDs must represent variants at the same chromosomal position
 *
 * @warning Memory and Performance:
 * - Input VCF records are modified in-place
 * - Processing time scales with number of variants at position (typically 1-2 variants expected)
 *
 * @see UpdateGenotypeRelatedFields For genotype field updates
 * @see MergeMultiAllelicRecords For multi-allelic record creation
 */
std::optional<s64> ReconcileGermlineDiploidRecords(const vec<VariantId>& vids,
                                                   const vec<io::VcfRecordPtr>& in_records,
                                                   const vec<io::VcfRecordPtr>& alt_wildcard_records,
                                                   const vec<Genotype>& ml_genotypes,
                                                   vec<io::VcfRecordPtr>& out_records);

/**
 * @brief Reconcile one or more ML-predicted genotypes of variant records at a single haploid chromosomal position for
 * the germline workflow.
 * @param vids Vector of variant IDs for all variants at this position. Must be non-empty.
 * @param in_records Vector of input VCF records corresponding to the variant IDs. Must match `vids` size.
 * @param ml_genotypes Vector of genotype predicted for each input record. Must match `vids` size.
 * @param out_records Vector to hold output VCF records
 *
 * This function is intended for variants in haploid regions (non-PAR of male XY chromosomes) only.
 *
 * @details Genotype and Filter Assignment Logic:
 * - Haploid genotypes are either GT=0 or GT=1, which require a binary classification.
 * - We can reuse the existing multi-class models originally intended for diploid chromosomes.
 * - We can simply treat the class labels differently to yield exactly 2 classes:
 *   - PASS: model predicted GT=1/1, assign GT=1
 *   - FAIL: model predicted either GT=0/0, GT=0/1, or GT=1/2, fail with "false_positive" filter
 *
 * @warning Input Validation:
 * - All input vectors must have matching sizes
 * - VCF records must contain required INFO/FORMAT fields for updates
 * - Variant IDs must represent variants at the same chromosomal position
 *
 * @warning Memory and Performance:
 * - Input VCF records are modified in-place
 * - Processing time scales with number of variants at position (typically 1-2 variants expected)
 *
 * @see UpdateGenotypeRelatedFields For genotype field updates
 */
void ReconcileGermlineHaploidRecords(const vec<VariantId>& vids,
                                     const vec<io::VcfRecordPtr>& in_records,
                                     const vec<Genotype>& ml_genotypes,
                                     vec<io::VcfRecordPtr>& out_records);

/**
 * @brief Reconcile a single variant record in a diploid region for the germline workflow.
 *
 * This function processes a single variant record at a diploid chromosomal position. Depending on the predicted
 * genotype, and if there are any wildcard ALT records, it updates the record accordingly.
 *
 * @param vids A vector of Variant IDs for the variant at this position. Should only be a single entry
 * @param in_records A vector of VCF records corresponding to the variant IDs. Should only be a single entry
 * @param ml_genotypes A vector of ML predicted genotypes for the records. Should only be a single entry
 * @param alt_wildcard_records A vector of VcfRecordPtr containing wildcard ALT records, which are used to handle
 * upstream deletions overlapping with the GT=1/2 and GT=1/1 variants at this position.
 * @param out_records A vector to store output VCF records. The function will append the updated record to this vector.
 * @param gt12_01_max_ref_pos An optional reference position is used to determine validity of wildcard ALT alleles at
 * downstream positions.
 *
 * @note This function is designed to be called directly by ReconcileGermlineDiploidRecords when there is only one
 * variant at the position.
 */
void ReconcileGermlineDiploidSingleRecord(const vec<VariantId>& vids,
                                          const vec<io::VcfRecordPtr>& in_records,
                                          const vec<Genotype>& ml_genotypes,
                                          const vec<io::VcfRecordPtr>& alt_wildcard_records,
                                          vec<io::VcfRecordPtr>& out_records,
                                          std::optional<s64>& gt12_01_max_ref_pos);

/**
 * @brief Reconcile two passing GT=1/2 records at a diploid chromosomal position for the germline workflow.
 *
 * This function processes two variant records at a diploid chromosomal position that both have GT=1/2. It combines them
 * into a multi-allelic record if possible and updates the output records accordingly. If they cannot be combined, it
 * fails the records with appropriate filters.
 *
 * @param vids A vector of Variant IDs for the variant at this position.
 * @param in_records A vector of VCF records corresponding to the variant IDs.
 * @param ml_genotypes A vector of ML predicted genotypes for the records.
 * @param out_records A vector to store output VCF records. The function will append the updated record to this vector.
 * @param gt12_01_max_ref_pos An optional reference position is used to determine validity of wildcard ALT alleles at
 * downstream positions.
 * @param gt12_indexes A vector of indexes for the GT=1/2 records, used to identify which records to process.
 *
 * @note This function is designed to be called directly by ReconcileGermlineDiploidRecords when there are two passing
 * GT=1/2 records at the position.
 */
void ReconcileGermlineDiploidTwoPassingGT12Records(const vec<VariantId>& vids,
                                                   const vec<io::VcfRecordPtr>& in_records,
                                                   const vec<Genotype>& ml_genotypes,
                                                   vec<io::VcfRecordPtr>& out_records,
                                                   std::optional<s64>& gt12_01_max_ref_pos,
                                                   vec<size_t>& gt12_indexes);

/**
 * @brief Reconcile a single passing GT=0/1 or GT=1/1 record from multiple variants at a diploid chromosomal position
 * for the germline workflow.
 *
 * This function processes multiple variant records at a diploid chromosomal position where only one variant has a
 * passing genotype of either GT=0/1 or GT=1/1. It assigns the appropriate genotype to the passing variant and fails all
 * other variants as False Positives
 *
 * @param vids A vector of Variant IDs for the variant at this position.
 * @param in_records A vector of VCF records corresponding to the variant IDs.
 * @param ml_genotypes A vector of ML predicted genotypes for the records.
 * @param out_records A vector to store output VCF records. The function will append the updated records to this vector.
 * @param gt12_01_max_ref_pos An optional reference position is used to determine validity of wildcard ALT alleles at
 * downstream positions.
 * @param gt01_11_indexes A vector of indexes for the GT=0/1 and GT=1/1 records, used to identify which records to
 * process.
 *
 * @note This function is designed to be called directly by ReconcileGermlineDiploidRecords when there is a single
 * passing GT=0/1 or GT=1/1 record at the position.
 */
void ReconcileGermlineDiploidSinglePassingFromMultiple(const vec<VariantId>& vids,
                                                       const vec<io::VcfRecordPtr>& in_records,
                                                       const vec<Genotype>& ml_genotypes,
                                                       vec<io::VcfRecordPtr>& out_records,
                                                       std::optional<s64>& gt12_01_max_ref_pos,
                                                       vec<size_t>& gt01_11_indexes);

/**
 * @brief Reconcile two passing variants, a deletion and an SNV/insertion, both with either GT=0/1 or GT=1/2 at a
 * diploid chromosomal position for the germline workflow.
 *
 * This function processes two variant records at a diploid chromosomal position where one is a deletion and the other
 * is an SNV/insertion, both having either GT=0/1 or GT=1/2. It combines them into a multi-allelic record if possible
 * and updates the output records accordingly. If there are any other variants at this position they are failed as False
 * Positives.
 *
 * @param vids A vector of Variant IDs for the variant at this position.
 * @param in_records A vector of VCF records corresponding to the variant IDs.
 * @param ml_genotypes A vector of ML predicted genotypes for the records.
 * @param out_records A vector to store output VCF records. The function will append the updated records to this vector.
 * @param gt12_01_max_ref_pos An optional reference position is used to determine validity of wildcard ALT alleles at
 * downstream positions.
 * @param gt01_12_indexes A vector of indexes for the GT=0/1 and GT=1/2 records, used to identify which records to
 * process.
 *
 * @note This function is designed to be called directly by ReconcileGermlineDiploidRecords when there are two passing
 * variants, a deletion and an SNV/insertion, both with either GT=0/1 or GT=1/2 at the position.
 */
void ReconcileGermlineDiploidTwoPassingMix(const vec<VariantId>& vids,
                                           const vec<io::VcfRecordPtr>& in_records,
                                           const vec<Genotype>& ml_genotypes,
                                           vec<io::VcfRecordPtr>& out_records,
                                           std::optional<s64>& gt12_01_max_ref_pos,
                                           vec<size_t>& gt01_12_indexes);

}  // namespace xoos::svc
