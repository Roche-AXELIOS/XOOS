#include "reconcile-germline.h"

#include <algorithm>
#include <utility>

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
#include "core/vcf-fields.h"
#include "util/seq-util.h"

namespace xoos::svc {

void GetGTIndexes(const vec<VariantId>& vids,
                  vec<size_t>& gt01_12_indexes,
                  vec<bool>& is_snv_ins_vec,
                  vec<size_t>& gt12_indexes,
                  u32& num_fail,
                  vec<size_t>& gt01_11_indexes,
                  const vec<Genotype>& ml_genotypes) {
  // Get the indexes of the variants based on their genotypes and store them in the relevant vectors. Needed for
  // downstream filtering
  using enum VariantType;
  using enum Genotype;
  auto num_variants = vids.size();
  is_snv_ins_vec.reserve(num_variants);
  for (size_t i = 0; i < num_variants; ++i) {
    const auto& vid = vids.at(i);
    is_snv_ins_vec.emplace_back(vid.type == kSNV || vid.type == kInsertion);
    switch (ml_genotypes.at(i)) {
      case kGT01:
        gt01_11_indexes.emplace_back(i);
        gt01_12_indexes.emplace_back(i);
        break;
      case kGT11:
        gt01_11_indexes.emplace_back(i);
        break;
      case kGT12:
        gt12_indexes.emplace_back(i);
        gt01_12_indexes.emplace_back(i);
        break;
      default:
        ++num_fail;
        break;
    }
  }
}

static void SetInfoFieldsForMultiAllelic(const io::VcfRecordPtr& record1, const io::VcfRecordPtr& record2) {
  // TODO: check VCF header for IDs with `Number=A` or `Number=R` and copy the record values systematically instead
  // INFO fields with ALT alleles
  static const vec<std::string> kIntInfoFields{kFieldHapcomp, kFieldMleac, kFieldRpa};
  for (const auto& field : kIntInfoFields) {
    const auto& values1 = record1->GetInfoFieldNoCheck<s32>(field);
    if (!values1.empty()) {
      const auto& values2 = record2->GetInfoFieldNoCheck<s32>(field);
      if (!values2.empty()) {
        if (field == kFieldRpa) {
          record1->SetInfoField<s32>(kFieldRpa, {values1[0], values1[1], values2[1]});
        } else {
          record1->SetInfoField<s32>(field, {values1[0], values2[0]});
        }
      }
    }
  }
  static const vec<std::string> kFloatInfoFields{kFieldHapdom, kFieldMleaf};
  for (const auto& field : kFloatInfoFields) {
    const auto& values1 = record1->GetInfoFieldNoCheck<f32>(field);
    if (!values1.empty()) {
      const auto& values2 = record2->GetInfoFieldNoCheck<f32>(field);
      if (!values2.empty()) {
        record1->SetInfoField<f32>(field, {values1[0], values2[0]});
      }
    }
  }
}

static void SetFormatFieldsForMultiAllelic(const io::VcfRecordPtr& record1, const io::VcfRecordPtr& record2) {
  // FORMAT fields with ALT alleles
  static const vec<std::string> kFloatFormatFields{
      kGermlineGnomadAFId, kGermlineRefAvgMapqId, kGermlineAltAvgMapqId, kGermlineRefAvgDistId, kGermlineAltAvgDistId};
  for (const auto& field : kFloatFormatFields) {
    const auto& values1 = record1->GetFormatFieldNoCheck<f32>(field);
    if (!values1.empty()) {
      const auto& values2 = record2->GetFormatFieldNoCheck<f32>(field);
      if (!values2.empty()) {
        record1->SetFormatField<f32>(field, {values1[0], values2[0]});
      }
    }
  }

  // FORMAT fields with REF and ALT alleles
  const auto& ad_values1 = record1->GetFormatFieldNoCheck<s32>(kFieldAd);
  if (!ad_values1.empty()) {
    const auto& ad_values2 = record2->GetFormatFieldNoCheck<s32>(kFieldAd);
    if (!ad_values2.empty()) {
      record1->SetFormatField<s32>(kFieldAd, {ad_values1[0], ad_values1[1], ad_values2[1]});
    }
  }
}

bool MergeMultiAllelicRecords(const io::VcfRecordPtr& record1, const io::VcfRecordPtr& record2) {
  // Merge two multi-allelic VCF records into one. First set the alleles and then set the info and format fields
  const auto& [ref1, alt1] = TrimVariant(record1->Allele(0), record1->Allele(1));
  const auto& [ref2, alt2] = TrimVariant(record2->Allele(0), record2->Allele(1));
  const auto& [ref_new, alt1_new, alt2_new] = FormatVariants(ref1, alt1, ref2, alt2);
  if (ref_new.empty() || alt1_new.empty() || alt2_new.empty()) {
    return false;
  }
  record1->SetAlleles({ref_new, alt1_new, alt2_new});
  UpdateGenotypeRelatedFields(record1, Genotype::kGT12);

  SetInfoFieldsForMultiAllelic(record1, record2);
  SetFormatFieldsForMultiAllelic(record1, record2);

  return true;
}

static void UpdateGenotypeRelatedFieldsHelper(const io::VcfRecordPtr& record,
                                              const Genotype genotype,
                                              const std::string& filter,
                                              const vec<s32>& ac,
                                              const vec<f32>& af,
                                              const vec<s32>& an) {
  record->AddFilter(filter);
  record->SetInfoField<s32>(kFieldAc, ac);
  record->SetInfoField<f32>(kFieldAf, af);
  record->SetInfoField<s32>(kFieldAn, an);
  record->SetFormatField<std::string>(kFieldGt, {GenotypeToString(genotype)});
}

void UpdateGenotypeRelatedFields(const io::VcfRecordPtr& record, const Genotype genotype) {
  using enum Genotype;
  switch (genotype) {
    case kGT01:
      UpdateGenotypeRelatedFieldsHelper(record, genotype, kFilteringPassId, {1}, {0.5}, {2});
      break;
    case kGT11:
      UpdateGenotypeRelatedFieldsHelper(record, genotype, kFilteringPassId, {2}, {1.0}, {2});
      break;
    case kGT12:
      UpdateGenotypeRelatedFieldsHelper(record, genotype, kFilteringPassId, {1, 1}, {0.5, 0.5}, {2});
      break;
    case kGT1:
      UpdateGenotypeRelatedFieldsHelper(record, genotype, kFilteringPassId, {1}, {1.0}, {1});
      break;
    case kGT0:
      // For consistency across failed haploid records, force GT=1, AC=1, AF=0.
      UpdateGenotypeRelatedFieldsHelper(record, kGT1, kFilteringFailId, {1}, {0}, {1});
      break;
    default:
      // For consistency across failed diploid records, force GT=0/1, AC=1, AF=0.
      // We cannot naively keep the original `GT` and `AC` values. For example, GT=1/2 split records with only one ALT
      // allele would not comply with GATK ValidateVariants.
      UpdateGenotypeRelatedFieldsHelper(record, kGT01, kFilteringFailId, {1}, {0}, {2});
  }
}

/**
 * @brief Fail a VCF record and update relevant fields.
 * @param record VCF record
 */
void FailRecord(const io::VcfRecordPtr& record, const std::string& fail_id, const Genotype& genotype) {
  using enum Genotype;
  UpdateGenotypeRelatedFields(record, genotype);
  record->AddFilter(fail_id);
}

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
                                      const vec<Genotype>& ml_genotypes) {
  using enum Genotype;
  for (size_t i = 0; i < in_records.size(); ++i) {
    if (std::ranges::find(skip_indexes, i) == skip_indexes.end()) {
      const auto& rec = in_records[i];
      UpdateGenotypeRelatedFields(rec, kGT00);
      if (ml_genotypes.at(i) == kGT00) {
        // classified as a false positive
        rec->AddFilter(kFilteringFalsePositiveId);
      } else {
        // classified as a positive, but force it to fail due to multiallele classification conflict
        rec->AddFilter(kFilteringMultialleleConflictId);
      }
      out_records.emplace_back(rec);
    }
  }
}

void ReconcileGermlineDiploidSingleRecord(const vec<VariantId>& vids,
                                          const vec<io::VcfRecordPtr>& in_records,
                                          const vec<Genotype>& ml_genotypes,
                                          const vec<io::VcfRecordPtr>& alt_wildcard_records,
                                          vec<io::VcfRecordPtr>& out_records,
                                          std::optional<s64>& gt12_01_max_ref_pos) {
  using enum Genotype;
  const auto& record = in_records.at(0);
  const auto& genotype = ml_genotypes.at(0);
  if ((genotype == kGT12 || genotype == kGT11) && !alt_wildcard_records.empty()) {
    // CASE 1(a): GT=1/2 or GT=1/1 and a "wildcard" variant.
    // Wildcard indicates an upstream deletion overlapping this position.
    // GT=0/1 is not allowed because it implies that both REF and ALT alleles have ~0.5 AF.
    // GT=1/2 or GT=1/1 implies that REF allele should have close to ~0.0 AF.
    UpdateGenotypeRelatedFields(record, kGT12);
    MergeMultiAllelicRecords(record, alt_wildcard_records.at(0));
    gt12_01_max_ref_pos = vids.at(0).GetMaxRefAllelePos();
  } else if (genotype == kGT12) {
    // CASE 1(b): GT=1/2 only
    // GT=1/2 is invalid without a partner variant allele.
    // Wildcard can be ignored and not written.
    UpdateGenotypeRelatedFields(record, kGT00);
    record->AddFilter(kFilteringMultiallelePartnerId);
  } else {
    // CASE 1(c): either GT=0/0, GT=0/1, or GT=1/1.
    UpdateGenotypeRelatedFields(record, genotype);
    if (genotype == kGT00) {
      record->AddFilter(kFilteringFalsePositiveId);
    } else if (genotype == kGT01) {
      gt12_01_max_ref_pos = vids.at(0).GetMaxRefAllelePos();
    }
  }
  out_records.emplace_back(record);
}

void ReconcileGermlineDiploidTwoPassingGT12Records(const vec<VariantId>& vids,
                                                   const vec<io::VcfRecordPtr>& in_records,
                                                   const vec<Genotype>& ml_genotypes,
                                                   vec<io::VcfRecordPtr>& out_records,
                                                   std::optional<s64>& gt12_01_max_ref_pos,
                                                   vec<size_t>& gt12_indexes) {
  using enum Genotype;
  // CASE 2: Two passing variants, both GT=1/2
  auto idx1 = gt12_indexes.at(0);
  auto idx2 = gt12_indexes.at(1);
  // for consistency, use the shorter ALT for `1` and the longer ALT for `2`
  if (HasLongerAlt(vids.at(idx1).ref, vids.at(idx1).alt, vids.at(idx2).ref, vids.at(idx2).alt)) {
    std::swap(idx1, idx2);
  }
  const auto& record1 = in_records.at(idx1);
  const auto& record2 = in_records.at(idx2);
  if (MergeMultiAllelicRecords(record1, record2)) {
    out_records.emplace_back(record1);
    gt12_01_max_ref_pos = std::max(vids.at(idx1).GetMaxRefAllelePos(), vids.at(idx2).GetMaxRefAllelePos());
  } else {
    // REF,ALT for the two indels cannot be formatted
    FailRecord(record1, kFilteringMultialleleFormatId, kGT00);
    FailRecord(record2, kFilteringMultialleleFormatId, kGT00);
    out_records.emplace_back(record1);
    out_records.emplace_back(record2);
  }
  // Since there are two variants with GT=1/2, there should be no other variants because diploid should only have
  // two alleles at each position. If other variants happen to exist, fail them.
  FailAndAddGermlineDiploidRecords(in_records, out_records, gt12_indexes, ml_genotypes);
}

void ReconcileGermlineDiploidSinglePassingFromMultiple(const vec<VariantId>& vids,
                                                       const vec<io::VcfRecordPtr>& in_records,
                                                       const vec<Genotype>& ml_genotypes,
                                                       vec<io::VcfRecordPtr>& out_records,
                                                       std::optional<s64>& gt12_01_max_ref_pos,
                                                       vec<size_t>& gt01_11_indexes) {
  using enum Genotype;
  // CASE 3: Single passing variant, either GT=0/1 or GT=1/1
  const auto idx = gt01_11_indexes.at(0);
  const auto& record = in_records.at(idx);
  const auto& score = ml_genotypes.at(idx);
  UpdateGenotypeRelatedFields(record, score);
  out_records.emplace_back(record);
  if (score == kGT01) {
    gt12_01_max_ref_pos = vids.at(idx).GetMaxRefAllelePos();
  }
  // fail other variants
  FailAndAddGermlineDiploidRecords(in_records, out_records, gt01_11_indexes, ml_genotypes);
}

void ReconcileGermlineDiploidTwoPassingMix(const vec<VariantId>& vids,
                                           const vec<io::VcfRecordPtr>& in_records,
                                           const vec<Genotype>& ml_genotypes,
                                           vec<io::VcfRecordPtr>& out_records,
                                           std::optional<s64>& gt12_01_max_ref_pos,
                                           vec<size_t>& gt01_12_indexes) {
  using enum Genotype;
  // CASE 4: Two passing variants, a deletion and an SNV/insertion, both with either GT=0/1 or GT=1/2
  auto idx1 = gt01_12_indexes.at(0);
  auto idx2 = gt01_12_indexes.at(1);
  // for consistency, use the shorter ALT for `1` and the longer ALT for `2`
  if (HasLongerAlt(vids.at(idx1).ref, vids.at(idx1).alt, vids.at(idx2).ref, vids.at(idx2).alt)) {
    std::swap(idx1, idx2);
  }
  const auto& record1 = in_records.at(idx1);
  const auto& record2 = in_records.at(idx2);
  MergeMultiAllelicRecords(record1, record2);
  out_records.emplace_back(record1);
  gt12_01_max_ref_pos = std::max(vids.at(idx1).GetMaxRefAllelePos(), vids.at(idx2).GetMaxRefAllelePos());
  // There should be no other variants because a diploid region is expected to have only 2 variants at a position.
  // If other variants happened to exist, fail them.
  FailAndAddGermlineDiploidRecords(in_records, out_records, gt01_12_indexes, ml_genotypes);
}

std::optional<s64> ReconcileGermlineDiploidRecords(const vec<VariantId>& vids,
                                                   const vec<io::VcfRecordPtr>& in_records,
                                                   const vec<io::VcfRecordPtr>& alt_wildcard_records,
                                                   const vec<Genotype>& ml_genotypes,
                                                   vec<io::VcfRecordPtr>& out_records) {
  // Note that `out_records` is a vector for output VCF records for the region being processed in parallel.
  // Reconciled records will be appended to this vector.

  // This variable stores the maximum REF allele position for `0/1` or `1/2` genotypes only.
  // It will be returned for reconciling wildcard ALT alleles at downstream positions.
  std::optional<s64> gt12_01_max_ref_pos;

  // Store the index of each variant for different genotype groups:
  // - GT=0/1 and GT=1/2:
  //   - for updating `gt12_01_max_ref_pos`
  //   - for rescuing misclassified GT=1/2 records involving a deletion and an SNV/insertion
  // - GT=1/2 only:
  //   - for merging unambiguous multi-allelic records
  // - GT=0/1 and GT=1/1:
  //   - for checking unambiguous `0/1` or `1/1` genotypes
  vec<size_t> gt01_12_indexes;
  vec<size_t> gt12_indexes;
  vec<size_t> gt01_11_indexes;
  u32 num_fail = 0;
  vec<bool> is_snv_ins_vec;
  const size_t num_variants = vids.size();
  // We should expect to see at most 2 variants at a position because this is a diploid region.
  // If not then this is likely an ill-formatted input VCF record. We will still process this position and let the
  // outcome of genotype reconciliation dictate the output VCF records.
  GetGTIndexes(vids, gt01_12_indexes, is_snv_ins_vec, gt12_indexes, num_fail, gt01_11_indexes, ml_genotypes);

  if (num_variants == 1) {
    // CASE 1: only 1 variant at this position
    ReconcileGermlineDiploidSingleRecord(
        vids, in_records, ml_genotypes, alt_wildcard_records, out_records, gt12_01_max_ref_pos);
  } else {
    // multiple variants at this position
    const auto num_gt12 = gt12_indexes.size();
    const auto num_gt01_11 = gt01_11_indexes.size();
    const auto num_pass = num_variants - num_fail;

    bool snv_ins_split = false;
    if (num_pass == 2 && gt01_12_indexes.size() == 2 && !is_snv_ins_vec.empty()) {
      snv_ins_split = ((is_snv_ins_vec.at(gt01_12_indexes.at(0)) && !is_snv_ins_vec.at(gt01_12_indexes.at(1))) ||
                       (!is_snv_ins_vec.at(gt01_12_indexes.at(0)) && is_snv_ins_vec.at(gt01_12_indexes.at(1))));
    }

    if (num_pass == 2 && num_gt12 == 2) {
      // CASE 2: Two passing variants, both GT=1/2
      ReconcileGermlineDiploidTwoPassingGT12Records(
          vids, in_records, ml_genotypes, out_records, gt12_01_max_ref_pos, gt12_indexes);
    } else if (num_pass == 1 && num_gt01_11 == 1) {
      // CASE 3: Single passing variant, either GT=0/1 or GT=1/1
      ReconcileGermlineDiploidSinglePassingFromMultiple(
          vids, in_records, ml_genotypes, out_records, gt12_01_max_ref_pos, gt01_11_indexes);
    } else if (num_pass == 2 && gt01_12_indexes.size() == 2 && snv_ins_split) {
      // CASE 4: Two passing variants, a deletion and an SNV/insertion, both with either GT=0/1 or GT=1/2
      // This is intended to rescue ML misclassification due to ref-allele features being offset by 1 between deletions
      // and SNV/insertion.
      // Only GT=0/1 or GT=1/2 are allowed because they imply half of the read depth supports the ALT allele.
      ReconcileGermlineDiploidTwoPassingMix(
          vids, in_records, ml_genotypes, out_records, gt12_01_max_ref_pos, gt01_12_indexes);
    } else {
      // CASE 5: GT=0/0 or conflicting genotypes.
      FailAndAddGermlineDiploidRecords(in_records, out_records, {}, ml_genotypes);
    }
  }

  return gt12_01_max_ref_pos;
}

void ReconcileGermlineHaploidRecords(const vec<VariantId>& vids,
                                     const vec<io::VcfRecordPtr>& in_records,
                                     const vec<Genotype>& ml_genotypes,
                                     vec<io::VcfRecordPtr>& out_records) {
  // Note that `out_records` is a vector for output VCF records for the region being processed in parallel.
  // Reconciled records will be appended to this vector.
  // Unlike ReconcileGermlineDiploidRecords, we do not track the max REF allele position because haploid regions are
  // not expected to have multiple variants at the same position.

  // find the index of the GT=1/1 record, if it exists
  const size_t num_variants = vids.size();
  vec<size_t> gt11_indexes;
  for (size_t i = 0; i < vids.size(); ++i) {
    const auto& genotype = ml_genotypes.at(i);
    if (genotype == Genotype::kGT11) {
      gt11_indexes.emplace_back(i);
    }
  }

  if (gt11_indexes.size() == 1) {
    // GT=1/1 diploid prediction is treated as GT=1 in haploid region.
    // Other genotypes are failed.
    const auto index = gt11_indexes.at(0);
    const auto& record = in_records.at(index);
    UpdateGenotypeRelatedFields(record, Genotype::kGT1);
    out_records.emplace_back(record);
    if (num_variants > 1) {
      // fail all other variants
      for (size_t i = 0; i < in_records.size(); ++i) {
        if (i != index) {
          const auto& rec = in_records.at(i);
          FailRecord(rec, kFilteringFalsePositiveId, Genotype::kGT0);
          out_records.emplace_back(rec);
        }
      }
    }
  } else {
    // fail all variants
    for (const auto& rec : in_records) {
      FailRecord(rec, kFilteringFalsePositiveId, Genotype::kGT0);
      out_records.emplace_back(rec);
    }
  }
}
}  // namespace xoos::svc
