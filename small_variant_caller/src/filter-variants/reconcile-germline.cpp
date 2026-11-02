#include "reconcile-germline.h"

#include <algorithm>
#include <unordered_map>
#include <utility>

#include <xoos/error/error.h>
#include <xoos/io/bed-region.h>
#include <xoos/io/fasta-reader.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/log/logging.h>
#include <xoos/types/vec.h>

#include "core/filtering.h"
#include "core/vcf-fields.h"
#include "util/seq-util.h"
#include "util/vcf-util.h"

namespace xoos::svc {

void GetGTIndexes(const vec<VariantId>& vids,
                  vec<size_t>& gt01_12_indexes,
                  std::deque<bool>& is_snv_ins_vec,
                  vec<size_t>& gt12_indexes,
                  u32& num_fail,
                  vec<size_t>& gt01_11_indexes,
                  const vec<PredictionScore>& ml_genotypes) {
  // Get the indexes of the variants based on their genotypes and store them in the relevant vectors. Needed for
  // downstream filtering
  using enum VariantType;
  using enum Genotype;
  const auto num_variants = vids.size();
  for (size_t i = 0; i < num_variants; ++i) {
    const auto& vid = vids.at(i);
    is_snv_ins_vec.emplace_back(vid.type == kSNV || vid.type == kInsertion);
    switch (ml_genotypes.at(i).genotype) {
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

namespace {
void UpdateGenotypeRelatedFieldsHelper(const io::VcfRecordPtr& record,
                                       const Genotype genotype,
                                       const std::string& filter,
                                       const vec<s32>& ac,
                                       const vec<f32>& af,
                                       const vec<s32>& an) {
  record->AddFilter(filter);
  record->SetInfoField(kFieldAc, ac);
  record->SetInfoField(kFieldAf, af);
  record->SetInfoField(kFieldAn, an);
  record->SetGTField(GenotypeToString(genotype));
}
}  // namespace

void UpdateGenotypeRelatedFields(const io::VcfRecordPtr& record, const Genotype genotype) {
  static const vec<s32> kDiploidAN{2};
  static const vec<s32> kHaploidAN{1};
  constexpr s32 kHomAltAlleleCount = 2;
  // Expected allele frequency for a heterozygous genotype.
  constexpr f64 kHetAlleleFreq{0.5};
  using enum Genotype;
  switch (genotype) {
    case kGT01:
      UpdateGenotypeRelatedFieldsHelper(record, genotype, kFilteringPassId, {1}, {kHetAlleleFreq}, kDiploidAN);
      break;
    case kGT11:
      UpdateGenotypeRelatedFieldsHelper(record, genotype, kFilteringPassId, {kHomAltAlleleCount}, {1.0}, kDiploidAN);
      break;
    case kGT12:
      UpdateGenotypeRelatedFieldsHelper(
          record, genotype, kFilteringPassId, {1, 1}, {kHetAlleleFreq, kHetAlleleFreq}, kDiploidAN);
      break;
    case kGT1:
      UpdateGenotypeRelatedFieldsHelper(record, genotype, kFilteringPassId, {1}, {1.0}, kHaploidAN);
      break;
    case kGT0:
      // For consistency across failed haploid records, force GT=., AC=1, AF=0.
      UpdateGenotypeRelatedFieldsHelper(record, kGTN, kFilteringFailId, {1}, {0}, kHaploidAN);
      break;
    default:
      // For consistency across failed diploid records, force GT=./., AC=1, AF=0.
      // We cannot naively keep the original `GT` and `AC` values. For example, GT=1/2 split records with only one ALT
      // allele would not comply with GATK ValidateVariants.
      UpdateGenotypeRelatedFieldsHelper(record, kGTNA, kFilteringFailId, {1}, {0}, kDiploidAN);
  }
}

/**
 * @brief Fail a VCF record and update relevant fields.
 * @param record VCF record
 * @param fail_id Filter ID to be added
 * @param genotype Genotype to be set
 */
void FailRecord(const io::VcfRecordPtr& record, const std::string& fail_id, const Genotype& genotype) {
  UpdateGenotypeRelatedFields(record, genotype);
  record->AddFilter(fail_id);
}

namespace {
/**
 * @brief Helper function to update the ML score, genotype quality, and variant quality fields of a VCF record based on
 * the predicted genotype score.
 * @param record VCF record to be updated
 * @param score PredictionScore struct containing the genotype, ML score, genotype quality, and variant quality
 */
void UpdateScoreRelatedFields(const io::VcfRecordPtr& record, const PredictionScore& score) {
  record->SetFormatField(kMachineLearningId, vec<f32>{static_cast<f32>(score.probability)});
  record->SetFormatField(kFieldGq, vec<s32>{score.genotype_quality});
  record->SetQuality(score.variant_quality);
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
                                      const vec<PredictionScore>& ml_genotypes) {
  using enum Genotype;
  for (size_t i = 0; i < in_records.size(); ++i) {
    if (std::ranges::find(skip_indexes, i) == skip_indexes.end()) {
      const auto& rec = in_records[i];
      UpdateGenotypeRelatedFields(rec, kGT00);
      if (ml_genotypes.at(i).genotype == kGT00) {
        // classified as a false positive
        rec->AddFilter(kFilteringFalsePositiveId);
      } else {
        // classified as a positive, but force it to fail due to multiallele classification conflict
        rec->AddFilter(kFilteringMultialleleConflictId);
      }
      // update score related fields for all failed records
      UpdateScoreRelatedFields(rec, ml_genotypes.at(i));
      out_records.emplace_back(rec);
    }
  }
}

/**
 * @brief Helper function to update the ML score, genotype quality, and variant quality fields of a VCF record based
 * on two predicted genotype scores. Used for cases where there are two GT=1/2 records at the same position, and we
 * need to determine which one to use for downstream filtering.
 * @pre At least one of the two predicted genotypes must be classified as GT=1/2, otherwise there would be a
 * classification conflict that should have been resolved in upstream filtering steps.
 * @param record VCF record to be updated
 * @param score1 First predicted genotype score
 * @param score2 Second predicted genotype score
 */
void UpdateScoreRelatedFields(const io::VcfRecordPtr& record,
                              const PredictionScore& score1,
                              const PredictionScore& score2) {
  const bool is_score1_gt12 = score1.genotype == Genotype::kGT12;
  const bool is_score2_gt12 = score2.genotype == Genotype::kGT12;
  // If only one of the two alleles is classified as GT=1/2, we use the PredictionScore from that allele.
  if (is_score1_gt12 && !is_score2_gt12) {
    UpdateScoreRelatedFields(record, score1);
    return;
  }
  if (!is_score1_gt12 && is_score2_gt12) {
    UpdateScoreRelatedFields(record, score2);
    return;
  }
  // If both alleles are classified as GT=1/2, we use the PredictionScore with the lower prediction probability, which
  // is more conservative and has a lower chance of overestimating the variant quality.
  UpdateScoreRelatedFields(record, score1.probability < score2.probability ? score1 : score2);
}
}  // namespace

void ReconcileGermlineDiploidSingleRecord(const vec<VariantId>& vids,
                                          const vec<io::VcfRecordPtr>& in_records,
                                          const vec<PredictionScore>& ml_genotypes,
                                          const vec<io::VcfRecordPtr>& alt_wildcard_records,
                                          vec<io::VcfRecordPtr>& out_records,
                                          std::optional<s64>& gt12_01_max_ref_pos,
                                          const vec<io::InfoFieldMetadata>& info_metadata,
                                          const vec<io::FormatFieldMetadata>& fmt_metadata) {
  using enum Genotype;
  const auto& record = in_records.at(0);
  const auto& genotype = ml_genotypes.at(0).genotype;
  if ((genotype == kGT12 || genotype == kGT11) && !alt_wildcard_records.empty()) {
    // CASE 1(a): GT=1/2 or GT=1/1 and a "wildcard" variant.
    // Wildcard indicates an upstream deletion overlapping this position.
    // GT=0/1 is not allowed because it implies that both REF and ALT alleles have ~0.5 AF.
    // GT=1/2 or GT=1/1 implies that REF allele should have close to ~0.0 AF.
    MergeMultiAllelicRecords(record, alt_wildcard_records.at(0), info_metadata, fmt_metadata);
    record->AddFilter(kFilteringPassId);
    gt12_01_max_ref_pos = vids.at(0).GetMaxRefAllelePos();
  } else if (genotype == kGT12) {
    // CASE 1(b): GT=1/2 without partner allele.
    // Use afq_gt (AF-based genotype from BAM features) to rescue or fail.
    const auto afq_gt = ml_genotypes.at(0).afq_gt;
    if (afq_gt == kGT00 || afq_gt == kGT12) {
      UpdateGenotypeRelatedFields(record, kGT00);
      record->AddFilter(kFilteringMultiallelePartnerId);
    } else {
      UpdateGenotypeRelatedFields(record, afq_gt);
      if (afq_gt == kGT01) {
        gt12_01_max_ref_pos = vids.at(0).GetMaxRefAllelePos();
      }
    }
  } else {
    // CASE 1(c): either GT=0/0, GT=0/1, or GT=1/1.
    UpdateGenotypeRelatedFields(record, genotype);
    if (genotype == kGT00) {
      record->AddFilter(kFilteringFalsePositiveId);
    } else {
      if (genotype == kGT01) {
        gt12_01_max_ref_pos = vids.at(0).GetMaxRefAllelePos();
      }
    }
  }
  UpdateScoreRelatedFields(record, ml_genotypes.at(0));
  out_records.emplace_back(record);
}

void ReconcileGermlineDiploidTwoPassingGT12Records(const vec<VariantId>& vids,
                                                   const vec<io::VcfRecordPtr>& in_records,
                                                   const vec<PredictionScore>& ml_genotypes,
                                                   vec<io::VcfRecordPtr>& out_records,
                                                   std::optional<s64>& gt12_01_max_ref_pos,
                                                   const vec<size_t>& gt12_indexes,
                                                   const vec<io::InfoFieldMetadata>& info_metadata,
                                                   const vec<io::FormatFieldMetadata>& fmt_metadata) {
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
  if (MergeMultiAllelicRecords(record1, record2, info_metadata, fmt_metadata)) {
    // Record 2 is merged into record 1
    record1->AddFilter(kFilteringPassId);
    UpdateScoreRelatedFields(record1, ml_genotypes.at(idx1), ml_genotypes.at(idx2));
    out_records.emplace_back(record1);
    gt12_01_max_ref_pos = std::max(vids.at(idx1).GetMaxRefAllelePos(), vids.at(idx2).GetMaxRefAllelePos());
  } else {
    // REF,ALT for the two indels cannot be formatted
    FailRecord(record1, kFilteringMultialleleFormatId, kGT00);
    FailRecord(record2, kFilteringMultialleleFormatId, kGT00);
    UpdateScoreRelatedFields(record1, ml_genotypes.at(idx1));
    UpdateScoreRelatedFields(record2, ml_genotypes.at(idx2));
    out_records.emplace_back(record1);
    out_records.emplace_back(record2);
  }
  // Since there are two variants with GT=1/2, there should be no other variants because diploid should only have
  // two alleles at each position. If other variants happen to exist, fail them.
  FailAndAddGermlineDiploidRecords(in_records, out_records, gt12_indexes, ml_genotypes);
}

void ReconcileGermlineDiploidSinglePassingFromMultiple(const vec<VariantId>& vids,
                                                       const vec<io::VcfRecordPtr>& in_records,
                                                       const vec<PredictionScore>& ml_genotypes,
                                                       vec<io::VcfRecordPtr>& out_records,
                                                       std::optional<s64>& gt12_01_max_ref_pos,
                                                       const vec<size_t>& gt01_11_indexes) {
  using enum Genotype;
  // CASE 3: Single passing variant, either GT=0/1 or GT=1/1
  const auto idx = gt01_11_indexes.at(0);
  const auto& record = in_records.at(idx);
  const auto& genotype = ml_genotypes.at(idx).genotype;
  UpdateGenotypeRelatedFields(record, genotype);
  UpdateScoreRelatedFields(record, ml_genotypes.at(idx));
  out_records.emplace_back(record);
  if (genotype == kGT01) {
    gt12_01_max_ref_pos = vids.at(idx).GetMaxRefAllelePos();
  }
  // fail other variants
  FailAndAddGermlineDiploidRecords(in_records, out_records, gt01_11_indexes, ml_genotypes);
}

void ReconcileGermlineDiploidTwoPassingMix(const vec<VariantId>& vids,
                                           const vec<io::VcfRecordPtr>& in_records,
                                           const vec<PredictionScore>& ml_genotypes,
                                           vec<io::VcfRecordPtr>& out_records,
                                           std::optional<s64>& gt12_01_max_ref_pos,
                                           const vec<size_t>& gt01_12_indexes,
                                           const vec<io::InfoFieldMetadata>& info_metadata,
                                           const vec<io::FormatFieldMetadata>& fmt_metadata) {
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
  MergeMultiAllelicRecords(record1, record2, info_metadata, fmt_metadata);
  record1->AddFilter(kFilteringPassId);
  UpdateScoreRelatedFields(record1, ml_genotypes.at(idx1), ml_genotypes.at(idx2));
  out_records.emplace_back(record1);
  gt12_01_max_ref_pos = std::max(vids.at(idx1).GetMaxRefAllelePos(), vids.at(idx2).GetMaxRefAllelePos());
  // There should be no other variants because a diploid region is expected to have only 2 variants at a position.
  // If other variants happened to exist, fail them.
  FailAndAddGermlineDiploidRecords(in_records, out_records, gt01_12_indexes, ml_genotypes);
}

bool ReconcileGermlineDiploidAfqGtRescue(const vec<VariantId>& vids,
                                         const vec<io::VcfRecordPtr>& in_records,
                                         const vec<PredictionScore>& ml_genotypes,
                                         vec<io::VcfRecordPtr>& out_records,
                                         std::optional<s64>& gt12_01_max_ref_pos,
                                         const vec<io::InfoFieldMetadata>& info_metadata,
                                         const vec<io::FormatFieldMetadata>& fmt_metadata) {
  using enum Genotype;
  constexpr size_t kDiploidVariantCount = 2;

  // Group variant indexes by their AF-based afq_gt call.
  vec<size_t> afq_gt01_11_indexes;
  vec<size_t> afq_gt12_indexes;
  for (size_t i = 0; i < vids.size(); ++i) {
    const auto afq_gt = ml_genotypes.at(i).afq_gt;
    if (afq_gt == kGT01 || afq_gt == kGT11) {
      afq_gt01_11_indexes.emplace_back(i);
    } else if (afq_gt == kGT12) {
      afq_gt12_indexes.emplace_back(i);
    }
  }

  if (afq_gt01_11_indexes.size() == 1 && afq_gt12_indexes.empty()) {
    // CASE 5a: Exactly one variant with afq_gt=0/1 or 1/1 — rescue like CASE 3.
    const auto idx = afq_gt01_11_indexes.at(0);
    const auto& record = in_records.at(idx);
    const auto afq_gt = ml_genotypes.at(idx).afq_gt;
    UpdateGenotypeRelatedFields(record, afq_gt);
    UpdateScoreRelatedFields(record, ml_genotypes.at(idx));
    out_records.emplace_back(record);
    if (afq_gt == kGT01) {
      gt12_01_max_ref_pos = vids.at(idx).GetMaxRefAllelePos();
    }
    FailAndAddGermlineDiploidRecords(in_records, out_records, afq_gt01_11_indexes, ml_genotypes);
    return true;
  }

  if (afq_gt12_indexes.size() == kDiploidVariantCount && afq_gt01_11_indexes.empty()) {
    // CASE 5b: Exactly two variants with afq_gt=1/2 — merge like CASE 2.
    auto idx1 = afq_gt12_indexes.at(0);
    auto idx2 = afq_gt12_indexes.at(1);
    if (HasLongerAlt(vids.at(idx1).ref, vids.at(idx1).alt, vids.at(idx2).ref, vids.at(idx2).alt)) {
      std::swap(idx1, idx2);
    }
    const auto& record1 = in_records.at(idx1);
    const auto& record2 = in_records.at(idx2);
    if (MergeMultiAllelicRecords(record1, record2, info_metadata, fmt_metadata)) {
      record1->AddFilter(kFilteringPassId);
      UpdateScoreRelatedFields(record1, ml_genotypes.at(idx1), ml_genotypes.at(idx2));
      out_records.emplace_back(record1);
      gt12_01_max_ref_pos = std::max(vids.at(idx1).GetMaxRefAllelePos(), vids.at(idx2).GetMaxRefAllelePos());
    } else {
      FailRecord(record1, kFilteringMultialleleFormatId, kGT00);
      FailRecord(record2, kFilteringMultialleleFormatId, kGT00);
      UpdateScoreRelatedFields(record1, ml_genotypes.at(idx1));
      UpdateScoreRelatedFields(record2, ml_genotypes.at(idx2));
      out_records.emplace_back(record1);
      out_records.emplace_back(record2);
    }
    FailAndAddGermlineDiploidRecords(in_records, out_records, afq_gt12_indexes, ml_genotypes);
    return true;
  }

  return false;
}

std::optional<s64> ReconcileGermlineDiploidRecords(const vec<VariantId>& vids,
                                                   const vec<io::VcfRecordPtr>& in_records,
                                                   const vec<io::VcfRecordPtr>& alt_wildcard_records,
                                                   const vec<PredictionScore>& ml_genotypes,
                                                   vec<io::VcfRecordPtr>& out_records,
                                                   const VcfFieldMetadata& field_metadata,
                                                   PerAlleleFieldSnapshot& span_deletion_fields) {
  const vec<io::InfoFieldMetadata>& info_metadata = field_metadata.info;
  const vec<io::FormatFieldMetadata>& fmt_metadata = field_metadata.format;
  // Note that `out_records` is a vector for output VCF records for the region being processed in parallel.
  // Reconciled records will be appended to this vector.

  // This variable stores the maximum REF allele position for `0/1` or `1/2` genotypes only.
  // It will be returned for reconciling wildcard ALT alleles at downstream positions.
  std::optional<s64> gt12_01_max_ref_pos;

  // Snapshot the per-allele FORMAT fields of every deletion record before dispatch, keyed by the reference position its
  // REF allele spans to. Only a deletion spans downstream positions, so only a deletion can later appear as a wildcard
  // (`*`) allele. Snapshots are taken here (pre-merge) because the CASE handlers below mutate records in-place when
  // merging multi-allelic records, which would otherwise overwrite the deletion's original bi-allelic counts.
  std::unordered_map<s64, PerAlleleFieldSnapshot> deletion_fields_by_max_ref_pos;
  for (size_t i = 0; i < vids.size(); ++i) {
    if (vids.at(i).type == VariantType::kDeletion) {
      deletion_fields_by_max_ref_pos.try_emplace(static_cast<s64>(vids.at(i).GetMaxRefAllelePos()),
                                                 SnapshotPerAlleleFormatFields(in_records.at(i), fmt_metadata));
    }
  }

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
  std::deque<bool> is_snv_ins_vec;
  const size_t num_variants = vids.size();
  // We should expect to see at most 2 variants at a position because this is a diploid region.
  // If not then this is likely an ill-formatted input VCF record. We will still process this position and let the
  // outcome of genotype reconciliation dictate the output VCF records.
  GetGTIndexes(vids, gt01_12_indexes, is_snv_ins_vec, gt12_indexes, num_fail, gt01_11_indexes, ml_genotypes);

  if (num_variants == 1) {
    // CASE 1: only 1 variant at this position
    ReconcileGermlineDiploidSingleRecord(vids,
                                         in_records,
                                         ml_genotypes,
                                         alt_wildcard_records,
                                         out_records,
                                         gt12_01_max_ref_pos,
                                         info_metadata,
                                         fmt_metadata);
  } else {
    // multiple variants at this position
    const auto num_gt12 = gt12_indexes.size();
    const auto num_gt01_11 = gt01_11_indexes.size();
    const auto num_pass = num_variants - num_fail;

    constexpr size_t kDiploidVariantCount = 2;
    bool snv_ins_split = false;
    if (num_pass == kDiploidVariantCount && gt01_12_indexes.size() == kDiploidVariantCount && !is_snv_ins_vec.empty()) {
      snv_ins_split = ((is_snv_ins_vec.at(gt01_12_indexes.at(0)) && !is_snv_ins_vec.at(gt01_12_indexes.at(1))) ||
                       (!is_snv_ins_vec.at(gt01_12_indexes.at(0)) && is_snv_ins_vec.at(gt01_12_indexes.at(1))));
    }

    if (num_pass == kDiploidVariantCount && num_gt12 == kDiploidVariantCount) {
      // CASE 2: Two passing variants, both GT=1/2
      ReconcileGermlineDiploidTwoPassingGT12Records(
          vids, in_records, ml_genotypes, out_records, gt12_01_max_ref_pos, gt12_indexes, info_metadata, fmt_metadata);
    } else if (num_pass == 1 && num_gt01_11 == 1) {
      // CASE 3: Single passing variant, either GT=0/1 or GT=1/1
      ReconcileGermlineDiploidSinglePassingFromMultiple(
          vids, in_records, ml_genotypes, out_records, gt12_01_max_ref_pos, gt01_11_indexes);
    } else if (num_pass == kDiploidVariantCount && gt01_12_indexes.size() == kDiploidVariantCount && snv_ins_split) {
      // CASE 4: Two passing variants, a deletion and an SNV/insertion, both with either GT=0/1 or GT=1/2
      // This is intended to rescue ML misclassification due to ref-allele features being offset by 1 between deletions
      // and SNV/insertion.
      // Only GT=0/1 or GT=1/2 are allowed because they imply half of the read depth supports the ALT allele.
      ReconcileGermlineDiploidTwoPassingMix(vids,
                                            in_records,
                                            ml_genotypes,
                                            out_records,
                                            gt12_01_max_ref_pos,
                                            gt01_12_indexes,
                                            info_metadata,
                                            fmt_metadata);
    } else {
      // CASE 5: GT=0/0 or conflicting genotypes.
      // Attempt rescue using afq_gt before failing all variants.
      const bool rescued = ReconcileGermlineDiploidAfqGtRescue(
          vids, in_records, ml_genotypes, out_records, gt12_01_max_ref_pos, info_metadata, fmt_metadata);
      if (!rescued) {
        // No rescue possible — fail all variants.
        FailAndAddGermlineDiploidRecords(in_records, out_records, {}, ml_genotypes);
      }
    }
  }

  // If a spanning deletion determined the returned max reference position, carry its per-allele field snapshot forward
  // so a downstream wildcard (`*`) allele can be emitted with the deletion's computed counts. Otherwise clear it.
  span_deletion_fields = {};
  if (gt12_01_max_ref_pos.has_value()) {
    const auto it = deletion_fields_by_max_ref_pos.find(gt12_01_max_ref_pos.value());
    if (it != deletion_fields_by_max_ref_pos.end()) {
      span_deletion_fields = it->second;
    }
  }

  return gt12_01_max_ref_pos;
}

void ReconcileGermlineHaploidRecords(const vec<VariantId>& vids,
                                     const vec<io::VcfRecordPtr>& in_records,
                                     const vec<PredictionScore>& ml_genotypes,
                                     vec<io::VcfRecordPtr>& out_records) {
  // Note that `out_records` is a vector for output VCF records for the region being processed in parallel.
  // Reconciled records will be appended to this vector.
  // Unlike ReconcileGermlineDiploidRecords, we do not track the max REF allele position because haploid regions are
  // not expected to have multiple variants at the same position.

  // find the index of the GT=1/1 record, if it exists
  const size_t num_variants = vids.size();
  vec<size_t> gt11_indexes;
  for (size_t i = 0; i < vids.size(); ++i) {
    const auto& genotype = ml_genotypes.at(i).genotype;
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
    UpdateScoreRelatedFields(record, ml_genotypes.at(index));
    out_records.emplace_back(record);
    if (num_variants > 1) {
      // fail all other variants
      for (size_t i = 0; i < in_records.size(); ++i) {
        if (i != index) {
          const auto& rec = in_records.at(i);
          FailRecord(rec, kFilteringFalsePositiveId, Genotype::kGT0);
          UpdateScoreRelatedFields(rec, ml_genotypes.at(i));
          out_records.emplace_back(rec);
        }
      }
    }
  } else {
    // fail all variants
    for (size_t i = 0; i < in_records.size(); ++i) {
      const auto& rec = in_records.at(i);
      FailRecord(rec, kFilteringFalsePositiveId, Genotype::kGT0);
      UpdateScoreRelatedFields(rec, ml_genotypes.at(i));
      out_records.emplace_back(rec);
    }
  }
}
}  // namespace xoos::svc
