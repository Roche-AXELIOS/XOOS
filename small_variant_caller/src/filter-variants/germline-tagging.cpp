#include "germline-tagging.h"

#include "core/filtering.h"
#include "core/vcf-fields.h"

namespace xoos::svc {

void UpdateGermlineTaggingRecord(const io::VcfRecordPtr& record,
                                 const Genotype genotype,
                                 const GenotypeScore& score,
                                 const s32 sample_idx) {
  // Number of samples; tumor and normal samples
  static constexpr size_t kNumSamples = 2;
  if (sample_idx < 0) {
    return;
  }
  // update GT field
  record->SetGTField(GenotypeToString(genotype), sample_idx);
  // update QUAL
  record->SetQuality(score.variant_quality);
  // update FILTER, INFO, and FORMAT fields
  switch (genotype) {
    case Genotype::kGT01:
    case Genotype::kGT11: {
      record->SetFilter(kFilteringPassId);
      record->AddInfoFieldFlag(kGermlineTaggingInfoGermlineId);
      // remove the SOMATIC INFO field, if present
      record->RemoveInfoFieldFlag(kGermlineTaggingInfoSomaticId);
      break;
    }
    default: {
      record->SetFilter(kFilteringFailId);
      // remove the GERMLINE and SOMATIC INFO fields, if present
      record->RemoveInfoFieldFlag(kGermlineTaggingInfoGermlineId);
      record->RemoveInfoFieldFlag(kGermlineTaggingInfoSomaticId);
    }
  }
  // update GQ field
  auto gq_vals = record->GetFormatFieldNoCheck<s32>(kFieldGq);
  if (gq_vals.size() != kNumSamples) {
    gq_vals = vec<s32>(kNumSamples, 0);
  }
  gq_vals[ToUnsigned(sample_idx)] = score.genotype_quality;
  record->SetFormatField<s32>(kFieldGq, gq_vals);
  // update ML field
  auto scores = record->GetFormatFieldNoCheck<float>(kMachineLearningId);
  if (scores.size() != kNumSamples) {
    scores = vec<f32>(kNumSamples, 0.0f);
  }
  scores[ToUnsigned(sample_idx)] = static_cast<f32>(score.probability);
  record->SetFormatField(kMachineLearningId, scores);
}

void ReconcileGermlineTaggingRecords(const vec<io::VcfRecordPtr>& in_records,
                                     const vec<GenotypeScore>& ml_genotypes,
                                     const bool is_haploid,
                                     const s32 normal_sample_idx,
                                     vec<io::VcfRecordPtr>& out_records) {
  using enum Genotype;
  // identify indexes of passing variants
  vec<size_t> pass_idx;
  const Genotype fail_genotype = is_haploid ? kGT0 : kGT00;
  for (size_t i = 0; i < ml_genotypes.size(); ++i) {
    const auto genotype = ml_genotypes.at(i).genotype;
    if ((!is_haploid && genotype == kGT01) || genotype == kGT11) {
      // pass GT=0/1 or GT=1/1 in diploid region
      // pass only GT=1/1 in haploid region
      pass_idx.emplace_back(i);
    }
  }
  if (pass_idx.size() == 1) {
    // pass the variant
    const size_t idx = pass_idx.at(0);
    // convert GT=1/1 to GT=1 in haploid region
    const auto pass_genotype = is_haploid ? kGT1 : ml_genotypes.at(idx).genotype;
    UpdateGermlineTaggingRecord(in_records.at(idx), pass_genotype, ml_genotypes.at(idx), normal_sample_idx);
    out_records.emplace_back(in_records.at(idx));
    // fail other variants
    for (size_t i = 0; i < in_records.size(); ++i) {
      if (i != idx) {
        UpdateGermlineTaggingRecord(in_records.at(i), fail_genotype, ml_genotypes.at(i), normal_sample_idx);
        out_records.emplace_back(in_records.at(i));
      }
    }
  } else {
    // fail all variants
    for (size_t i = 0; i < in_records.size(); ++i) {
      UpdateGermlineTaggingRecord(in_records.at(i), fail_genotype, ml_genotypes.at(i), normal_sample_idx);
      out_records.emplace_back(in_records.at(i));
    }
  }
}

}  // namespace xoos::svc
