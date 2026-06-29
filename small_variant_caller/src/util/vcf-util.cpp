#include "vcf-util.h"

#include "core/genotype.h"
#include "core/vcf-fields.h"
#include "seq-util.h"

namespace xoos::svc {

bool IsMultiAllelicRecord(const io::VcfRecordPtr& record) {
  // A multi-allelic record has 3 or more alleles (1 REF and 2+ ALT)
  static const s32 kMinNumAlleles = 3;
  return record->NumAlleles() >= kMinNumAlleles;
}

/**
 * @brief Update a field (INFO or FORMAT) for a specific alternate allele index.
 * @tparam Metadata Metadata type (InfoFieldMetadata or FormatFieldMetadata)
 * @tparam GetFieldFn Function type to get the field values
 * @tparam SetFieldFn Function type to set the field values
 * @param record VCF record
 * @param metadata Field metadata
 * @param alt_idx Alternate allele index (1-based)
 * @param get_field Function to get the field values
 * @param set_field Function to set the field values
 */
template <typename Metadata, typename GetFieldFn, typename SetFieldFn>
static void UpdateFieldForAltAllele(const io::VcfRecordPtr& record,
                                    const Metadata& metadata,
                                    const s32 alt_idx,
                                    GetFieldFn get_field,
                                    SetFieldFn set_field) {
  if (alt_idx < 1) {
    return;
  }
  const auto& values = get_field(record, metadata.id);
  if (metadata.number == io::kNumberEachAllele) {
    const auto val_idx = ToUnsigned(alt_idx - 1);
    if (val_idx < values.size()) {
      set_field(record, metadata.id, {values.at(val_idx)});
    }
  } else if (metadata.number == io::kNumberR) {
    const auto val_idx = ToUnsigned(alt_idx);
    if (val_idx <= values.size()) {
      set_field(record, metadata.id, {values.front(), values.at(val_idx)});
    }
  }
}

/**
 * @brief Update an INFO field for a specific alternate allele index.
 * @tparam T Field value type
 * @param record VCF record
 * @param metadata INFO field metadata
 * @param alt_idx Alternate allele index (1-based)
 */
template <typename T>
static void UpdateInfoFieldForAltAllele(const io::VcfRecordPtr& record,
                                        const io::InfoFieldMetadata& metadata,
                                        const s32 alt_idx) {
  UpdateFieldForAltAllele(
      record,
      metadata,
      alt_idx,
      [](const io::VcfRecordPtr& rec, const std::string& id) { return rec->GetInfoFieldNoCheck<T>(id); },
      [](const io::VcfRecordPtr& rec, const std::string& id, const vec<T>& vals) { rec->SetInfoField<T>(id, vals); });
}

/**
 * @brief Update a FORMAT field for a specific alternate allele index.
 * @tparam T Field value type
 * @param record VCF record
 * @param metadata FORMAT field metadata
 * @param alt_idx Alternate allele index (1-based)
 */
template <typename T>
static void UpdateFormatFieldForAltAllele(const io::VcfRecordPtr& record,
                                          const io::FormatFieldMetadata& metadata,
                                          const s32 alt_idx) {
  UpdateFieldForAltAllele(
      record,
      metadata,
      alt_idx,
      [](const io::VcfRecordPtr& rec, const std::string& id) { return rec->GetFormatFieldNoCheck<T>(id); },
      [](const io::VcfRecordPtr& rec, const std::string& id, const vec<T>& vals) { rec->SetFormatField<T>(id, vals); });
}

template void UpdateInfoFieldForAltAllele<s32>(const io::VcfRecordPtr& record,
                                               const io::InfoFieldMetadata& metadata,
                                               s32 alt_idx);
template void UpdateInfoFieldForAltAllele<f32>(const io::VcfRecordPtr& record,
                                               const io::InfoFieldMetadata& metadata,
                                               s32 alt_idx);

template void UpdateFormatFieldForAltAllele<s32>(const io::VcfRecordPtr& record,
                                                 const io::FormatFieldMetadata& metadata,
                                                 s32 alt_idx);
template void UpdateFormatFieldForAltAllele<f32>(const io::VcfRecordPtr& record,
                                                 const io::FormatFieldMetadata& metadata,
                                                 s32 alt_idx);

vec<io::VcfRecordPtr> SplitMultiAllelicRecord(const io::VcfRecordPtr& record,
                                              const io::VcfHeaderPtr& new_header,
                                              const vec<io::InfoFieldMetadata>& info_metadata,
                                              const vec<io::FormatFieldMetadata>& fmt_metadata,
                                              const bool trim_variant) {
  vec<io::VcfRecordPtr> result;

  // Iterate over each alternate allele and create a new VcfRecord for it based on the original record.
  // Update the alleles, INFO fields, and FORMAT fields in the new VcfRecord using the metadata from the VCF header.
  const std::string ref_allele = record->Allele(0);
  for (s32 alt_idx = 1; alt_idx < record->NumAlleles(); ++alt_idx) {
    const auto& new_record = record->Clone(new_header);

    // Set alleles: REF and the current ALT allele
    const std::string alt_allele = record->Allele(alt_idx);
    if (trim_variant) {
      const auto& [trimmed_ref, trimmed_alt] = TrimVariant(ref_allele, alt_allele);
      new_record->SetAlleles({trimmed_ref, trimmed_alt});
    } else {
      new_record->SetAlleles({ref_allele, alt_allele});
    }

    // Update INFO fields
    for (const auto& metadata : info_metadata) {
      if (!record->HasInfoFieldNoCheck(metadata.id) || metadata.number == io::kNumberOne ||
          metadata.number == io::kNumberZero || metadata.number == io::kNumberDot) {
        // if field is not present, nothing needs to be done
        // if Number={1, 0, .}, the fields were copied correctly when the record was cloned
        continue;
      }
      if (metadata.type == io::FieldType::kInteger) {
        UpdateInfoFieldForAltAllele<s32>(new_record, metadata, alt_idx);
      } else if (metadata.type == io::FieldType::kFloat) {
        UpdateInfoFieldForAltAllele<f32>(new_record, metadata, alt_idx);
      }
      // Flags do not have values to set
      // String and Character INFO fields are not typically per-allele, so they are not handled here
    }

    // Update FORMAT fields
    for (const auto& metadata : fmt_metadata) {
      if (!record->HasFormatFieldNoCheck(metadata.id) || metadata.number == io::kNumberOne ||
          metadata.number == io::kNumberZero || metadata.number == io::kNumberDot) {
        // if field is not present, nothing needs to be done
        // if Number={1, 0, .}, the fields were copied correctly when the record was cloned
        continue;
      }
      if (metadata.type == io::FieldType::kInteger) {
        UpdateFormatFieldForAltAllele<s32>(new_record, metadata, alt_idx);
      } else if (metadata.type == io::FieldType::kFloat) {
        UpdateFormatFieldForAltAllele<f32>(new_record, metadata, alt_idx);
      }
      // Flags do not have values to set
      // String and Character INFO fields are not typically per-allele, so they are not handled here
    }

    result.push_back(new_record);
  }

  return result;
}

/**
 * @brief Merge a multi-allelic field (INFO or FORMAT) from two records into the first record.
 * @tparam RecordPtr Pointer type to VCF record
 * @tparam Metadata Metadata type (InfoFieldMetadata or FormatFieldMetadata)
 * @tparam GetFieldFn Function type to get the field values
 * @tparam SetFieldFn Function type to set the field values
 * @param record1 First VCF record (to be updated)
 * @param record2 Second VCF record
 * @param metadata Field metadata
 * @param get_field Function to get the field values
 * @param set_field Function to set the field values
 * @return true if the field was successfully merged, false otherwise
 */
template <typename RecordPtr, typename Metadata, typename GetFieldFn, typename SetFieldFn>
static void MergeMultiAllelicField(const RecordPtr& record1,
                                   const RecordPtr& record2,
                                   const Metadata& metadata,
                                   GetFieldFn get_field,
                                   SetFieldFn set_field) {
  const auto& field = metadata.id;
  const auto& values1 = get_field(record1, field);
  const auto& values2 = get_field(record2, field);
  if (metadata.number == io::kNumberR) {
    // For Number=R, both records should have 2 values: [REF, ALT]
    static const size_t kNumAlleles = 2;
    if (values1.size() >= kNumAlleles && values2.size() >= kNumAlleles) {
      // Merge as [REF, ALT1, ALT2]
      set_field(record1, field, {values1.at(0), values1.at(1), values2.at(1)});
    }
  } else if (metadata.number == io::kNumberEachAllele && !values1.empty() && !values2.empty()) {
    // For Number=A, both records should have one value
    // Merge as [ALT1, ALT2]
    set_field(record1, field, {values1.at(0), values2.at(0)});
  }
}

/**
 * @brief Merge a multi-allelic INFO field from two records into the first record.
 * @tparam T Field value type
 * @param record1 First VCF record (to be updated)
 * @param record2 Second VCF record
 * @param metadata INFO field metadata
 */
template <typename T>
static void MergeMultiAllelicInfoField(const io::VcfRecordPtr& record1,
                                       const io::VcfRecordPtr& record2,
                                       const io::InfoFieldMetadata& metadata) {
  MergeMultiAllelicField(
      record1,
      record2,
      metadata,
      [](const io::VcfRecordPtr& rec, const std::string& id) { return rec->GetInfoFieldNoCheck<T>(id); },
      [](const io::VcfRecordPtr& rec, const std::string& id, const vec<T>& vals) { rec->SetInfoField<T>(id, vals); });
}

/**
 * @brief Merge a multi-allelic FORMAT field from two records into the first record.
 * @tparam T Field value type
 * @param record1 First VCF record (to be updated)
 * @param record2 Second VCF record
 * @param metadata FORMAT field metadata
 */
template <typename T>
static void MergeMultiAllelicFormatField(const io::VcfRecordPtr& record1,
                                         const io::VcfRecordPtr& record2,
                                         const io::FormatFieldMetadata& metadata) {
  MergeMultiAllelicField(
      record1,
      record2,
      metadata,
      [](const io::VcfRecordPtr& rec, const std::string& id) { return rec->GetFormatFieldNoCheck<T>(id); },
      [](const io::VcfRecordPtr& rec, const std::string& id, const vec<T>& vals) { rec->SetFormatField<T>(id, vals); });
}

template void MergeMultiAllelicInfoField<s32>(const io::VcfRecordPtr& record1,
                                              const io::VcfRecordPtr& record2,
                                              const io::InfoFieldMetadata& metadata);
template void MergeMultiAllelicInfoField<f32>(const io::VcfRecordPtr& record1,
                                              const io::VcfRecordPtr& record2,
                                              const io::InfoFieldMetadata& metadata);

template void MergeMultiAllelicFormatField<s32>(const io::VcfRecordPtr& record1,
                                                const io::VcfRecordPtr& record2,
                                                const io::FormatFieldMetadata& metadata);
template void MergeMultiAllelicFormatField<f32>(const io::VcfRecordPtr& record1,
                                                const io::VcfRecordPtr& record2,
                                                const io::FormatFieldMetadata& metadata);

bool MergeMultiAllelicRecords(const io::VcfRecordPtr& record1,
                              const io::VcfRecordPtr& record2,
                              const vec<io::InfoFieldMetadata>& info_metadata,
                              const vec<io::FormatFieldMetadata>& fmt_metadata) {
  // Merge allele representations for the multi-allelic record
  const auto& [ref1, alt1] = TrimVariant(record1->Allele(0), record1->Allele(1));
  const auto& [ref2, alt2] = TrimVariant(record2->Allele(0), record2->Allele(1));
  const auto& [ref_new, alt1_new, alt2_new] = FormatVariants(ref1, alt1, ref2, alt2);
  if (ref_new.empty() || alt1_new.empty() || alt2_new.empty()) {
    return false;
  }
  record1->SetAlleles({ref_new, alt1_new, alt2_new});

  // Merge INFO fields as needed
  for (const auto& metadata : info_metadata) {
    if (!record1->HasInfoFieldNoCheck(metadata.id) || !record2->HasInfoFieldNoCheck(metadata.id)) {
      // if field is not present in both records, nothing needs to be done
      continue;
    }
    if (metadata.number == io::kNumberOne || metadata.number == io::kNumberZero || metadata.number == io::kNumberDot) {
      // if Number={1, 0, .}, the fields were copied correctly when the record was cloned
      continue;
    }
    if (metadata.type == io::FieldType::kInteger) {
      MergeMultiAllelicInfoField<s32>(record1, record2, metadata);
    } else if (metadata.type == io::FieldType::kFloat) {
      MergeMultiAllelicInfoField<f32>(record1, record2, metadata);
    }
    // Flags do not have values to set
    // String and Character INFO fields are not typically per-allele, so they are not handled here
  }

  // Merge FORMAT fields as needed
  for (const auto& metadata : fmt_metadata) {
    if (!record1->HasFormatFieldNoCheck(metadata.id) || !record2->HasFormatFieldNoCheck(metadata.id)) {
      // if field is not present in both records, nothing needs to be done
      continue;
    }
    if (metadata.number == io::kNumberOne || metadata.number == io::kNumberZero || metadata.number == io::kNumberDot) {
      // if Number={1, 0, .}, the fields were copied correctly when the record was cloned
      continue;
    }
    if (metadata.type == io::FieldType::kInteger) {
      MergeMultiAllelicFormatField<s32>(record1, record2, metadata);
    } else if (metadata.type == io::FieldType::kFloat) {
      MergeMultiAllelicFormatField<f32>(record1, record2, metadata);
    }
    // Flags do not have values to set
    // String and Character FORMAT fields are not typically per-allele, so they are not handled here
  }

  // Update AC, AF, AN, and GT field values for multi-allelic record
  static const s32 kMultiallelicAc = 1;
  static const f32 kMultiallelicAf = 0.5;
  static const s32 kMultiallelicAn = 2;
  record1->SetInfoField<s32>(kFieldAc, {kMultiallelicAc, kMultiallelicAc});
  record1->SetInfoField<f32>(kFieldAf, {kMultiallelicAf, kMultiallelicAf});
  record1->SetInfoField<s32>(kFieldAn, {kMultiallelicAn});
  record1->SetGTField(GenotypeToString(Genotype::kGT12));

  return true;
}

}  // namespace xoos::svc
