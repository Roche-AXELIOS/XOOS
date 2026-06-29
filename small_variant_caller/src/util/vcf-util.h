#pragma once

#include "xoos/io/vcf/vcf-record.h"

namespace xoos::svc {

/**
 * @brief Check if a VCF record is multi-allelic, i.e. 2 or more ALT alleles.
 * @param record VCF record
 * @return True if the record is multi-allelic, false otherwise
 */
bool IsMultiAllelicRecord(const io::VcfRecordPtr& record);

/**
 * @brief Split a multi-allelic VCF record into multiple bi-allelic records. If any INFO or FORMAT fields have
 * unexpected number of values according to the VCF header metadata, these fields are not split.
 * @param record Original multi-allelic VCF record
 * @param new_header VCF header for the new VCF records
 * @param info_metadata Metadata for INFO fields from the VCF header for the original record
 * @param fmt_metadata Metadata for FORMAT fields from the VCF header for the original record
 * @param trim_variant Whether to trim the variants when splitting
 * @return Vector of bi-allelic VCF records
 */
vec<io::VcfRecordPtr> SplitMultiAllelicRecord(const io::VcfRecordPtr& record,
                                              const io::VcfHeaderPtr& new_header,
                                              const vec<io::InfoFieldMetadata>& info_metadata,
                                              const vec<io::FormatFieldMetadata>& fmt_metadata,
                                              bool trim_variant);

/**
 * @brief Creates a multi-allelic variant record by combining two single-allelic VCF records.
 *
 * @details This function takes two VCF records that represent different alleles at the same genomic position
 * and combines them into a single multi-allelic record based on the VCF header metadata for INFO and FORMAT fields.
 * The records must share the same chromosome and position but have different alternative alleles.
 * Each INFO or FORMAT field from the two VCF records is merged only if both records contain the same field and the
 * field type is either "R" or "A" in the VCF header metadata.
 *
 * @param record1 Pointer to the first VCF record containing one allele
 * @param record2 Pointer to the second VCF record containing another allele at the same position
 * @param info_metadata Metadata for INFO fields from the VCF header
 * @param fmt_metadata Metadata for FORMAT fields from the VCF header
 *
 * @return true if the multi-allelic record was successfully created, false otherwise
 *
 * @note The function assumes both input records are valid and non-null
 */
bool MergeMultiAllelicRecords(const io::VcfRecordPtr& record1,
                              const io::VcfRecordPtr& record2,
                              const vec<io::InfoFieldMetadata>& info_metadata,
                              const vec<io::FormatFieldMetadata>& fmt_metadata);
}  // namespace xoos::svc
