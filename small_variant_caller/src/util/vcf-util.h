#pragma once

#include <functional>

#include <xoos/types/fs.h>
#include <xoos/types/str-container.h>

#include "core/variant-info.h"
#include "xoos/io/vcf/vcf-record.h"

namespace xoos::svc {

/**
 * @brief Snapshot of the per-allele (Number=R or Number=A) FORMAT field values of a bi-allelic VCF record.
 *
 * @details Used to carry a spanning (upstream) deletion's computed per-allele counts (e.g. `AD`, `ADC`, `ADS`, `ADD`,
 * `ADL`) forward to a downstream position where the deletion appears as a wildcard (`*`) allele. The snapshot copies
 * the field value vectors so it remains valid even if the source record is later mutated (e.g. merged into a
 * multi-allelic record).
 */
struct PerAlleleFieldSnapshot {
  StrMap<vec<s32>> int_fields;
  StrMap<vec<f32>> float_fields;
};

/**
 * @brief Non-owning view of a VCF header's INFO and FORMAT field metadata, bundled together.
 * @details These two vectors are always sourced and passed together (e.g. from
 * `io::VcfHeader::GetFieldMetadata()`). Bundling them keeps function signatures compact. The members are
 * `std::reference_wrapper`s rather than owning vectors so that constructing this view does not copy the
 * metadata vectors on every call. The referenced vectors must outlive every use of the view.
 */
struct VcfFieldMetadata {
  std::reference_wrapper<const vec<io::InfoFieldMetadata>> info;
  std::reference_wrapper<const vec<io::FormatFieldMetadata>> format;
};

/**
 * @brief Check if a VCF record is multi-allelic, i.e. 2 or more ALT alleles.
 * @param record VCF record
 * @return True if the record is multi-allelic, false otherwise
 */
bool IsMultiAllelicRecord(const io::VcfRecordPtr& record);

/**
 * @brief Copy the per-allele (Number=R and Number=A) FORMAT field values of a bi-allelic record into a snapshot.
 * @param record Bi-allelic VCF record to snapshot (e.g. a spanning deletion record)
 * @param fmt_metadata Metadata for all FORMAT fields declared in the VCF header
 * @return Snapshot holding a copy of each present per-allele integer and float FORMAT field
 */
PerAlleleFieldSnapshot SnapshotPerAlleleFormatFields(const io::VcfRecordPtr& record,
                                                     const vec<io::FormatFieldMetadata>& fmt_metadata);

/**
 * @brief Overwrite a record's per-allele FORMAT fields with the values held in a snapshot.
 *
 * @details Applied to a wildcard (`*`) split record so that, after merging, the `*` allele slot carries the spanning
 * deletion's computed counts instead of stale or missing values. Only fields present in the snapshot are set.
 *
 * @param record VCF record to update (typically a bi-allelic `[REF, *]` record)
 * @param snapshot Per-allele field values copied from the spanning deletion record
 */
void ApplyPerAlleleFormatFields(const io::VcfRecordPtr& record, const PerAlleleFieldSnapshot& snapshot);

/**
 * @brief Set the non-model per-sample FORMAT data fields on a VCF record from BAM features.
 *
 * @details Populates the allele depth (`AD`), total depth (`DP`), and, when the sequencing protocol is duplex, the
 * duplex count fields (concordant, simplex, discordant, low-base-quality). These are the data fields shared by the
 * germline and germline-mitochondria workflows; model-related FORMAT fields are intentionally not set here.
 *
 * @param record VCF record to be updated
 * @param ref_feat Reference-allele BAM features
 * @param var_feat Variant-allele BAM features
 * @param is_duplex_protocol Whether the sequencing protocol is duplex
 */
void SetRecordDataFormatFields(const io::VcfRecordPtr& record,
                               const UnifiedReferenceFeature& ref_feat,
                               const UnifiedVariantFeature& var_feat,
                               bool is_duplex_protocol);

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

/**
 * @brief Whether a VCF file holds at least one record.
 *
 * @details Opens its own reader and reads a single record, so it does not disturb a reader the caller
 * already holds. Call this before planning work from the VCF index: a valid VCF with no records indexes to
 * an index that covers no sequence, and `VcfReader::GetContigIndexes` rejects that as invalid because it
 * cannot tell it apart from an index that does not belong to the file.
 *
 * @param vcf_file Path to the VCF file.
 * @return true if the file holds at least one record.
 */
bool VcfHasRecords(const fs::path& vcf_file);
}  // namespace xoos::svc
