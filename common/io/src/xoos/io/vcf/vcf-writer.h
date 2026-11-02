#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include <htslib/bgzf.h>

#include "vcf-header.h"
#include "vcf-record.h"

namespace xoos::io {

namespace fs = std::filesystem;

// for unknown reference or sample name (or other string fields)
const constexpr char* const kUnspecified = "unspecified";

// S6045: VcfXxxFields are public type aliases used across the codebase; adding transparent comparators
// is a breaking API change that requires all downstream callers to update their local map types.
using VcfIntegerFields = std::unordered_map<std::string, std::vector<int>>;         // NOSONAR(S6045)
using VcfFloatFields = std::unordered_map<std::string, std::vector<float>>;         // NOSONAR(S6045)
using VcfStringFields = std::unordered_map<std::string, std::vector<std::string>>;  // NOSONAR(S6045)

// Transient parameter bundle: constructed at a call site and consumed immediately
// by WriteRecord/ApplyTypedFields, never stored or copied. The const-reference
// members deliberately alias the caller's field maps to avoid copying them.
// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
struct TypedVcfFields {
  const VcfIntegerFields& integer_fields;
  const VcfFloatFields& float_fields;
  const VcfStringFields& string_fields;
  std::optional<std::vector<std::string>> field_order = std::nullopt;
};

// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

/// Tag type to construct a VcfWriter that does not build an index on close.
struct NoIndex {};

/**
 * @brief Writes VCF/BCF files using htslib.
 *
 * Reference: https://github.com/EBIvariation/vcf-validator
 *
 * For .vcf.gz outputs the default constructor builds a .csi index when the file is closed.
 * Use the NoIndex overload for intermediate files that will be concatenated or do not need an index.
 */
class VcfWriter {
 public:
  /**
   * @brief Construct a VcfWriter that builds a .csi index on close for .vcf.gz outputs.
   * @param vcf_file_path Output file path. Files ending in ".gz" are BGZF-compressed.
   */
  explicit VcfWriter(const fs::path& vcf_file_path);

  /**
   * @brief Construct a VcfWriter that builds a .csi index on close for .vcf.gz outputs.
   * @param vcf_file_path Output file path. Files ending in ".gz" are BGZF-compressed.
   * @param header VCF header to use. If nullptr, an empty header is created via VcfHeader::Create().
   */
  VcfWriter(const fs::path& vcf_file_path, VcfHeaderPtr header);

  /**
   * @brief Construct a VcfWriter without automatic index building.
   *
   * Use this overload when writing temporary or intermediate files (e.g., header-only files
   * for BGZF concatenation) where indexing is unnecessary or will be done separately.
   *
   * @param vcf_file_path Output file path. Files ending in ".gz" are BGZF-compressed.
   * @param header VCF header to use. If nullptr, an empty header is created via VcfHeader::Create().
   */
  VcfWriter(const fs::path& vcf_file_path, VcfHeaderPtr header, NoIndex no_index);

  /**
   * @brief Write the VCF header that was provided at construction time.
   *
   * Triggers bcf_hdr_sync which rebuilds the header's internal dictionary. Must be called
   * before any operations that read the dictionary (e.g., GetFieldMetadata).
   */
  void WriteHeader() const;

  /**
   * @brief Populate the header with the given metadata and write it.
   *
   * Adds all provided metadata lines, sample name, and field definitions to the header,
   * then writes it to the output file.
   *
   * @param custom_meta_data_lines Free-form VCF meta-information lines (e.g., "##source=...").
   * @param filter_lines FILTER field definitions.
   * @param info_lines INFO field definitions.
   * @param format_lines FORMAT field definitions.
   * @param contig_lines Contig definitions with name and length.
   * @param sample_name Sample name to add to the header.
   */
  void WriteHeader(const std::vector<std::string>& custom_meta_data_lines,
                   const std::vector<FilterFieldMetadata>& filter_lines,
                   const std::vector<InfoFieldMetadata>& info_lines,
                   const std::vector<FormatFieldMetadata>& format_lines,
                   const std::vector<ContigMetadata>& contig_lines,
                   const std::string& sample_name) const;

  /**
   * @brief Create a new VCF record populated with the given fields.
   *
   * The record is associated with this writer's header. Field values are set in the order
   * specified by TypedVcfFields::field_order if present, otherwise in unordered-map iteration order.
   *
   * @param chromosome Chromosome/contig name (e.g., "chr1").
   * @param position Genomic position (1-based, as in the VCF POS column).
   * @param id Variant identifier (e.g., rsID) or "." for missing.
   * @param alleles Allele strings; the first element is the REF allele.
   * @param quality Variant quality score (QUAL column). std::nullopt sets QUAL to missing (".").
   * @param filter_name Filter status string (e.g., "PASS").
   * @param info_fields INFO field values grouped by type (integer, float, string).
   * @param format_fields FORMAT/sample field values grouped by type (integer, float, string).
   * @return Shared pointer to the created VcfRecord.
   */
  VcfRecordPtr CreateRecord(const std::string& chromosome,
                            s32 position,
                            const std::string& id,
                            const std::vector<std::string>& alleles,
                            const std::optional<f32>& quality,
                            const std::string& filter_name,
                            const TypedVcfFields& info_fields,
                            const TypedVcfFields& format_fields);

  /**
   * @brief Write a VCF record to the output file using bcf_write.
   * @param record The record to write. Must have been created with a compatible header.
   * @throws error::Error if the write fails.
   */
  void WriteRecord(const VcfRecordPtr& record) const;

  /**
   * @brief Flush any buffered output to disk.
   */
  void Flush() const;

  /**
   * @brief Format a VCF record as text and write it directly to a BGZF handle.
   *
   * Uses vcf_format + bgzf_write instead of bcf_write. Does not require bcf_hdr_write
   * to have been called on the handle, making it safe for writing to headerless data files
   * used in BGZF block-level concatenation.
   *
   * @param bgzf_file Open BGZF handle to write to.
   * @param hdr VCF header used to format the record fields.
   * @param record The record to format and write.
   * @throws error::Error if formatting or writing fails.
   */
  static void WriteBgzfRecord(BGZF* bgzf_file, const VcfHeaderPtr& hdr, const VcfRecordPtr& record);

 private:
  HtsFileSharedPtr _file = nullptr;
  VcfHeaderPtr _hdr;
  std::string _file_path;
};

}  // namespace xoos::io
