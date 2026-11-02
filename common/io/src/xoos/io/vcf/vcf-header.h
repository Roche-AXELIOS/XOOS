#pragma once

#include <gsl/gsl>
#include <map>
#include <memory>
#include <mutex>
#include <string>

#include <htslib/tbx.h>
#include <htslib/vcf.h>

#include <xoos/io/htslib-util/htslib-bcf-ptr.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>

namespace xoos::io {

using HtsFileSharedPtr = std::shared_ptr<htsFile>;

enum class FieldType : u8 {
  kInteger,
  kFloat,
  kCharacter,
  kString,
  kFlag
};

// Number of alleles for each alternate allele, excluding REF
const std::string kNumberEachAllele = "A";
// Number of alleles for each allele, including REF
const std::string kNumberR = "R";
const std::string kNumberOne = "1";
const std::string kNumberZero = "0";
const std::string kNumberDot = ".";

std::string FieldTypeToString(FieldType type);

struct InfoFieldMetadata {
  std::string id;
  std::string description;
  std::string number;
  FieldType type;
};

struct FormatFieldMetadata {
  std::string id;
  std::string description;
  std::string number;
  FieldType type;
};

struct FilterFieldMetadata {
  std::string id;
  std::string description;
};

struct ContigMetadata {
  std::string id;
  u32 length;
};

class VcfHeader;
using VcfHeaderPtr = std::shared_ptr<VcfHeader>;
using BcfHeaderPtr = std::shared_ptr<bcf_hdr_t>;

struct VcfHeaderState;

struct VcfHeaderTranslationContext {
  /**
   * @brief Snapshot a source header for one source/destination layout pair.
   * @param source_header Source header to duplicate for htslib translation.
   * @param source_revision Revision of the source header snapshot.
   * @param destination_revision Revision of the destination header snapshot.
   */
  VcfHeaderTranslationContext(const bcf_hdr_t* source_header, u64 source_revision, u64 destination_revision);

  BcfHdrPtr source_header_snapshot;
  /// Revisions prevent reuse after either aliased header changes its layout.
  u64 source_revision_snapshot;
  u64 destination_revision_snapshot;
  /// htslib updates translation internals during bcf_translate; serialize that access.
  std::mutex mutex;
};

using VcfHeaderTranslationContextPtr = std::shared_ptr<VcfHeaderTranslationContext>;

class VcfHeader {
 public:
  /**
   * @brief Wrap an htslib VCF header.
   * @param hdr Shared header pointer to wrap. Wrappers sharing this pointer also share
   *        revision tracking and synchronization state.
   */
  explicit VcfHeader(BcfHeaderPtr hdr);

  static VcfHeaderPtr Read(const HtsFileSharedPtr& input_file);
  static VcfHeaderPtr Create();
  VcfHeaderPtr Clone() const;
  void Write(const HtsFileSharedPtr& output_file) const;
  // S5817: these methods mutate the underlying bcf_hdr_t via bcf_hdr_append/bcf_hdr_add_sample;
  // the shared_ptr member is unchanged but the pointed-to htslib header is modified.
  void AddCustomMetaDataLine(const std::string& line);         // NOSONAR(S5817)
  void AddInfoLine(const InfoFieldMetadata& info_line);        // NOSONAR(S5817)
  void AddFilterLine(const FilterFieldMetadata& info_filter);  // NOSONAR(S5817)
  void AddFormatLine(const FormatFieldMetadata& format_line);  // NOSONAR(S5817)
  void AddContigLine(const ContigMetadata& line);              // NOSONAR(S5817)

  /**
   * @brief Remove one INFO meta-information line from the header by field ID.
   *
   * Drops only the matching `##INFO=<ID=...>` definition while leaving all
   * other INFO lines and other header lines intact. Removing an absent field is
   * a no-op. Call `Sync()` after remove/add edits before reading metadata.
   */
  void RemoveInfoLine(const std::string& id) const;

  /**
   * @brief Remove all INFO meta-information lines from the header.
   *
   * Drops every `##INFO` definition while leaving all other meta lines
   * (`##FORMAT`, `##contig`, reference, `##FILTER`, custom lines) and samples
   * intact. Call `Sync()` after re-adding any desired INFO lines so the header's
   * internal dictionary is rebuilt.
   */
  void RemoveAllInfoLines() const;

  /**
   * @brief Remove all FILTER meta-information lines from the header.
   *
   * Drops every `##FILTER` definition while leaving all other meta lines
   * (`##INFO`, `##FORMAT`, `##contig`, reference, custom lines) and samples
   * intact. Call `Sync()` after re-adding any desired FILTER lines so the
   * header's internal dictionary is rebuilt.
   */
  void RemoveAllFilterLines() const;

  /**
   * @brief Remove all FORMAT meta-information lines from the header.
   *
   * Drops every `##FORMAT` definition while leaving all other meta lines
   * (`##INFO`, `##contig`, reference, `##FILTER`, custom lines) and samples
   * intact. Call `Sync()` after re-adding any desired FORMAT lines so the
   * header's internal dictionary is rebuilt.
   */
  void RemoveAllFormatLines() const;

  void AddSample(const std::string& sample_name);  // NOSONAR(S5817)
  s32 GetContigId(const std::string& ctg) const;
  StrMap<s32> GetContigIndexes() const;
  StrMap<u64> GetContigLengths() const;
  char* GetValue(const char* key) const;
  StrMap<s32> GetSampleIndexes() const;
  bool HasInfoField(const std::string& name) const;
  bool HasFilterField(const std::string& name) const;

  /**
   * @brief Check whether a FORMAT field is defined in the VCF header.
   * @param name FORMAT field ID to look up.
   * @return True when the header contains a `##FORMAT` definition with this ID.
   */
  bool HasFormatField(const std::string& name) const;

  /**
   * @brief Retrieve metadata for INFO and FORMAT fields in the VCF header.
   *
   * This method extracts and returns the metadata for all INFO and FORMAT
   * fields defined in the VCF header. The metadata includes field ID,
   * description, number, and type.
   *
   * @return A tuple containing two vectors:
   *         - First vector contains InfoFieldMetadata for INFO fields.
   *         - Second vector contains FormatFieldMetadata for FORMAT fields.
   */
  std::tuple<vec<InfoFieldMetadata>, vec<FormatFieldMetadata>> GetFieldMetadata() const;
  vec<FilterFieldMetadata> GetFilterMetadata() const;

  /**
   * @brief Translate a record into a destination header using a cached source-header context.
   *
   * The source header is duplicated once per destination header and reused for subsequent records. The
   * translation itself is serialized per context because htslib updates its translation cache on the
   * source header.
   * @param destination Header whose dictionary should receive the translated record.
   * @param record Record to translate in place.
   * @throws error::Error if htslib cannot translate the record.
   * @thread_safety Safe to call concurrently with other translations and synchronized header mutations.
   */
  void TranslateRecord(const VcfHeaderPtr& destination, bcf1_t* record) const;

  /**
   * @brief Synchronize the header data structures after modifications.
   *
   * This method must be called after making structural changes to the header
   * (such as adding samples, INFO fields, FORMAT fields, etc.) to ensure
   * internal data structures are properly updated and consistent.
   *
   * @throws error::Error if synchronization fails
   */
  void Sync() const;

  s32 GetNumSamples() const {
    return bcf_hdr_nsamples(_hdr);
  }

  // using to give access to `Get()`
  friend class VcfRecord;
  friend class VcfWriter;

 private:
  VcfHeaderTranslationContextPtr GetTranslationContext(const VcfHeaderPtr& destination,
                                                       u64 source_revision,
                                                       u64 destination_revision) const;
  void InvalidateTranslationContexts() const;

  BcfHeaderPtr _hdr;
  /// Shared with other VcfHeader wrappers that alias the same underlying bcf_hdr_t.
  std::shared_ptr<VcfHeaderState> _state;
  mutable std::mutex _translation_contexts_mutex;
  mutable std::map<std::weak_ptr<VcfHeader>, VcfHeaderTranslationContextPtr, std::owner_less<std::weak_ptr<VcfHeader>>>
      _translation_contexts;
};
}  // namespace xoos::io
