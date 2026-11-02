#pragma once

#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include <htslib/vcf.h>

#include <xoos/io/htslib-util/kstring.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>

#include "vcf-header.h"

namespace xoos::io {
const constexpr char* const kGT = "GT";

inline bool IsGtStandardFormatField(const std::string_view field) {
  return field == kGT;
}

class VcfRecord;
using VcfRecordPtr = std::shared_ptr<VcfRecord>;
using BcfRecordPtr = std::shared_ptr<bcf1_t>;

// S1448: VcfRecord is the primary VCF record API; splitting the class would fragment the interface.
class VcfRecord {  // NOSONAR(S1448)
 public:
  explicit VcfRecord(const VcfHeaderPtr& hdr);
  static VcfRecordPtr CreateFromHeader(const VcfHeaderPtr& hdr, s32 unpack);
  static VcfRecordPtr CreateFromHeader(const VcfHeaderPtr& hdr);
  static VcfRecordPtr ReadFromFile(const VcfHeaderPtr& hdr, const HtsFileSharedPtr& input_vcf_name, s32 unpack);
  static VcfRecordPtr ReadFromFile(const VcfHeaderPtr& hdr, const HtsFileSharedPtr& input_vcf_name);
  static VcfRecordPtr ReadFromRegion(const VcfHeaderPtr& hdr,
                                     const HtsFileSharedPtr& input_vcf_fp,
                                     tbx_t* tbx_idx,
                                     hts_itr_t* hts_itr,
                                     Kstring& kstr,
                                     s32 unpack);

  std::string Chromosome() const;

  /// Returns the 0-based position from htslib's internal representation (bcf1_t::pos).
  /// This is the same coordinate system as BED start positions.
  /// Add 1 when outputting to VCF format (which uses 1-based positions).
  hts_pos_t Position() const;

  bool IsSnp() const;
  bool IsPass() const;
  std::string Id() const;
  std::string Allele(s32 which) const;
  s32 NumAlleles() const;

  bool HasInfoFieldNoCheck(const std::string& field) const;
  template <typename T>
  std::vector<T> GetInfoFieldNoCheck(const std::string& field) const;
  std::string GetInfoFieldStringNoCheck(const std::string& field) const;
  bool HasFormatFieldNoCheck(const std::string& field) const;
  template <typename T>
  std::vector<T> GetFormatFieldNoCheck(const std::string& field) const;
  std::string GetFormatFieldStringNoCheck(const std::string& field) const;

  void ConfirmFieldInHeader(const std::string& field, const std::string& field_type) const;
  void ConfirmInfoFieldInHeader(const std::string& field) const;
  void ConfirmFormatFieldInHeader(const std::string& field) const;

  template <typename T>
  std::vector<T> GetInfoField(const std::string& field) const;
  // S5817: these methods mutate the underlying bcf1_t data through a raw pointer; the pointer
  // member is unchanged but the pointed-to htslib record is modified, so const is wrong here.
  void SetInfoField(const std::string& field, const std::vector<f32>& values);          // NOSONAR(S5817)
  void SetInfoField(const std::string& field, const std::vector<s32>& values);          // NOSONAR(S5817)
  void SetInfoField(const std::string& field, const std::vector<std::string>& values);  // NOSONAR(S5817)
  void AddInfoFieldFlag(const std::string& field);                                      // NOSONAR(S5817)
  void RemoveInfoFieldFlag(const std::string& field);                                   // NOSONAR(S5817)

  template <typename T>
  std::vector<T> GetFormatField(const std::string& field) const;
  void SetFormatField(const std::string& field, const std::vector<f32>& values);  // NOSONAR(S5817)
  void SetFormatField(const std::string& field, const std::vector<s32>& values);  // NOSONAR(S5817)
  void SetFormatField(const std::string& field, const std::vector<const char*>& values);
  void SetFormatField(const std::string& field, const std::vector<std::string>& values);
  std::string GetGTField() const;
  std::string GetGTField(s32 sample_index) const;
  void SetGTField(const std::string& value);
  void SetGTField(const std::string& value, s32 sample_index);

  void SetChromosome(const std::string& chromosome) const;
  /// Set the 0-based position (htslib internal). Pass a 0-based coordinate, not a VCF 1-based POS.
  void SetPosition(s32 position) const;
  void SetId(const std::string& id) const;
  void SetAlleles(const std::vector<std::string>& alleles) const;
  std::optional<f32> GetQuality() const;
  f32 GetQuality(f32 default_value) const;
  void SetQuality(const std::optional<f32>& quality) const;
  void SetFilter(const char* filter) const;
  void SetFilter(const std::string& filter) const;
  void SetFilter(std::string_view filter) const;
  void AddFilter(const char* filter) const;
  void AddFilter(const std::string& filter) const;
  void AddFilter(std::string_view filter) const;
  void ClearFilters() const;
  std::vector<std::string> GetFilters() const;
  VcfRecordPtr Clone(const VcfHeaderPtr& new_header) const;

  // used to expose _record to VcfWriter
  friend class VcfWriter;

 private:
  // Parse a VCF GT string (e.g. "0/1" or "1|0") into a vector of htslib GT integers.
  // Throws if the string does not match the expected genotype format.
  static std::vector<s32> ParseGenotypeString(const std::string& value);
  // Build the full multi-sample GT buffer for bcf_update_genotypes, writing the given alleles for
  // the target sample while preserving the other samples' existing genotypes.
  std::vector<s32> BuildGenotypeBuffer(const std::vector<s32>& alleles, size_t sample_idx) const;

  // variables
  VcfHeaderPtr _header;
  BcfHeaderPtr _hdr;
  BcfRecordPtr _record;
  static constexpr char kPhasedAlleleSeparator = '|';
  static constexpr char kUnphasedAlleleSeparator = '/';
  static constexpr char kMissingAllele = '.';
  // Matches a VCF GT string: a single allele (`.` or integer), optionally followed by more
  // alleles joined by `/` (unphased) or `|` (phased), but not a mix of both separators.
  static const constexpr char* const kGenotypeRegex = R"(^(?:\.|\d+)(?:(?:/(?:\.|\d+))*|(?:\|(?:\.|\d+))*)$)";
};
}  // namespace xoos::io
