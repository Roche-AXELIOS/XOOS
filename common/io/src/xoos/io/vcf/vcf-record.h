#pragma once

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <htslib/vcf.h>

#include <xoos/io/htslib-util/kstring.h>

#include "vcf-header.h"

namespace xoos::io {
constexpr auto kGT = "GT";

inline bool IsGtStandardFormatField(const std::string& field) {
  return field == kGT;
}

class VcfRecord;
using VcfRecordPtr = std::shared_ptr<VcfRecord>;
using BcfRecordPtr = std::shared_ptr<bcf1_t>;

class VcfRecord {
 public:
  explicit VcfRecord(const VcfHeaderPtr& hdr);
  static VcfRecordPtr CreateFromHeader(const VcfHeaderPtr& hdr, int unpack);
  static VcfRecordPtr CreateFromHeader(const VcfHeaderPtr& hdr);
  static VcfRecordPtr ReadFromFile(const VcfHeaderPtr& hdr, const HtsFileSharedPtr& input_vcf_name, int unpack);
  static VcfRecordPtr ReadFromFile(const VcfHeaderPtr& hdr, const HtsFileSharedPtr& input_vcf_name);
  static VcfRecordPtr ReadFromRegion(const VcfHeaderPtr& hdr,
                                     const HtsFileSharedPtr& input_vcf_fp,
                                     tbx_t* tbx_idx,
                                     hts_itr_t* hts_itr,
                                     Kstring& kstr,
                                     int unpack);

  std::string Chromosome() const;
  hts_pos_t Position() const;
  bool IsSnp() const;
  std::string Id() const;
  std::string Allele(int which) const;
  int NumAlleles() const;

  bool HasInfoFieldNoCheck(const std::string& field) const;
  template <typename T>
  std::vector<T> GetInfoFieldNoCheck(const std::string& field) const;
  std::string GetInfoFieldStringNoCheck(const std::string& field) const;
  bool HasFormatFieldNoCheck(const std::string& field) const;
  template <typename T>
  std::vector<T> GetFormatFieldNoCheck(const std::string& field) const;
  std::string GetFormatFieldStringNoCheck(const std::string& field) const;

  template <typename T>
  std::vector<T> GetInfoField(const std::string& field) const;
  template <typename T>
  void SetInfoField(const std::string& field, const std::vector<T>& values);

  template <typename T>
  std::vector<T> GetFormatField(const std::string& field) const;
  template <typename T>
  void SetFormatField(const std::string& field, const std::vector<T>& values);
  std::string GetGTField() const;
  void SetGTField(const std::string& value);

  void SetChromosome(const std::string& chromosome);
  void SetPosition(int position);
  void SetId(const std::string& id);
  void SetAlleles(const std::vector<std::string>& alleles);
  std::optional<float> GetQuality() const;
  float GetQuality(float default_value) const;
  void SetQuality(const std::optional<float>& quality);
  void SetFilter(const std::string& filter);
  void AddFilter(const std::string& filter);
  void ClearFilters();
  std::vector<std::string> GetFilters() const;
  VcfRecordPtr Clone(const VcfHeaderPtr& new_header) const;

  // used to expose _record to VcfWriter
  friend class VcfWriter;

 private:
  // variables
  BcfHeaderPtr _hdr;
  BcfRecordPtr _record;
  static constexpr auto kGenotypeRegex = R"(^(\.|\d+)(([/|])(\.|\d+))?$)";
};
}  // namespace xoos::io
