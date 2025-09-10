#pragma once

#include <gsl/gsl>
#include <map>
#include <memory>
#include <string>
#include <utility>

#include <htslib/tbx.h>
#include <htslib/vcf.h>

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/int.h>

namespace xoos::io {

using HtsFileSharedPtr = std::shared_ptr<htsFile>;

enum class FieldType {
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
  uint32_t length;
};

class VcfHeader;
using VcfHeaderPtr = std::shared_ptr<VcfHeader>;
using BcfHeaderPtr = std::shared_ptr<bcf_hdr_t>;

class VcfHeader {
 public:
  explicit VcfHeader(BcfHeaderPtr hdr) : _hdr(std::move(hdr)) {
  }

  static VcfHeaderPtr Read(const HtsFileSharedPtr& input_file);
  static VcfHeaderPtr Create();
  VcfHeaderPtr Clone() const;
  void Write(const HtsFileSharedPtr& output_file);
  void AddCustomMetaDataLine(const std::string& line);
  void AddInfoLine(const InfoFieldMetadata& info_line);
  void AddFilterLine(const FilterFieldMetadata& info_filter);
  void AddFormatLine(const FormatFieldMetadata& format_line);
  void AddContigLine(const ContigMetadata& line);
  void AddSample(const std::string& sample_name);
  int GetContigId(const std::string& ctg);
  std::map<std::string, s32> GetContigIndexes();
  std::map<std::string, u64> GetContigLengths();
  char* GetValue(const char* key);
  std::map<std::string, int> GetSampleIndexes();
  bool HasInfoField(const std::string& name);

  /**
   * @brief Synchronize the header data structures after modifications.
   *
   * This method must be called after making structural changes to the header
   * (such as adding samples, INFO fields, FORMAT fields, etc.) to ensure
   * internal data structures are properly updated and consistent.
   *
   * @throws std::runtime_error if synchronization fails
   */
  void Sync();

  int GetNumSamples() const {
    return bcf_hdr_nsamples(_hdr);
  }

  // using to give access to `Get()`
  friend class VcfRecord;
  friend class VcfWriter;

 private:
  BcfHeaderPtr _hdr;
};
}  // namespace xoos::io
