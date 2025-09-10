#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "vcf-header.h"
#include "vcf-record.h"

namespace xoos::io {

namespace fs = std::filesystem;

// for unknown reference or sample name (or other string fields)
constexpr auto kUnspecified = "unspecified";

using VcfIntegerFields = std::unordered_map<std::string, std::vector<int>>;
using VcfFloatFields = std::unordered_map<std::string, std::vector<float>>;
using VcfStringFields = std::unordered_map<std::string, std::vector<std::string>>;

struct TypedVcfFields {
  const VcfIntegerFields& integer_fields;
  const VcfFloatFields& float_fields;
  const VcfStringFields& string_fields;
  std::optional<std::vector<std::string>> field_order = std::nullopt;
};

/**
 * @brief Using https://github.com/EBIvariation/vcf-validator as a reference for writing VCF files
 */
class VcfWriter {
 public:
  explicit VcfWriter(const fs::path& vcf_file_path, VcfHeaderPtr header = VcfHeaderPtr());
  void WriteHeader() const;
  void WriteHeader(const std::vector<std::string>& custom_meta_data_lines,
                   const std::vector<FilterFieldMetadata>& filter_lines,
                   const std::vector<InfoFieldMetadata>& info_lines,
                   const std::vector<FormatFieldMetadata>& format_lines,
                   const std::vector<ContigMetadata>& contig_lines,
                   const std::string& sample_name) const;
  VcfRecordPtr CreateRecord(const std::string& chromosome,
                            int position,
                            const std::string& id,
                            const std::vector<std::string>& alleles,
                            const std::optional<float>& quality,
                            const std::string& filter_name,
                            const TypedVcfFields& info_fields,
                            const TypedVcfFields& format_fields);
  void WriteRecord(const VcfRecordPtr& record) const;
  void Flush();

 private:
  HtsFileSharedPtr _file;
  VcfHeaderPtr _hdr;
  std::string _file_path;
};

}  // namespace xoos::io
