#include "xoos/io/vcf/vcf-writer.h"

#include <string>
#include <utility>
#include <vector>

#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/vcf.h>

#include <xoos/error/error.h>
#include <xoos/io/htslib-util/htslib-bcf-util.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/log/logging.h>
#include <xoos/util/adopt.h>

namespace xoos::io {

namespace {

/**
 * @brief Apply the typed values in a TypedVcfFields group to a record via the given setter.
 *
 * When `field_order` is present, fields are applied in that order (the first matching typed map wins);
 * otherwise all integer, float, then string fields are applied in map iteration order.
 *
 * @param fields The typed INFO or FORMAT field values.
 * @param set Callable invoked as `set(field_name, value)` for each present field; overloaded on value type.
 */
template <typename SetField>
void ApplyTypedFields(const TypedVcfFields& fields, const SetField& set) {
  if (fields.field_order.has_value()) {
    for (const auto& field : fields.field_order.value()) {
      if (fields.integer_fields.contains(field)) {
        set(field, fields.integer_fields.at(field));
      } else if (fields.float_fields.contains(field)) {
        set(field, fields.float_fields.at(field));
      } else if (fields.string_fields.contains(field)) {
        set(field, fields.string_fields.at(field));
      }
    }
    return;
  }
  for (const auto& [name, value] : fields.integer_fields) {
    set(name, value);
  }
  for (const auto& [name, value] : fields.float_fields) {
    set(name, value);
  }
  for (const auto& [name, value] : fields.string_fields) {
    set(name, value);
  }
}

}  // namespace

VcfWriter::VcfWriter(const fs::path& vcf_file_path) : VcfWriter(vcf_file_path, VcfHeaderPtr()) {
}

VcfWriter::VcfWriter(const fs::path& vcf_file_path, VcfHeaderPtr header) : _hdr(std::move(header)) {
  if (vcf_file_path.string().ends_with(".gz")) {
    // Capture the path by value: this deleter runs when the shared_ptr is destroyed
    // (at VcfWriter destruction), which outlives the constructor's reference parameter.
    const auto close_and_index = [path = vcf_file_path](htsFile* const fp) {
      HtsFileDeleter{}(fp);
      if (!BcfIndexBuild(path, 0)) {
        Logging::Warn("Failed to build index for '{}'", path.string());
      }
    };
    // HtsOpen throws if the file cannot be opened; adopt into the shared handle with the indexing deleter.
    _file = AdoptShared(HtsOpen(vcf_file_path, "wz").release(), close_and_index);  // NOSONAR(M23_259)
  } else {
    _file = AdoptShared(HtsOpen(vcf_file_path, "w"));
  }
  if (_hdr == nullptr) {
    _hdr = VcfHeader::Create();
  }
  _file_path = vcf_file_path.string();
}

VcfWriter::VcfWriter(const fs::path& vcf_file_path, VcfHeaderPtr header, const NoIndex /*no_index*/)
    : _hdr(std::move(header)) {
  const char* const mode = vcf_file_path.string().ends_with(".gz") ? "wz" : "w";
  // HtsOpen throws if the file cannot be opened.
  _file = AdoptShared(HtsOpen(vcf_file_path, mode));
  if (_hdr == nullptr) {
    _hdr = VcfHeader::Create();
  }
  _file_path = vcf_file_path.string();
}

void VcfWriter::WriteHeader(const std::vector<std::string>& custom_meta_data_lines,
                            const std::vector<FilterFieldMetadata>& filter_lines,
                            const std::vector<InfoFieldMetadata>& info_lines,
                            const std::vector<FormatFieldMetadata>& format_lines,
                            const std::vector<ContigMetadata>& contig_lines,
                            const std::string& sample_name) const {
  _hdr->AddSample(sample_name);
  for (const auto& line : custom_meta_data_lines) {
    _hdr->AddCustomMetaDataLine(line);
  }
  for (const auto& field : info_lines) {
    _hdr->AddInfoLine(field);
  }
  for (const auto& line : filter_lines) {
    _hdr->AddFilterLine(line);
  }
  for (const auto& field : format_lines) {
    _hdr->AddFormatLine(field);
  }
  for (const auto& line : contig_lines) {
    _hdr->AddContigLine(line);
  }
  WriteHeader();
}

void VcfWriter::WriteHeader() const {
  _hdr->Write(_file);
}

VcfRecordPtr VcfWriter::CreateRecord(const std::string& chromosome,
                                     const s32 position,
                                     const std::string& id,
                                     const std::vector<std::string>& alleles,
                                     const std::optional<f32>& quality,
                                     const std::string& filter_name,
                                     const TypedVcfFields& info_fields,
                                     const TypedVcfFields& format_fields) {
  auto record = VcfRecord::CreateFromHeader(_hdr);

  record->SetChromosome(chromosome);
  record->SetPosition(position);
  record->SetId(id);
  record->SetAlleles(alleles);
  record->SetQuality(quality);
  record->SetFilter(filter_name);

  ApplyTypedFields(info_fields,
                   [&record](const std::string& name, const auto& value) { record->SetInfoField(name, value); });
  ApplyTypedFields(format_fields,
                   [&record](const std::string& name, const auto& value) { record->SetFormatField(name, value); });
  return record;
}

void VcfWriter::WriteRecord(const VcfRecordPtr& record) const {
  BcfWrite(_file.get(), _hdr->_hdr.get(), record->_record.get());
}

void VcfWriter::Flush() const {
  BcfFlush(_file.get());
}

void VcfWriter::WriteBgzfRecord(BGZF* const bgzf_file, const VcfHeaderPtr& hdr, const VcfRecordPtr& record) {
  const std::string formatted = VcfFormat(hdr->_hdr.get(), record->_record.get());
  // S5356: bgzf_write() requires const void* at its C API boundary; the implicit const char* → const void*
  // conversion is both safe and unavoidable.
  if (bgzf_write(bgzf_file, formatted.data(), formatted.size()) < 0) {  // NOSONAR(S5356)
    throw error::Error("Failed to write VCF record to BGZF");
  }
}

}  // namespace xoos::io
