#include "xoos/io/vcf/vcf-writer.h"

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <htslib/vcf.h>

#include <xoos/log/logging.h>

namespace xoos::io {

VcfWriter::VcfWriter(const fs::path& vcf_file_path, VcfHeaderPtr header) : _file(nullptr), _hdr(std::move(header)) {
  if (vcf_file_path.string().ends_with(".gz")) {
    const auto close_and_index = [&vcf_file_path](const auto& fp) {
      hts_close(fp);
      if (bcf_index_build(vcf_file_path.c_str(), 0) != 0) {
        Logging::Warn("Failed to build index for '{}'", vcf_file_path.string());
      }
    };
    _file = HtsFileSharedPtr{hts_open(vcf_file_path.c_str(), "wz"), close_and_index};
  } else {
    _file = HtsFileSharedPtr{hts_open(vcf_file_path.c_str(), "w"), hts_close};
  }
  if (_hdr == nullptr) {
    _hdr = VcfHeader::Create();
  }
  _file_path = vcf_file_path.string();

  if (_file == nullptr) {
    throw std::runtime_error("Could not open VCF file " + vcf_file_path.string());
  }
}

/**
 * @brief Writes a VCF header to the internal header object file.
 * @param custom_meta_data_lines
 * @param filter_lines
 * @param info_lines
 * @param format_lines
 * @param contig_lines
 * @param sample_name
 */
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

/**
 * @brief Create a new VCF record with the specified fields.
 * @param chromosome Chromosome name
 * @param position Genomic position (1-based)
 * @param id Variant ID (e.g., rsID)
 * @param alleles Vector of alleles, including the reference allele as the first element
 * @param quality Optional quality score; if std::nullopt, sets the quality to missing
 * @param filter_name Filter name to set for the record
 * @param info_fields Info fields to set in the record, including integer, float, and string fields
 * @param format_fields Format fields to set in the record, including integer, float, and string fields
 * @return A VcfRecordPtr to the created VcfRecord
 */
VcfRecordPtr VcfWriter::CreateRecord(const std::string& chromosome,
                                     const int position,
                                     const std::string& id,
                                     const std::vector<std::string>& alleles,
                                     const std::optional<float>& quality,
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

  if (info_fields.field_order.has_value()) {
    for (const auto& field : info_fields.field_order.value()) {
      if (info_fields.integer_fields.contains(field)) {
        record->SetInfoField(field, info_fields.integer_fields.at(field));
      } else if (info_fields.float_fields.contains(field)) {
        record->SetInfoField(field, info_fields.float_fields.at(field));
      } else if (info_fields.string_fields.contains(field)) {
        record->SetInfoField(field, info_fields.string_fields.at(field));
      }
    }
  } else {
    for (const auto& [fst, snd] : info_fields.integer_fields) {
      record->SetInfoField(fst, snd);
    }
    for (const auto& [fst, snd] : info_fields.float_fields) {
      record->SetInfoField(fst, snd);
    }
    for (const auto& [fst, snd] : info_fields.string_fields) {
      record->SetInfoField(fst, snd);
    }
  }

  if (format_fields.field_order.has_value()) {
    for (const auto& field : format_fields.field_order.value()) {
      if (format_fields.integer_fields.contains(field)) {
        record->SetFormatField(field, format_fields.integer_fields.at(field));
      } else if (format_fields.float_fields.contains(field)) {
        record->SetFormatField(field, format_fields.float_fields.at(field));
      } else if (format_fields.string_fields.contains(field)) {
        record->SetFormatField(field, format_fields.string_fields.at(field));
      }
    }
  } else {
    for (const auto& [fst, snd] : format_fields.integer_fields) {
      record->SetFormatField(fst, snd);
    }
    for (const auto& [fst, snd] : format_fields.float_fields) {
      record->SetFormatField(fst, snd);
    }
    for (const auto& [fst, snd] : format_fields.string_fields) {
      record->SetFormatField(fst, snd);
    }
  }
  return record;
}

void VcfWriter::WriteRecord(const VcfRecordPtr& record) const {
  if (bcf_write(_file.get(), _hdr->_hdr.get(), record->_record.get()) < 0) {
    throw std::runtime_error("Failed to write VCF record");
  }
}

void VcfWriter::Flush() {
  bcf_flush(_file.get());
}

}  // namespace xoos::io
