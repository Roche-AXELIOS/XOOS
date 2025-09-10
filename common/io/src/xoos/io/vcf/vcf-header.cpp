#include "xoos/io/vcf/vcf-header.h"

#include <map>
#include <memory>

#include <fmt/format.h>

namespace xoos::io {

/**
 * @brief Convert FieldType to string representation.
 * @param type FieldType
 * @return String representation of the field type
 */
std::string FieldTypeToString(const FieldType type) {
  switch (type) {
    case FieldType::kInteger:
      return "Integer";
    case FieldType::kFloat:
      return "Float";
    case FieldType::kCharacter:
      return "Character";
    case FieldType::kString:
      return "String";
    case FieldType::kFlag:
      return "Flag";
  }
  return "Unknown";
}

/**
 * @brief Read a VCF header from the input file.
 * @param input_file HTS file handle
 */
VcfHeaderPtr VcfHeader::Read(const HtsFileSharedPtr& input_file) {
  auto* const vcf_hdr = bcf_hdr_read(input_file.get());
  if (vcf_hdr == nullptr) {
    throw std::runtime_error("Failed to read VCF header");
  }
  auto hdr = BcfHeaderPtr(vcf_hdr, bcf_hdr_destroy);
  return std::make_shared<VcfHeader>(hdr);
}

/**
 * @brief Create a new VCF header.
 * @return New VCF header
 */
VcfHeaderPtr VcfHeader::Create() {
  auto hdr = BcfHeaderPtr(bcf_hdr_init("w"), bcf_hdr_destroy);
  return std::make_shared<VcfHeader>(hdr);
}

/**
 * @brief Write the VCF header to the output file.
 * @param output_file HTS file handle
 */
void VcfHeader::Write(const HtsFileSharedPtr& output_file) {
  if (bcf_hdr_write(output_file.get(), _hdr.get()) < 0) {
    throw std::runtime_error("Failed to write VCF header");
  }
}

/**
 * @brief Clone the VCF header.
 * @return Cloned VCF header
 */
VcfHeaderPtr VcfHeader::Clone() const {
  auto* const vcf_hdr = bcf_hdr_dup(_hdr.get());
  if (vcf_hdr == nullptr) {
    throw std::runtime_error("Failed to clone VCF header");
  }
  auto hdr = BcfHeaderPtr(vcf_hdr, bcf_hdr_destroy);
  return std::make_shared<VcfHeader>(hdr);
}

/**
 * @brief Add a custom meta data line to the VCF header.
 * @param line Custom meta data line
 */
void VcfHeader::AddCustomMetaDataLine(const std::string& line) {
  if (bcf_hdr_append(_hdr.get(), line.c_str()) < 0) {
    throw std::runtime_error("Failed to add custom meta data line to VCF header");
  }
}

/**
 * @brief Add a sample to the VCF header.
 * @param sample_name Sample name
 */
void VcfHeader::AddSample(const std::string& sample_name) {
  if (bcf_hdr_add_sample(_hdr.get(), sample_name.c_str()) < 0) {
    throw std::runtime_error("Failed to add sample to VCF header");
  }
}

/**
 * @brief Add an INFO line to the VCF header.
 * @param info_line INFO field metadata
 */
void VcfHeader::AddInfoLine(const InfoFieldMetadata& info_line) {
  const std::string header_line = fmt::format("##INFO=<ID={},Number={},Type={},Description=\"{}\">",
                                              info_line.id,
                                              info_line.number,
                                              FieldTypeToString(info_line.type),
                                              info_line.description);
  if (bcf_hdr_append(_hdr.get(), header_line.c_str()) < 0) {
    throw std::runtime_error("Failed to add INFO line to VCF header");
  }
}

/**
 * @brief Add a FILTER line to the VCF header.
 * @param info_filter Filter field metadata
 */
void VcfHeader::AddFilterLine(const FilterFieldMetadata& info_filter) {
  const std::string header_line =
      fmt::format("##FILTER=<ID={},Description=\"{}\">", info_filter.id, info_filter.description);
  if (bcf_hdr_append(_hdr.get(), header_line.c_str()) < 0) {
    throw std::runtime_error("Failed to add FILTER line to VCF header");
  }
}

/**
 * @brief Add a FORMAT line to the VCF header.
 * @param format_line FORMAT field metadata
 */
void VcfHeader::AddFormatLine(const FormatFieldMetadata& format_line) {
  const std::string header_line = fmt::format("##FORMAT=<ID={},Number={},Type={},Description=\"{}\">",
                                              format_line.id,
                                              format_line.number,
                                              FieldTypeToString(format_line.type),
                                              format_line.description);
  if (bcf_hdr_append(_hdr.get(), header_line.c_str()) < 0) {
    throw std::runtime_error("Failed to add FORMAT line to VCF header");
  }
}

/**
 * @brief Add a contig line to the VCF header.
 * @param line Contig metadata
 */
void VcfHeader::AddContigLine(const ContigMetadata& line) {
  const std::string header_line = fmt::format("##contig=<ID={},length={}>", line.id, line.length);
  if (bcf_hdr_append(_hdr.get(), header_line.c_str()) < 0) {
    throw std::runtime_error("Failed to add contig line to VCF header");
  }
}

/**
 * Convert contig name to numeric ID.
 * @param ctg Contig name
 * @return Contig ID
 */
int VcfHeader::GetContigId(const std::string& ctg) {
  return bcf_hdr_name2id(_hdr.get(), ctg.c_str());
}

/**
 * @brief Generate a map of contig names and their order index.
 * @return Map of contig indexes
 */
std::map<std::string, s32> VcfHeader::GetContigIndexes() {
  std::map<std::string, s32> indexes;
  bcf_idpair_t* ctg = _hdr->id[BCF_DT_CTG];
  for (s32 i = 0; i < _hdr->n[BCF_DT_CTG]; ++i) {
    indexes[ctg[i].key] = i;
  }
  return indexes;
}

/**
 * @brief Generate a map of contig names and their lengths.
 * @return Map of contig lengths
 */
std::map<std::string, u64> VcfHeader::GetContigLengths() {
  std::map<std::string, u64> lengths{};
  bcf_idpair_t* ctg = _hdr->id[BCF_DT_CTG];
  for (s32 i = 0; i < _hdr->n[BCF_DT_CTG]; ++i) {
    lengths[ctg[i].key] = ctg[i].val->info[0];
  }
  return lengths;
}

/**
 * @brief Generate a map of sample names and their indexes.
 * @return Map of sample name indexes
 */
std::map<std::string, int> VcfHeader::GetSampleIndexes() {
  std::map<std::string, int> indexes{};
  bcf_idpair_t* samples = _hdr->id[BCF_DT_SAMPLE];
  for (s32 i = 0; i < _hdr->n[BCF_DT_SAMPLE]; ++i) {
    indexes[samples[i].key] = samples[i].val->id;
  }
  return indexes;
}

/**
 * @brief Return the value of specified key from the VCF header.
 * @param key The key
 * @return The value
 */
char* VcfHeader::GetValue(const char* key) {
  bcf_hrec_t* rec = bcf_hdr_get_hrec(_hdr.get(), BCF_HL_GEN, key, nullptr, nullptr);
  if (rec == nullptr) {
    return nullptr;
  }
  return rec->value;
}

bool VcfHeader::HasInfoField(const std::string& name) {
  return bcf_hdr_id2int(_hdr.get(), BCF_DT_ID, name.c_str()) >= 0;
}

/**
 * @brief Synchronize the header data structures after modifications.
 * This is necessary to ensure that the header is in a consistent state.
 */
void VcfHeader::Sync() {
  if (bcf_hdr_sync(_hdr.get()) < 0) {
    throw std::runtime_error("Failed to synchronize VCF header");
  }
}

}  // namespace xoos::io
