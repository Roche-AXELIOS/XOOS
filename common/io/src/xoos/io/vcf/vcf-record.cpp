#include "xoos/io/vcf/vcf-record.h"

#include <numeric>
#include <regex>
#include <stdexcept>

#include <htslib/vcf.h>

#include <xoos/io/malloc-ptr.h>

namespace xoos::io {

VcfRecord::VcfRecord(const VcfHeaderPtr& hdr) : _hdr(hdr->_hdr) {
}

/**
 * Create a new record from the header at the specified unpack level.
 * @param hdr Header
 * @param unpack Unpack level
 * @return New record
 */
VcfRecordPtr VcfRecord::CreateFromHeader(const VcfHeaderPtr& hdr, int unpack) {
  auto vcf_record = std::make_shared<VcfRecord>(hdr);
  const auto bcf_record = BcfRecordPtr(bcf_init(), bcf_destroy);
  if (bcf_record == nullptr) {
    throw std::runtime_error("Failed to initialize record");
  }
  vcf_record->_record = bcf_record;
  bcf_unpack(vcf_record->_record.get(), unpack);
  return vcf_record;
}

/**
 * Create a new record from the header at the maximum unpack level.
 * @param hdr Header
 * @return New record
 */
VcfRecordPtr VcfRecord::CreateFromHeader(const VcfHeaderPtr& hdr) {
  return CreateFromHeader(hdr, BCF_UN_ALL);
}

/**
 * Read record from file at the specified unpack level.
 * @param hdr Header
 * @param input_vcf_name HTS file handle
 * @param unpack Unpack level
 * @return Next record
 */
VcfRecordPtr VcfRecord::ReadFromFile(const VcfHeaderPtr& hdr, const HtsFileSharedPtr& input_vcf_name, int unpack) {
  const auto record = BcfRecordPtr{bcf_init(), bcf_destroy};
  if (record == nullptr) {
    throw std::runtime_error("Failed to initialize record");
  }
  if (const auto result = bcf_read(input_vcf_name.get(), hdr->_hdr.get(), record.get()); result < 0) {
    return nullptr;
  }
  auto shared_record = std::make_shared<VcfRecord>(hdr);
  shared_record->_record = record;
  bcf_unpack(record.get(), unpack);
  return shared_record;
}

/**
 * Read entire record from file.
 * @param hdr Header
 * @param input_vcf_name HTS file handle
 * @return Next record
 */
VcfRecordPtr VcfRecord::ReadFromFile(const VcfHeaderPtr& hdr, const HtsFileSharedPtr& input_vcf_name) {
  return ReadFromFile(hdr, input_vcf_name, BCF_UN_ALL);
}

/**
 * @brief Read the next record in the target region.
 * @param hdr Header
 * @param input_vcf_fp HTS file handle
 * @param tbx_idx Tabix index
 * @param hts_itr HTS iterator
 * @param kstr kstring struct
 * @param unpack Unpack level
 * @return Next record
 */
VcfRecordPtr VcfRecord::ReadFromRegion(const VcfHeaderPtr& hdr,
                                       const HtsFileSharedPtr& input_vcf_fp,
                                       tbx_t* tbx_idx,
                                       hts_itr_t* hts_itr,
                                       Kstring& kstr,
                                       const int unpack) {
  const auto record = BcfRecordPtr{bcf_init(), bcf_destroy};
  if (record == nullptr) {
    throw std::runtime_error("Failed to initialize record");
  }
  if (input_vcf_fp->format.format == bcf) {
    const auto bcf_itr_result = bcf_itr_next(input_vcf_fp.get(), hts_itr, record.get());
    if (bcf_itr_result < 0) {
      return nullptr;
    }
  } else {
    const auto tbx_itr_result = tbx_itr_next(input_vcf_fp.get(), tbx_idx, hts_itr, kstr.Get());
    if (tbx_itr_result < 0) {
      return nullptr;
    }
    if (const auto result = vcf_parse(kstr.Get(), hdr->_hdr.get(), record.get()); result < 0) {
      return nullptr;
    }
    kstr.Clear();
  }
  auto shared_record = std::make_shared<VcfRecord>(hdr);
  shared_record->_record = record;
  bcf_unpack(record.get(), unpack);
  return shared_record;
}

std::string VcfRecord::Chromosome() const {
  return {bcf_hdr_id2name(_hdr.get(), _record->rid)};
}

hts_pos_t VcfRecord::Position() const {
  return _record->pos;
}

bool VcfRecord::IsSnp() const {
  return bcf_is_snp(_record.get()) != 0;
}

std::string VcfRecord::Id() const {
  return {_record->d.id};
}

std::string VcfRecord::Allele(const int which) const {
  return {_record->d.allele[which]};
}

int VcfRecord::NumAlleles() const {
  return _record->n_allele;
}

void VcfRecord::SetChromosome(const std::string& chromosome) {
  _record->rid = bcf_hdr_name2id(_hdr.get(), chromosome.c_str());
}

void VcfRecord::SetPosition(const int position) {
  _record->pos = position;
}

void VcfRecord::SetId(const std::string& id) {
  if (bcf_update_id(_hdr.get(), _record.get(), id.c_str()) < 0) {
    throw std::runtime_error("Failed to update id " + id);
  }
}

/**
 * Set the alleles for the record.
 * @param alleles Alleles to set, including the reference allele
 */
void VcfRecord::SetAlleles(const std::vector<std::string>& alleles) {
  std::vector<const char*> tmp;
  tmp.reserve(alleles.size());
  for (const auto& value : alleles) {
    tmp.emplace_back(value.c_str());
  }
  if (bcf_update_alleles(_hdr.get(), _record.get(), &tmp[0], static_cast<int>(tmp.size())) < 0) {
    throw std::runtime_error("Failed to update alleles");
  }
}

/**
 * @brief Set the Quality score for the record.
 * @param quality Quality score, if std::nullopt, sets the quality to missing.
 */
void VcfRecord::SetQuality(const std::optional<float>& quality) {
  if (quality) {
    _record->qual = *quality;
  } else {
    bcf_float_set_missing(_record->qual);
  }
}

std::optional<float> VcfRecord::GetQuality() const {
  return (bcf_float_is_missing(_record->qual) != 0) ? std::nullopt : std::make_optional(_record->qual);
}

float VcfRecord::GetQuality(float default_value) const {
  return GetQuality().value_or(default_value);
}

/**
 * @brief Check whether an INFO field is in the record, without checking header.
 * @param field Field name
 * @return Whether field is found
 */
bool VcfRecord::HasInfoFieldNoCheck(const std::string& field) const {
  return bcf_get_info_id(_record.get(), bcf_hdr_id2int(_hdr.get(), BCF_DT_ID, field.c_str())) != nullptr;
}

/**
 * @brief Extract the float values of an INFO field, without checking header.
 * @param field Field name
 * @return Vector of float values
 */
template <>
std::vector<float> VcfRecord::GetInfoFieldNoCheck(const std::string& field) const {
  // look at record without checking header
  int num_entries = 0;
  MallocPtr<float> buffer = nullptr;
  if (bcf_get_info_values(
          _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_REAL) < 0) {
    return {};  // field not found
  }
  std::vector values(buffer.get(), buffer.get() + num_entries);
  return values;
}

/**
 * @brief Extract the int values of an INFO field, without checking header.
 * @param field Field name
 * @return Vector of int values
 */
template <>
std::vector<int> VcfRecord::GetInfoFieldNoCheck(const std::string& field) const {
  // look at record without checking header
  int num_entries = 0;
  MallocPtr<int> buffer = nullptr;
  if (bcf_get_info_values(
          _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_INT) < 0) {
    return {};  // field not found
  }
  std::vector values(buffer.get(), buffer.get() + num_entries);
  return values;
}

/**
 * @brief Extract the value string of an INFO field, without checking header.
 * @param field Field name
 * @return value string
 */
std::string VcfRecord::GetInfoFieldStringNoCheck(const std::string& field) const {
  // look at record without checking header
  int num_entries = 0;
  MallocPtr<char> buffer = nullptr;
  if (bcf_get_info_values(
          _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_STR) < 0) {
    return {};  // field not found
  }
  std::string value(buffer.get(), num_entries);

  // when we use bcf_read to set the bcf1_t record, it seems to add a trailing '\0' to format fields. However, when we
  // explicitly set the format fields, it does not add the trailing '\0'. So this may be a hacky way to check if there
  // is a '\0' at the end of the string and remove it if there is.
  if (!value.empty() && value.back() == '\0') {
    value.pop_back();
  }
  return value;
}

/**
 * @brief Check whether a FORMAT field is in the record, without checking header.
 * @param field Field name
 * @return Vector of int values
 */
bool VcfRecord::HasFormatFieldNoCheck(const std::string& field) const {
  return bcf_get_fmt_id(_record.get(), bcf_hdr_id2int(_hdr.get(), BCF_DT_ID, field.c_str())) != nullptr;
}

/**
 * @brief Extract the float values of a FORMAT field, without checking header.
 * @param field Field name
 * @return Vector of float values
 */
template <>
std::vector<float> VcfRecord::GetFormatFieldNoCheck(const std::string& field) const {
  // look at record without checking header
  int num_entries = 0;
  MallocPtr<float> buffer = nullptr;
  if (bcf_get_format_values(
          _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_REAL) < 0) {
    return {};
  }
  std::vector values(buffer.get(), buffer.get() + num_entries);
  return values;
}

/**
 * @brief Extract the int values of a FORMAT field, without checking header.
 * @param field Field name
 * @return Vector of int values
 */
template <>
std::vector<int> VcfRecord::GetFormatFieldNoCheck(const std::string& field) const {
  // look at record without checking header
  int num_entries = 0;
  MallocPtr<int> buffer = nullptr;
  if (bcf_get_format_values(
          _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_INT) < 0) {
    return {};
  }
  std::vector values(buffer.get(), buffer.get() + num_entries);
  return values;
}

/**
 * @brief Extract the value string of a FORMAT field, without checking header.
 * @param field Field name
 * @return value string
 */
std::string VcfRecord::GetFormatFieldStringNoCheck(const std::string& field) const {
  // look at record without checking header
  if (IsGtStandardFormatField(field)) {
    return {GetGTField()};
  }

  int num_entries = 0;
  MallocPtr<char> buffer = nullptr;
  if (bcf_get_format_values(
          _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_STR) < 0) {
    return {};
  }
  std::string value(buffer.get(), num_entries);

  // when we use bcf_read to set the bcf1_t record, it seems to add a trailing '\0' to format fields. However, when we
  // explicitly set the format fields, it does not add the trailing '\0'. So this may be a hacky way to check if there
  // is a '\0' at the end of the string and remove it if there is.
  if (!value.empty() && value.back() == '\0') {
    value.pop_back();
  }
  return value;
}

/**
 * @brief Extract the float values of an INFO field.
 * @param field Field name
 */
template <>
std::vector<float> VcfRecord::GetInfoField(const std::string& field) const {
  if (bcf_hdr_id2int(_hdr.get(), BCF_DT_ID, field.c_str()) < 0) {
    throw std::runtime_error("Info field " + field + " not found in header");
  }
  int num_entries = 0;
  MallocPtr<float> buffer = nullptr;
  if (bcf_get_info_values(
          _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_REAL) < 0) {
    throw std::runtime_error("Failed to get info field " + field);
  }
  std::vector values(buffer.get(), buffer.get() + num_entries);
  return values;
}

/**
 * @brief Extract the int values of an INFO field.
 * @param field Field name
 */
template <>
std::vector<int> VcfRecord::GetInfoField(const std::string& field) const {
  if (bcf_hdr_id2int(_hdr.get(), BCF_DT_ID, field.c_str()) < 0) {
    throw std::runtime_error("Info field " + field + " not found in header");
  }
  int num_entries = 0;
  MallocPtr<int> buffer = nullptr;
  if (bcf_get_info_values(
          _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_INT) < 0) {
    throw std::runtime_error("Failed to get info field " + field);
  }
  std::vector values(buffer.get(), buffer.get() + num_entries);
  return values;
}

/**
 * @brief Extract the string values of an INFO field.
 * @param field Field name
 */
template <>
std::vector<std::string> VcfRecord::GetInfoField(const std::string& field) const {
  if (bcf_hdr_id2int(_hdr.get(), BCF_DT_ID, field.c_str()) < 0) {
    throw std::runtime_error("Info field " + field + " not found in header");
  }
  int num_entries = 0;
  MallocPtr<char> buffer = nullptr;
  if (bcf_get_info_values(
          _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_STR) < 0) {
    throw std::runtime_error("Failed to get info field " + field);
  }
  std::string value(buffer.get(), num_entries);

  // when we use bcf_read to set the bcf1_t record, it seems to add a trailing '\0' to format fields. However, when we
  // explicitly set the format fields, it does not add the trailing '\0'. So this may be a hacky way to check if there
  // is a '\0' at the end of the string and remove it if there is.
  if (!value.empty() && value.back() == '\0') {
    value.pop_back();
  }
  return {value};
}

/**
 * @brief Extract the bool value of an INFO field.
 * For compatibility reasons, this will return a vector containing exactly one bool
 * @param field Field name
 */
template <>
std::vector<bool> VcfRecord::GetInfoField(const std::string& field) const {
  if (bcf_hdr_id2int(_hdr.get(), BCF_DT_ID, field.c_str()) < 0) {
    throw std::runtime_error("Info field " + field + " not found in header");
  }
  int num_entries = 0;
  MallocPtr<char> buffer = nullptr;  // this will be unused for retrieving bool
  // from htslib documentation:
  // Returns negative value on error or the number of values (including missing
  // values) put in *dst on success.
  // ...
  // bcf_get_info_flag() does not store anything in *dst but returns 1 if the
  // flag is set or 0 if not.
  return {(bcf_get_info_values(
               _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_FLAG) >
           0)};
}

/**
 * @brief Set the float values of an INFO field.
 * @param field Field name
 * @param values Vector of float values
 */
template <>
void VcfRecord::SetInfoField(const std::string& field, const std::vector<float>& values) {
  if (bcf_hdr_id2int(_hdr.get(), BCF_DT_ID, field.c_str()) < 0) {
    throw std::runtime_error("Info field " + field + " not found in header");
  }
  if (bcf_update_info_float(_hdr.get(), _record.get(), field.c_str(), &values[0], values.size()) < 0) {
    throw std::runtime_error("Failed to update info field " + field);
  }
}

/**
 * @brief Set the int values of an INFO field.
 * @param field Field name
 * @param values Vector of int values
 */
template <>
void VcfRecord::SetInfoField(const std::string& field, const std::vector<int>& values) {
  if (bcf_hdr_id2int(_hdr.get(), BCF_DT_ID, field.c_str()) < 0) {
    throw std::runtime_error("Info field " + field + " not found in header");
  }
  if (bcf_update_info_int32(_hdr.get(), _record.get(), field.c_str(), &values[0], values.size()) < 0) {
    throw std::runtime_error("Failed to update info field " + field);
  }
}

/**
 * @brief Set the string values of an INFO field.
 * @param field Field name
 * @param values Vector of string values
 */
template <>
void VcfRecord::SetInfoField(const std::string& field, const std::vector<std::string>& values) {
  const auto value = std::accumulate(values.begin(), values.end(), std::string{});
  if (bcf_hdr_id2int(_hdr.get(), BCF_DT_ID, field.c_str()) < 0) {
    throw std::runtime_error("Info field " + field + " not found in header");
  }
  if (bcf_update_info_string(_hdr.get(), _record.get(), field.c_str(), value.c_str()) < 0) {
    throw std::runtime_error("Failed to update info field " + field);
  }
}

/**
 * @brief Get the float value of an FORMAT field.
 * @param field Field name
 * @return Vector of float values
 */
template <>
std::vector<float> VcfRecord::GetFormatField(const std::string& field) const {
  if (bcf_hdr_id2int(_hdr.get(), BCF_DT_ID, field.c_str()) < 0) {
    throw std::runtime_error("Format field " + field + " not found in header");
  }
  int num_entries = 0;
  MallocPtr<float> buffer = nullptr;
  if (bcf_get_format_values(
          _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_REAL) < 0) {
    throw std::runtime_error("Failed to get format field " + field);
  }
  std::vector values(buffer.get(), buffer.get() + num_entries);
  return values;
}

/**
 * @brief Get the int value of an FORMAT field.
 * @param field Field name
 * @return Vector of int values
 */
template <>
std::vector<int> VcfRecord::GetFormatField(const std::string& field) const {
  if (bcf_hdr_id2int(_hdr.get(), BCF_DT_ID, field.c_str()) < 0) {
    throw std::runtime_error("Format field " + field + " not found in header");
  }
  int num_entries = 0;
  MallocPtr<int> buffer = nullptr;
  if (bcf_get_format_values(
          _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_INT) < 0) {
    throw std::runtime_error("Failed to get format field " + field);
  }
  std::vector values(buffer.get(), buffer.get() + num_entries);
  return values;
}

/**
 * @brief Get the string value of an FORMAT field.
 * @param field Field name
 * @return Vector of string values
 */
template <>
std::vector<std::string> VcfRecord::GetFormatField(const std::string& field) const {
  if (bcf_hdr_id2int(_hdr.get(), BCF_DT_ID, field.c_str()) < 0) {
    throw std::runtime_error("Format field " + field + " not found in header");
  }

  if (IsGtStandardFormatField(field)) {
    return {GetGTField()};
  }

  int num_entries = 0;
  MallocPtr<char> buffer = nullptr;
  if (bcf_get_format_values(
          _hdr.get(), _record.get(), field.c_str(), reinterpret_cast<void**>(&buffer), &num_entries, BCF_HT_STR) < 0) {
    throw std::runtime_error("Failed to get format field " + field);
  }
  std::string value(buffer.get(), num_entries);

  // when we use bcf_read to set the bcf1_t record, it seems to add a trailing '\0' to format fields. However, when we
  // explicitly set the format fields, it does not add the trailing '\0'. So this may be a hacky way to check if there
  // is a '\0' at the end of the string and remove it if there is.
  if (!value.empty() && value.back() == '\0') {
    value.pop_back();
  }
  return {value};
}

/**
 * @brief Set the float values of a FORMAT field.
 * @param field Field name
 * @param values Vector of float values
 */
template <>
void VcfRecord::SetFormatField(const std::string& field, const std::vector<float>& values) {
  if (bcf_hdr_id2int(_hdr.get(), BCF_DT_ID, field.c_str()) < 0) {
    throw std::runtime_error("Format field " + field + " not found in header");
  }
  if (bcf_update_format_float(_hdr.get(), _record.get(), field.c_str(), &values[0], values.size()) < 0) {
    throw std::runtime_error("Failed to update format field " + field);
  }
}

/**
 * @brief Set the int values of a FORMAT field.
 * @param field Field name
 * @param values Vector of int values
 */
template <>
void VcfRecord::SetFormatField(const std::string& field, const std::vector<int>& values) {
  if (bcf_hdr_id2int(_hdr.get(), BCF_DT_ID, field.c_str()) < 0) {
    throw std::runtime_error("Format field " + field + " not found in header");
  }
  if (bcf_update_format_int32(_hdr.get(), _record.get(), field.c_str(), &values[0], values.size()) < 0) {
    throw std::runtime_error("Failed to update format field " + field);
  }
}

/**
 * @brief Set the string values of a FORMAT field.
 * @param field Field name
 * @param values Vector of C-string values
 * @note If the field is GT, it will be set using SetGTField instead.
 */
template <>
void VcfRecord::SetFormatField(const std::string& field, const std::vector<const char*>& values) {
  if (bcf_hdr_id2int(_hdr.get(), BCF_DT_ID, field.c_str()) < 0) {
    throw std::runtime_error("Format field " + field + " not found in header");
  }

  if (IsGtStandardFormatField(field)) {
    SetGTField(values[0]);
    return;
  }

  std::vector tmp(values);
  if (bcf_update_format_string(_hdr.get(), _record.get(), field.c_str(), &tmp[0], static_cast<int>(tmp.size())) < 0) {
    throw std::runtime_error("Failed to update format field " + field);
  }
}

/**
 * @brief Set the string values of a FORMAT field.
 * @param field Field name
 * @param values Vector of string values
 * @note If the field is GT, it will be set using SetGTField instead.
 */
template <>
void VcfRecord::SetFormatField(const std::string& field, const std::vector<std::string>& values) {
  if (bcf_hdr_id2int(_hdr.get(), BCF_DT_ID, field.c_str()) < 0) {
    throw std::runtime_error("Format field " + field + " not found in header");
  }

  if (IsGtStandardFormatField(field)) {
    SetGTField(values[0]);
    return;
  }

  std::vector<const char*> tmp;
  tmp.reserve(values.size());
  for (const auto& value : values) {
    tmp.emplace_back(value.c_str());
  }
  if (bcf_update_format_string(_hdr.get(), _record.get(), field.c_str(), &tmp[0], static_cast<int>(tmp.size())) < 0) {
    throw std::runtime_error("Failed to update format field " + field);
  }
}

/**
 * @brief Set the genotype field (GT) in the record.
 * @param value Genotype string, e.g. "0/1" or "1|0"
 * @throws std::runtime_error if the value is invalid
 */
void VcfRecord::SetGTField(const std::string& value) {
  const std::regex pattern(kGenotypeRegex);
  std::smatch match;

  if (!std::regex_search(value, match, pattern)) {
    throw std::runtime_error("Invalid genotype field " + value);
  }

  std::vector<int> genotypes;

  const bool is_phased = match[3].str() == "|";

  // set the first genotype
  if (match[1].str() == ".") {
    genotypes.push_back(bcf_gt_missing);
  } else {
    genotypes.push_back(is_phased ? bcf_gt_phased(std::stoi(match[1].str()))
                                  : bcf_gt_unphased(std::stoi(match[1].str())));
  }

  // set the second genotype if it exists
  if (match[3].matched) {
    if (match[4].str() == ".") {
      genotypes.push_back(bcf_gt_missing);
    } else {
      genotypes.push_back(is_phased ? bcf_gt_phased(std::stoi(match[4].str()))
                                    : bcf_gt_unphased(std::stoi(match[4].str())));
    }
  }

  if (bcf_update_genotypes(_hdr.get(), _record.get(), &genotypes[0], genotypes.size()) < 0) {
    throw std::runtime_error("Failed to update format field GT");
  }
}

/**
 * @brief Get the genotype field (GT) from the record.
 * @return Genotype string, e.g. "0/1" or "1|0"
 * @throws std::runtime_error if the GT field cannot be retrieved
 */
std::string VcfRecord::GetGTField() const {
  int num_entries = 0;
  MallocPtr<int> buffer = nullptr;
  if (bcf_get_genotypes(_hdr.get(), _record.get(), &buffer, &num_entries) < 0) {
    throw std::runtime_error("Failed to get format field GT");
  }

  std::string value;
  if (bcf_gt_is_missing(buffer.get()[0])) {
    value = ".";
  } else {
    value = std::to_string(bcf_gt_allele(buffer.get()[0]));
  }

  if (num_entries > 1) {
    if (bcf_gt_is_missing(buffer.get()[1])) {
      value += "/.";
    } else if (bcf_gt_is_phased(buffer.get()[1])) {
      value += "|" + std::to_string(bcf_gt_allele(buffer.get()[1]));
    } else {
      value += "/" + std::to_string(bcf_gt_allele(buffer.get()[1]));
    }
  }
  return value;
}

/**
 * @brief Set the exact filter field in the record.
 * @param filter Filter name
 * @throws std::runtime_error if the filter is not found in the header or if updating fails
 */
void VcfRecord::SetFilter(const std::string& filter) {
  int index = bcf_hdr_id2int(_hdr.get(), BCF_DT_ID, filter.c_str());
  if (index < 0) {
    throw std::runtime_error("Filter " + filter + " not found in header");
  }
  if (bcf_update_filter(_hdr.get(), _record.get(), &index, 1) < 0) {
    throw std::runtime_error("Failed to update filter " + filter);
  }
}

/**
 * @brief Add a filter to the record.
 * @param filter Filter name
 * @throws std::runtime_error if the filter is not found in the header or if adding fails
 */
void VcfRecord::AddFilter(const std::string& filter) {
  const int index = bcf_hdr_id2int(_hdr.get(), BCF_DT_ID, filter.c_str());
  if (index < 0) {
    throw std::runtime_error("Filter " + filter + " not found in header");
  }
  if (bcf_add_filter(_hdr.get(), _record.get(), index) < 0) {
    throw std::runtime_error("Failed to add filter " + filter);
  }
}

/**
 * @brief Remove all filters from the record.
 * @throws std::runtime_error if clearing filters fails
 */
void VcfRecord::ClearFilters() {
  if (bcf_update_filter(_hdr.get(), _record.get(), nullptr, 0) < 0) {
    throw std::runtime_error("Failed to clear filters");
  }
}

/**
 * @brief Get the filters applied to the record.
 * @return Vector of filter names
 */
std::vector<std::string> VcfRecord::GetFilters() const {
  std::vector<std::string> filters;

  if (_record->d.n_flt > 0) {
    filters.reserve(_record->d.n_flt);
    for (int i = 0; i < _record->d.n_flt; i++) {
      filters.emplace_back(bcf_hdr_int2id(_hdr, BCF_DT_ID, _record->d.flt[i]));
    }
  }
  return filters;
}

/**
 * @brief Clones record with a new header
 * @param new_header
 * @return
 */
VcfRecordPtr VcfRecord::Clone(const VcfHeaderPtr& new_header) const {
  const auto record = BcfRecordPtr{bcf_dup(_record.get()), bcf_destroy};
  if (record == nullptr) {
    throw std::runtime_error("Failed to clone record");
  }
  auto shared_record = std::make_shared<VcfRecord>(new_header);
  shared_record->_record = record;
  return shared_record;
}

}  // namespace xoos::io
