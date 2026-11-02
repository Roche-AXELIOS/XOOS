#include "xoos/io/vcf/vcf-record.h"

#include <numeric>
#include <regex>

#include <htslib/vcf.h>

#include <xoos/error/error.h>
#include <xoos/io/htslib-util/htslib-bcf-util.h>
#include <xoos/util/adopt.h>

#include "xoos/io/vcf/constants.h"

namespace xoos::io {

namespace {
/// Adopt a wrapper-owned Bcf1Ptr (unique_ptr) into the module's shared BcfRecordPtr.
BcfRecordPtr ToSharedRecord(Bcf1Ptr record) {
  return AdoptShared(std::move(record));
}
}  // namespace

VcfRecord::VcfRecord(const VcfHeaderPtr& hdr) : _header(hdr), _hdr(hdr->_hdr) {
}

/**
 * Create a new record from the header at the specified unpack level.
 * @param hdr Header
 * @param unpack Unpack level
 * @return New record
 */
VcfRecordPtr VcfRecord::CreateFromHeader(const VcfHeaderPtr& hdr, const s32 unpack) {
  auto vcf_record = std::make_shared<VcfRecord>(hdr);
  vcf_record->_record = ToSharedRecord(BcfInit());
  BcfUnpack(vcf_record->_record.get(), unpack);
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
VcfRecordPtr VcfRecord::ReadFromFile(const VcfHeaderPtr& hdr,
                                     const HtsFileSharedPtr& input_vcf_name,
                                     const s32 unpack) {
  const auto record = ToSharedRecord(BcfInit());
  // BcfRead returns false at EOF and throws on a parse error (result < -1).
  if (!BcfRead(input_vcf_name.get(), hdr->_hdr.get(), record.get())) {
    return nullptr;
  }
  auto shared_record = std::make_shared<VcfRecord>(hdr);
  shared_record->_record = record;
  BcfUnpack(record.get(), unpack);
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
                                       tbx_t* const tbx_idx,
                                       hts_itr_t* const hts_itr,
                                       Kstring& kstr,
                                       const s32 unpack) {
  const auto record = ToSharedRecord(BcfInit());
  if (input_vcf_fp->format.format == bcf) {
    // BcfItrNext returns false at end-of-iterator and throws on a read/decode error.
    if (!BcfItrNext(input_vcf_fp.get(), hts_itr, record.get())) {
      return nullptr;
    }
  } else {
    // TbxItrNext returns false at end-of-iterator and throws on a read/decode error.
    if (!TbxItrNext(input_vcf_fp.get(), tbx_idx, hts_itr, kstr.Get())) {
      return nullptr;
    }
    if (!VcfParse(kstr.Get(), hdr->_hdr.get(), record.get())) {
      return nullptr;
    }
    kstr.Clear();
  }
  auto shared_record = std::make_shared<VcfRecord>(hdr);
  shared_record->_record = record;
  BcfUnpack(record.get(), unpack);
  return shared_record;
}

std::string VcfRecord::Chromosome() const {
  return {BcfHdrId2Name(_hdr.get(), _record->rid)};
}

hts_pos_t VcfRecord::Position() const {
  return _record->pos;
}

bool VcfRecord::IsPass() const {
  return BcfHasFilter(_hdr.get(), _record.get(), kFilterPASS) == 1 ||
         BcfHasFilter(_hdr.get(), _record.get(), kFilterDOT) == 1;
}

bool VcfRecord::IsSnp() const {
  return BcfIsSnp(_record.get());
}

std::string VcfRecord::Id() const {
  return {_record->d.id};
}

std::string VcfRecord::Allele(const s32 which) const {
  // htslib stores alleles in a C array; indexing it requires pointer arithmetic.
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  return {_record->d.allele[which]};
}

s32 VcfRecord::NumAlleles() const {
  // Returns s32 to match the signed integer convention used throughout the codebase.
  return static_cast<s32>(_record->n_allele);
}

void VcfRecord::SetChromosome(const std::string& chromosome) const {
  _record->rid = BcfHdrName2Id(_hdr.get(), chromosome);
}

void VcfRecord::SetPosition(const s32 position) const {
  _record->pos = position;
}

void VcfRecord::SetId(const std::string& id) const {
  BcfUpdateId(_hdr.get(), _record.get(), id);
}

/**
 * Set the alleles for the record.
 * @param alleles Alleles to set, including the reference allele
 */
void VcfRecord::SetAlleles(const std::vector<std::string>& alleles) const {
  BcfUpdateAlleles(_hdr.get(), _record.get(), alleles);
}

/**
 * @brief Set the Quality score for the record.
 * @param quality Quality score, if std::nullopt, sets the quality to missing.
 */
void VcfRecord::SetQuality(const std::optional<f32>& quality) const {
  if (quality.has_value()) {
    _record->qual = *quality;
  } else {
    bcf_float_set_missing(_record->qual);
  }
}

std::optional<f32> VcfRecord::GetQuality() const {
  return (bcf_float_is_missing(_record->qual) != 0) ? std::nullopt : std::make_optional(_record->qual);
}

f32 VcfRecord::GetQuality(f32 default_value) const {
  return GetQuality().value_or(default_value);
}

/**
 * @brief Check whether an INFO field is in the record, without checking header.
 * @param field Field name
 * @return Whether field is found
 */
bool VcfRecord::HasInfoFieldNoCheck(const std::string& field) const {
  return BcfRecordHasInfo(_record.get(), BcfHdrId2Int(_hdr.get(), BCF_DT_ID, field));
}

/**
 * @brief Extract the float values of an INFO field, without checking header.
 * @param field Field name
 * @return Vector of float values
 */
// S5536: explicit template specialization; must live in one .cpp to avoid multiple
// definitions — called by downstream modules that link this library.
template <>  // NOSONAR(S5536)
std::vector<f32> VcfRecord::GetInfoFieldNoCheck(const std::string& field) const {
  // look at record without checking header
  return BcfGetInfoValues<f32>(_hdr.get(), _record.get(), field).value_or(std::vector<f32>{});
}

/**
 * @brief Extract the int values of an INFO field, without checking header.
 * @param field Field name
 * @return Vector of int values
 */
// S5536: explicit template specialization; must live in one .cpp to avoid multiple
// definitions — called by downstream modules that link this library.
template <>  // NOSONAR(S5536)
std::vector<s32> VcfRecord::GetInfoFieldNoCheck(const std::string& field) const {
  // look at record without checking header
  return BcfGetInfoValues<s32>(_hdr.get(), _record.get(), field).value_or(std::vector<s32>{});
}

/**
 * @brief Extract the value string of an INFO field, without checking header.
 * @param field Field name
 * @return value string
 */
std::string VcfRecord::GetInfoFieldStringNoCheck(const std::string& field) const {
  // look at record without checking header. BcfGetInfoString already strips any trailing NUL byte.
  return BcfGetInfoString(_hdr.get(), _record.get(), field).value_or(std::string{});
}

/**
 * @brief Check whether a FORMAT field is in the record, without checking header.
 * @param field Field name
 * @return Vector of int values
 */
bool VcfRecord::HasFormatFieldNoCheck(const std::string& field) const {
  return BcfRecordHasFormat(_record.get(), BcfHdrId2Int(_hdr.get(), BCF_DT_ID, field));
}

/**
 * @brief Extract the float values of a FORMAT field, without checking header.
 * @param field Field name
 * @return Vector of float values
 */
// S5536: explicit template specialization; must live in one .cpp to avoid multiple
// definitions — called by downstream modules that link this library.
template <>  // NOSONAR(S5536)
std::vector<f32> VcfRecord::GetFormatFieldNoCheck(const std::string& field) const {
  // look at record without checking header
  return BcfGetFormatValues<f32>(_hdr.get(), _record.get(), field).value_or(std::vector<f32>{});
}

/**
 * @brief Extract the int values of a FORMAT field, without checking header.
 * @param field Field name
 * @return Vector of int values
 */
// S5536: explicit template specialization; must live in one .cpp to avoid multiple
// definitions — called by downstream modules that link this library.
template <>  // NOSONAR(S5536)
std::vector<s32> VcfRecord::GetFormatFieldNoCheck(const std::string& field) const {
  // look at record without checking header
  return BcfGetFormatValues<s32>(_hdr.get(), _record.get(), field).value_or(std::vector<s32>{});
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

  // BcfGetFormatString already strips any trailing NUL byte that htslib appends.
  return BcfGetFormatString(_hdr.get(), _record.get(), field).value_or(std::string{});
}

/**
 * @brief Confirm that a field is in the header.
 * @param field Field name
 * @param field_type Field type, e.g. INFO or FORMAT
 */
void VcfRecord::ConfirmFieldInHeader(const std::string& field, const std::string& field_type) const {
  if (BcfHdrId2Int(_hdr.get(), BCF_DT_ID, field) < 0) {
    throw error::Error("{} field {} not found in header", field_type, field);
  }
}

/**
 * @brief Confirm that an INFO field is in the header.
 * @param field Field name
 */
void VcfRecord::ConfirmInfoFieldInHeader(const std::string& field) const {
  ConfirmFieldInHeader(field, "INFO");
}

/**
 * @brief Confirm that a FORMAT field is in the header.
 * @param field Field name
 */
void VcfRecord::ConfirmFormatFieldInHeader(const std::string& field) const {
  ConfirmFieldInHeader(field, "FORMAT");
}

/**
 * @brief Extract the float values of an INFO field.
 * @param field Field name
 */
// S5536: explicit template specialization; must live in one .cpp to avoid multiple
// definitions — called by downstream modules that link this library.
template <>  // NOSONAR(S5536)
std::vector<f32> VcfRecord::GetInfoField(const std::string& field) const {
  ConfirmInfoFieldInHeader(field);
  auto values = BcfGetInfoValues<f32>(_hdr.get(), _record.get(), field);
  if (!values) {
    throw error::Error("Failed to get INFO field {}", field);
  }
  return *values;
}

/**
 * @brief Extract the int values of an INFO field.
 * @param field Field name
 */
// S5536: explicit template specialization; must live in one .cpp to avoid multiple
// definitions — called by downstream modules that link this library.
template <>  // NOSONAR(S5536)
std::vector<s32> VcfRecord::GetInfoField(const std::string& field) const {
  ConfirmInfoFieldInHeader(field);
  auto values = BcfGetInfoValues<s32>(_hdr.get(), _record.get(), field);
  if (!values) {
    throw error::Error("Failed to get INFO field {}", field);
  }
  return *values;
}

/**
 * @brief Extract the string values of an INFO field.
 * @param field Field name
 */
// S5536: explicit template specialization; must live in one .cpp to avoid multiple
// definitions — called by downstream modules that link this library.
template <>  // NOSONAR(S5536)
std::vector<std::string> VcfRecord::GetInfoField(const std::string& field) const {
  ConfirmInfoFieldInHeader(field);
  // BcfGetInfoString already strips any trailing NUL byte that htslib appends.
  auto value = BcfGetInfoString(_hdr.get(), _record.get(), field);
  if (!value) {
    throw error::Error("Failed to get INFO field {}", field);
  }
  return {*value};
}

/**
 * @brief Extract the bool value of an INFO field.
 * For compatibility reasons, this will return a vector containing exactly one bool
 * @param field Field name
 */
template <>  // NOSONAR(M23_357) — std::vector<bool> bit-packing is intentional here
std::vector<bool> VcfRecord::GetInfoField(const std::string& field) const {
  ConfirmInfoFieldInHeader(field);
  // bcf_get_info_flag() returns 1 if the flag is set or 0 if not; BcfGetInfoFlag wraps this.
  return {BcfGetInfoFlag(_hdr.get(), _record.get(), field)};
}

/**
 * @brief Set the float values of an INFO field.
 * @param field Field name
 * @param values Vector of float values
 */
void VcfRecord::SetInfoField(const std::string& field, const std::vector<f32>& values) {
  ConfirmInfoFieldInHeader(field);
  BcfUpdateInfoFloat(_hdr.get(), _record.get(), field, values);
}

/**
 * @brief Set the int values of an INFO field.
 * @param field Field name
 * @param values Vector of int values
 */
void VcfRecord::SetInfoField(const std::string& field, const std::vector<s32>& values) {
  ConfirmInfoFieldInHeader(field);
  BcfUpdateInfoInt(_hdr.get(), _record.get(), field, values);
}

/**
 * @brief Set the string values of an INFO field.
 * @param field Field name
 * @param values Vector of string values
 */
void VcfRecord::SetInfoField(const std::string& field, const std::vector<std::string>& values) {
  ConfirmInfoFieldInHeader(field);
  const auto value = std::accumulate(values.begin(), values.end(), std::string{});
  BcfUpdateInfoString(_hdr.get(), _record.get(), field, value);
}

/**
 * @brief Add a flag INFO field to the record.
 * @param field Field name
 */
void VcfRecord::AddInfoFieldFlag(const std::string& field) {
  ConfirmInfoFieldInHeader(field);
  // set to true to add this flag
  BcfUpdateInfoFlag(_hdr.get(), _record.get(), field, true);
}

/**
 * @brief Remove a flag INFO field from the record.
 * @param field Field name
 */
void VcfRecord::RemoveInfoFieldFlag(const std::string& field) {
  ConfirmInfoFieldInHeader(field);
  // set to false to remove this flag
  BcfUpdateInfoFlag(_hdr.get(), _record.get(), field, false);
}

/**
 * @brief Get the float value of an FORMAT field.
 * @param field Field name
 * @return Vector of float values
 */
// S5536: explicit template specialization; must live in one .cpp to avoid multiple
// definitions — called by downstream modules that link this library.
template <>  // NOSONAR(S5536)
std::vector<f32> VcfRecord::GetFormatField(const std::string& field) const {
  ConfirmFormatFieldInHeader(field);
  auto values = BcfGetFormatValues<f32>(_hdr.get(), _record.get(), field);
  if (!values) {
    throw error::Error("Failed to get FORMAT field {}", field);
  }
  return *values;
}

/**
 * @brief Get the int value of an FORMAT field.
 * @param field Field name
 * @return Vector of int values
 */
// S5536: explicit template specialization; must live in one .cpp to avoid multiple
// definitions — called by downstream modules that link this library.
template <>  // NOSONAR(S5536)
std::vector<s32> VcfRecord::GetFormatField(const std::string& field) const {
  ConfirmFormatFieldInHeader(field);
  auto values = BcfGetFormatValues<s32>(_hdr.get(), _record.get(), field);
  if (!values) {
    throw error::Error("Failed to get FORMAT field {}", field);
  }
  return *values;
}

/**
 * @brief Get the string value of an FORMAT field.
 * @param field Field name
 * @return Vector of string values
 */
// S5536: explicit template specialization; must live in one .cpp to avoid multiple
// definitions — called by downstream modules that link this library.
template <>  // NOSONAR(S5536)
std::vector<std::string> VcfRecord::GetFormatField(const std::string& field) const {
  ConfirmFormatFieldInHeader(field);

  if (IsGtStandardFormatField(field)) {
    return {GetGTField()};
  }

  // BcfGetFormatString already strips any trailing NUL byte that htslib appends.
  auto value = BcfGetFormatString(_hdr.get(), _record.get(), field);
  if (!value) {
    throw error::Error("Failed to get FORMAT field {}", field);
  }
  return {*value};
}

/**
 * @brief Set the float values of a FORMAT field.
 * @param field Field name
 * @param values Vector of float values
 */
void VcfRecord::SetFormatField(const std::string& field, const std::vector<f32>& values) {
  ConfirmFormatFieldInHeader(field);
  BcfUpdateFormatFloat(_hdr.get(), _record.get(), field, values);
}

/**
 * @brief Set the int values of a FORMAT field.
 * @param field Field name
 * @param values Vector of int values
 */
void VcfRecord::SetFormatField(const std::string& field, const std::vector<s32>& values) {
  ConfirmFormatFieldInHeader(field);
  BcfUpdateFormatInt(_hdr.get(), _record.get(), field, values);
}

/**
 * @brief Set the string values of a FORMAT field.
 * @param field Field name
 * @param values Vector of C-string values
 * @note If the field is GT, it will be set using SetGTField instead.
 */
void VcfRecord::SetFormatField(const std::string& field, const std::vector<const char*>& values) {
  ConfirmFormatFieldInHeader(field);

  if (IsGtStandardFormatField(field)) {
    SetGTField(values[0]);
    return;
  }

  BcfUpdateFormatString(_hdr.get(), _record.get(), field, std::vector<std::string>(values.begin(), values.end()));
}

/**
 * @brief Set the string values of a FORMAT field.
 * @param field Field name
 * @param values Vector of string values
 * @note If the field is GT, it will be set using SetGTField instead.
 */
void VcfRecord::SetFormatField(const std::string& field, const std::vector<std::string>& values) {
  ConfirmFormatFieldInHeader(field);

  if (IsGtStandardFormatField(field)) {
    SetGTField(values[0]);
    return;
  }

  BcfUpdateFormatString(_hdr.get(), _record.get(), field, values);
}

/**
 * @brief Parse a VCF GT string into a vector of htslib GT integers.
 * @param value Genotype string, e.g. "0/1" or "1|0"
 * @return Vector of htslib-encoded allele integers (one per allele in the genotype)
 * @throws error::Error if the value does not match the expected genotype format
 */
std::vector<s32> VcfRecord::ParseGenotypeString(const std::string& value) {
  static const std::regex kPattern(kGenotypeRegex);
  if (!std::regex_match(value, kPattern)) {
    throw error::Error("Invalid genotype field value: {}", value);
  }

  // Detect separator: phased (|) or unphased (/)
  const bool is_phased = value.contains(kPhasedAlleleSeparator);
  const char sep = is_phased ? kPhasedAlleleSeparator : kUnphasedAlleleSeparator;

  // Convert a single allele token to an htslib GT integer.
  auto parse_allele = [&is_phased](const std::string& allele_str) {
    if (allele_str == ".") {
      return bcf_gt_missing;
    }
    const s32 allele = std::stoi(allele_str);
    return is_phased ? bcf_gt_phased(allele) : bcf_gt_unphased(allele);
  };

  // Walk the GT string, splitting on the detected separator to build the allele integer vector.
  std::vector<s32> alleles;
  std::string token;
  for (const char ch : value) {
    if (ch == sep) {
      alleles.push_back(parse_allele(token));
      token.clear();
    } else {
      token += ch;
    }
  }
  alleles.push_back(parse_allele(token));
  return alleles;
}

/**
 * @brief Build the full multi-sample GT buffer for bcf_update_genotypes.
 * @param alleles The target sample's new alleles (htslib-encoded integers)
 * @param sample_idx Target sample index (0-based), assumed already validated
 * @return Flat genotype buffer sized num_samples * buf_ploidy, ready for bcf_update_genotypes
 * @note Preserves other samples' existing genotypes and pads each sample to the buffer ploidy
 *       with bcf_int32_vector_end so htslib serializes shorter genotypes correctly.
 */
std::vector<s32> VcfRecord::BuildGenotypeBuffer(const std::vector<s32>& alleles, const size_t sample_idx) const {
  const auto num_samples_u = static_cast<size_t>(BcfHdrNumSamples(_hdr.get()));
  const size_t new_ploidy = alleles.size();

  // Read the existing GT buffer (needed to preserve other samples' genotypes).
  const auto existing_gt_values = BcfGetGenotypes(_hdr.get(), _record.get());
  const auto existing_num_gts = existing_gt_values ? static_cast<s32>(existing_gt_values->size()) : 0;
  const bool has_existing = existing_num_gts > 0;

  // Use the larger of the existing ploidy and the new ploidy as the buffer ploidy. This avoids
  // silently truncating other samples' genotype data when the new GT has fewer alleles.
  const size_t existing_ploidy = has_existing ? static_cast<size_t>(existing_num_gts) / num_samples_u : new_ploidy;
  const size_t buf_ploidy = std::max(existing_ploidy, new_ploidy);
  const size_t total_gts = num_samples_u * buf_ploidy;

  // Initialize all positions to vector_end so padding slots beyond each sample's true ploidy
  // are handled correctly by htslib (they serialize as absent rather than as missing alleles).
  std::vector<s32> genotypes(total_gts, bcf_int32_vector_end);

  // Copy existing genotypes for all samples except the target, preserving their full allele vectors.
  if (has_existing) {
    for (size_t s = 0; s < num_samples_u; ++s) {
      if (s == sample_idx) {
        continue;
      }
      for (size_t a = 0; a < existing_ploidy; ++a) {
        genotypes[(s * buf_ploidy) + a] = (*existing_gt_values)[(s * existing_ploidy) + a];
      }
      // Positions [existing_ploidy, buf_ploidy) remain bcf_int32_vector_end.
    }
  }

  // Write the target sample's new alleles.
  for (size_t i = 0; i < new_ploidy; ++i) {
    genotypes[(sample_idx * buf_ploidy) + i] = alleles[i];
  }
  // Place a vector_end sentinel after the last allele when the new ploidy is shorter than the
  // buffer ploidy, so htslib serializes this sample as e.g. "0/1" instead of "0/1/.".
  if (new_ploidy < buf_ploidy) {
    genotypes[(sample_idx * buf_ploidy) + new_ploidy] = bcf_int32_vector_end;
  }

  return genotypes;
}

/**
 * @brief Set the genotype field (GT) in the record for a specific sample.
 * @param value Genotype string, e.g. "0/1" or "1|0"
 * @param sample_index Sample index (0-based)
 * @throws error::Error if the value is invalid or sample_index is out of range
 */
// Setter mutates *_record via htslib bcf_update_genotypes; the pointee mutation is invisible to the
// analyzer because it is reached through the const _record.get() accessor, hence the suppression.
void VcfRecord::SetGTField(const std::string& value, const s32 sample_index) {  // NOSONAR(S5817)
  const std::vector<s32> alleles = ParseGenotypeString(value);

  const s32 num_samples = BcfHdrNumSamples(_hdr.get());
  if (sample_index < 0 || sample_index >= num_samples) {
    throw error::Error("Sample index out of range: {}", sample_index);
  }

  const std::vector<s32> genotypes = BuildGenotypeBuffer(alleles, static_cast<size_t>(sample_index));

  BcfUpdateGenotypes(_hdr.get(), _record.get(), genotypes);
}

/**
 * @brief Set the genotype field (GT) in the record.
 * @param value Genotype string, e.g. "0/1" or "1|0"
 * @throws error::Error if the value is invalid
 */
void VcfRecord::SetGTField(const std::string& value) {
  SetGTField(value, 0);
}

/**
 * @brief Get the genotype field (GT) from the record.
 * @return Genotype string, e.g. "0/1" or "1|0"
 * @throws error::Error if the GT field cannot be retrieved
 */
std::string VcfRecord::GetGTField() const {
  return GetGTField(0);
}

/**
 * @brief Get the genotype field (GT) for a specific sample.
 * @param sample_index Sample index (0-based)
 * @return Genotype string, e.g. "0/1" or "1|0"
 * @throws error::Error if the GT field cannot be retrieved or sample_index is out of range
 */
std::string VcfRecord::GetGTField(const s32 sample_index) const {
  const s32 num_samples = BcfHdrNumSamples(_hdr.get());
  if (sample_index < 0 || sample_index >= num_samples) {
    throw error::Error("Sample index out of range: {}", sample_index);
  }

  auto genotypes = BcfGetGenotypes(_hdr.get(), _record.get());
  if (!genotypes) {
    throw error::Error("Failed to get FORMAT field GT");
  }
  const auto& gts = *genotypes;
  const auto num_entries = static_cast<s32>(gts.size());

  const auto ploidy = static_cast<size_t>(num_entries / num_samples);
  const auto offset = static_cast<size_t>(sample_index) * ploidy;

  std::string value;
  for (size_t i = 0; i < ploidy; ++i) {
    if (gts.at(offset + i) == bcf_int32_vector_end) {
      break;
    }
    if (i > 0) {
      value += bcf_gt_is_phased(gts.at(offset + i)) ? kPhasedAlleleSeparator : kUnphasedAlleleSeparator;
    }
    if (static_cast<bool>(bcf_gt_is_missing(gts.at(offset + i)))) {
      value += kMissingAllele;
    } else {
      value += std::to_string(bcf_gt_allele(gts.at(offset + i)));
    }
  }
  return value;
}

/**
 * @brief Get the genotype field (GT) for a specific sample.
 * @param filter Filter name
 * @throws error::Error if the filter is not found in the header or if updating fails
 */
void VcfRecord::SetFilter(const char* filter) const {
  const s32 index = BcfHdrId2Int(_hdr.get(), BCF_DT_ID, filter);
  if (index < 0) {
    throw error::Error("FILTER {} not found in header", filter);
  }
  BcfUpdateFilter(_hdr.get(), _record.get(), {index});
}

/**
 * @brief Set the filter for the record.
 * @param filter Filter name
 * @throws error::Error if the filter is not found in the header or if updating fails
 */
void VcfRecord::SetFilter(const std::string& filter) const {
  SetFilter(filter.c_str());
}

/**
 * @brief Set the filter for the record.
 * @param filter Filter name
 * @throws error::Error if the filter is not found in the header or if updating fails
 */
void VcfRecord::SetFilter(const std::string_view filter) const {
  const std::string filter_str(filter);
  SetFilter(filter_str.c_str());
}

/**
 * @brief Add a filter to the record.
 * @param filter Filter name
 * @throws error::Error if the filter is not found in the header or if adding fails
 */
void VcfRecord::AddFilter(const char* filter) const {
  const s32 index = BcfHdrId2Int(_hdr.get(), BCF_DT_ID, filter);
  if (index < 0) {
    throw error::Error("FILTER {} not found in header", filter);
  }
  BcfAddFilter(_hdr.get(), _record.get(), index);
}

/**
 * @brief Add a filter to the record.
 * @param filter Filter name
 * @throws error::Error if the filter is not found in the header or if adding fails
 */
void VcfRecord::AddFilter(const std::string& filter) const {
  AddFilter(filter.c_str());
}

/**
 * @brief Add a filter to the record.
 * @param filter Filter name
 * @throws error::Error if the filter is not found in the header or if adding fails
 */
void VcfRecord::AddFilter(const std::string_view filter) const {
  const std::string filter_str(filter);
  AddFilter(filter_str.c_str());
}

/**
 * @brief Remove all filters from the record.
 * @throws error::Error if clearing filters fails
 */
void VcfRecord::ClearFilters() const {
  BcfUpdateFilter(_hdr.get(), _record.get(), {});
}

/**
 * @brief Get the filters applied to the record.
 * @return Vector of filter names
 */
std::vector<std::string> VcfRecord::GetFilters() const {
  std::vector<std::string> filters;

  if (_record->d.n_flt > 0) {
    filters.reserve(static_cast<size_t>(_record->d.n_flt));
    for (s32 i = 0; i < _record->d.n_flt; ++i) {
      // htslib stores filter ids in a C array; indexing it requires pointer arithmetic.
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      filters.emplace_back(BcfHdrInt2Id(_hdr.get(), BCF_DT_ID, _record->d.flt[i]));
    }
  }
  return filters;
}

/**
 * @brief Clones record with a new header
 * @param new_header
 * @return Cloned record
 */
VcfRecordPtr VcfRecord::Clone(const VcfHeaderPtr& new_header) const {
  auto shared_record = std::make_shared<VcfRecord>(new_header);
  shared_record->_record = ToSharedRecord(BcfDup(_record.get()));
  _header->TranslateRecord(new_header, shared_record->_record.get());
  return shared_record;
}

}  // namespace xoos::io
