#include "xoos/io/vcf/vcf-header.h"

#include <strings.h>

#include <memory>
#include <mutex>
#include <shared_mutex>
#include <string>
#include <tuple>
#include <utility>

#include <fmt/format.h>
#include <htslib/vcf.h>

#include <xoos/error/error.h>
#include <xoos/io/htslib-util/htslib-bcf-util.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>
#include <xoos/util/adopt.h>

namespace xoos::io {

/**
 * @brief Synchronization state shared by all wrappers around one htslib header.
 *
 * VcfHeader accepts shared BCF header pointers, so multiple VcfHeader instances can
 * alias the same mutable bcf_hdr_t. The revision identifies the header layout used
 * by cached translation contexts; the mutex prevents a translation from observing a
 * header while it is being modified.
 */
struct VcfHeaderState {
  std::shared_mutex mutex;
  u64 revision = 0;
};

namespace {
/// Adopt a wrapper-owned BcfHdrPtr (unique_ptr) into the module's shared BcfHeaderPtr.
BcfHeaderPtr ToSharedHeader(BcfHdrPtr hdr) {
  return AdoptShared(std::move(hdr));
}

/**
 * @brief Return the synchronization state associated with an underlying BCF header.
 *
 * Weak keys avoid keeping headers alive solely because they were seen by this cache.
 * owner_less makes aliases that share the same control block use one state.
 * @param header Underlying BCF header whose shared state is requested.
 * @return Shared synchronization state for @p header.
 */
std::shared_ptr<VcfHeaderState> GetHeaderState(const BcfHeaderPtr& header) {
  using HeaderStateMap =
      std::map<std::weak_ptr<bcf_hdr_t>, std::weak_ptr<VcfHeaderState>, std::owner_less<std::weak_ptr<bcf_hdr_t>>>;
  static std::mutex states_mutex;
  static HeaderStateMap states;

  std::scoped_lock lock(states_mutex);
  std::erase_if(states, [](const auto& entry) { return entry.first.expired(); });
  const std::weak_ptr<bcf_hdr_t> key{header};
  const auto it = states.find(key);
  if (it != states.end()) {
    if (auto state = it->second.lock()) {
      return state;
    }
  }

  auto state = std::make_shared<VcfHeaderState>();
  states.insert_or_assign(key, state);
  return state;
}
}  // namespace

/**
 * @brief Construct a VCF header wrapper and attach shared synchronization state.
 * @param hdr Shared htslib header to wrap.
 */
VcfHeader::VcfHeader(BcfHeaderPtr hdr) : _hdr(std::move(hdr)), _state(GetHeaderState(_hdr)) {
}

VcfHeaderTranslationContext::VcfHeaderTranslationContext(const bcf_hdr_t* const source_header,
                                                         const u64 source_revision,
                                                         const u64 destination_revision)
    : source_header_snapshot(BcfHdrDup(source_header)),
      source_revision_snapshot(source_revision),
      destination_revision_snapshot(destination_revision) {
}

/**
 * @brief Convert FieldType to string representation.
 * @param type FieldType
 * @return String representation of the field type
 */
std::string FieldTypeToString(const FieldType type) {
  using enum FieldType;
  switch (type) {
    case kInteger:
      return "Integer";
    case kFloat:
      return "Float";
    case kCharacter:
      return "Character";
    case kString:
      return "String";
    case kFlag:
      return "Flag";
    default:
      break;
  }
  return "Unknown";
}

/**
 * @brief Read a VCF header from the input file.
 * @param input_file HTS file handle
 */
VcfHeaderPtr VcfHeader::Read(const HtsFileSharedPtr& input_file) {
  auto hdr = ToSharedHeader(BcfHdrRead(input_file.get()));
  return std::make_shared<VcfHeader>(hdr);
}

/**
 * @brief Create a new VCF header.
 * @return New VCF header
 */
VcfHeaderPtr VcfHeader::Create() {
  auto hdr = ToSharedHeader(BcfHdrInit("w"));
  return std::make_shared<VcfHeader>(hdr);
}

/**
 * @brief Write the VCF header to the output file.
 * @param output_file HTS file handle
 */
void VcfHeader::Write(const HtsFileSharedPtr& output_file) const {
  BcfHdrWrite(output_file.get(), _hdr.get());
}

/**
 * @brief Clone the VCF header.
 * @return Cloned VCF header
 */
VcfHeaderPtr VcfHeader::Clone() const {
  auto hdr = ToSharedHeader(BcfHdrDup(_hdr.get()));
  return std::make_shared<VcfHeader>(hdr);
}

void VcfHeader::TranslateRecord(const VcfHeaderPtr& destination, bcf1_t* const record) const {
  // Records already use the same htslib dictionary when both wrappers share a header state.
  // Skip translation and its cache/synchronization overhead in this common clone path.
  if (_state == destination->_state) {
    return;
  }

  // Shared locks prevent a context from being built against a layout that is being changed.
  const auto translate = [&]() {
    const auto context = GetTranslationContext(destination, _state->revision, destination->_state->revision);
    std::scoped_lock context_lock(context->mutex);
    if (bcf_translate(destination->_hdr.get(), context->source_header_snapshot.get(), record) != 0) {
      throw error::Error("Unable to translate VCF record to destination header");
    }
  };

  std::shared_lock source_lock(_state->mutex, std::defer_lock);
  std::shared_lock destination_lock(destination->_state->mutex, std::defer_lock);
  std::lock(source_lock, destination_lock);
  translate();
}

VcfHeaderTranslationContextPtr VcfHeader::GetTranslationContext(const VcfHeaderPtr& destination,
                                                                const u64 source_revision,
                                                                const u64 destination_revision) const {
  std::scoped_lock lock(_translation_contexts_mutex);
  std::erase_if(_translation_contexts, [](const auto& entry) { return entry.first.expired(); });
  const std::weak_ptr<VcfHeader> destination_key{destination};
  auto it = _translation_contexts.find(destination_key);
  if (it == _translation_contexts.end() || it->second->source_revision_snapshot != source_revision ||
      it->second->destination_revision_snapshot != destination_revision) {
    it = _translation_contexts
             .insert_or_assign(
                 destination_key,
                 std::make_shared<VcfHeaderTranslationContext>(_hdr.get(), source_revision, destination_revision))
             .first;
  }
  return it->second;
}

void VcfHeader::InvalidateTranslationContexts() const {
  // Mutators hold the state's exclusive lock before reaching this function. Incrementing
  // before changing the htslib header ensures any context created afterward observes a new
  // source/destination revision.
  ++_state->revision;
  std::scoped_lock lock(_translation_contexts_mutex);
  _translation_contexts.clear();
}

/**
 * @brief Add a custom meta data line to the VCF header.
 * @param line Custom meta data line
 */
void VcfHeader::AddCustomMetaDataLine(const std::string& line) {
  std::unique_lock lock(_state->mutex);
  InvalidateTranslationContexts();
  BcfHdrAppend(_hdr.get(), line);
}

/**
 * @brief Add a sample to the VCF header.
 * @param sample_name Sample name
 */
void VcfHeader::AddSample(const std::string& sample_name) {
  std::unique_lock lock(_state->mutex);
  InvalidateTranslationContexts();
  BcfHdrAddSample(_hdr.get(), sample_name);
}

/**
 * @brief Add an INFO line to the VCF header.
 * @param info_line INFO field metadata
 */
void VcfHeader::AddInfoLine(const InfoFieldMetadata& info_line) {
  std::unique_lock lock(_state->mutex);
  InvalidateTranslationContexts();
  const std::string header_line = fmt::format("##INFO=<ID={},Number={},Type={},Description=\"{}\">",
                                              info_line.id,
                                              info_line.number,
                                              FieldTypeToString(info_line.type),
                                              info_line.description);
  BcfHdrAppend(_hdr.get(), header_line);
}

/**
 * @brief Add a FILTER line to the VCF header.
 * @param info_filter Filter field metadata
 */
void VcfHeader::AddFilterLine(const FilterFieldMetadata& info_filter) {
  std::unique_lock lock(_state->mutex);
  InvalidateTranslationContexts();
  const std::string header_line =
      fmt::format("##FILTER=<ID={},Description=\"{}\">", info_filter.id, info_filter.description);
  BcfHdrAppend(_hdr.get(), header_line);
}

/**
 * @brief Add a FORMAT line to the VCF header.
 * @param format_line FORMAT field metadata
 */
void VcfHeader::AddFormatLine(const FormatFieldMetadata& format_line) {
  std::unique_lock lock(_state->mutex);
  InvalidateTranslationContexts();
  const std::string header_line = fmt::format("##FORMAT=<ID={},Number={},Type={},Description=\"{}\">",
                                              format_line.id,
                                              format_line.number,
                                              FieldTypeToString(format_line.type),
                                              format_line.description);
  BcfHdrAppend(_hdr.get(), header_line);
}

/**
 * @brief Add a contig line to the VCF header.
 * @param line Contig metadata
 */
void VcfHeader::AddContigLine(const ContigMetadata& line) {
  std::unique_lock lock(_state->mutex);
  InvalidateTranslationContexts();
  const std::string header_line = fmt::format("##contig=<ID={},length={}>", line.id, line.length);
  BcfHdrAppend(_hdr.get(), header_line);
}

/**
 * @brief Remove one INFO meta-information line from the header by field ID.
 * @param id INFO field ID to remove
 */
void VcfHeader::RemoveInfoLine(const std::string& id) const {
  std::unique_lock lock(_state->mutex);
  InvalidateTranslationContexts();
  BcfHdrRemove(_hdr.get(), BCF_HL_INFO, id.c_str());
}

/**
 * @brief Remove all INFO meta-information lines from the header.
 */
void VcfHeader::RemoveAllInfoLines() const {
  std::unique_lock lock(_state->mutex);
  InvalidateTranslationContexts();
  BcfHdrRemoveAll(_hdr.get(), BCF_HL_INFO);
}

/**
 * @brief Remove all FORMAT meta-information lines from the header.
 */
void VcfHeader::RemoveAllFormatLines() const {
  std::unique_lock lock(_state->mutex);
  InvalidateTranslationContexts();
  BcfHdrRemoveAll(_hdr.get(), BCF_HL_FMT);
}

/**
 * @brief Remove all FILTER meta-information lines from the header.
 */
void VcfHeader::RemoveAllFilterLines() const {
  std::unique_lock lock(_state->mutex);
  InvalidateTranslationContexts();
  BcfHdrRemoveAll(_hdr.get(), BCF_HL_FLT);
}

/**
 * Convert contig name to numeric ID.
 * @param ctg Contig name
 * @return Contig ID
 */
s32 VcfHeader::GetContigId(const std::string& ctg) const {
  return BcfHdrName2Id(_hdr.get(), ctg);
}

/**
 * @brief Generate a map of contig names and their order index.
 * @return Map of contig indexes
 */
StrMap<s32> VcfHeader::GetContigIndexes() const {
  StrMap<s32> indexes;
  const bcf_idpair_t* const ctg = _hdr->id[BCF_DT_CTG];
  for (s32 i = 0; i < _hdr->n[BCF_DT_CTG]; ++i) {
    // Indexing into htslib's contig id-pair array.
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    indexes[ctg[i].key] = i;
  }
  return indexes;
}

/**
 * @brief Generate a map of contig names and their lengths.
 * @return Map of contig lengths
 */
StrMap<u64> VcfHeader::GetContigLengths() const {
  StrMap<u64> lengths{};
  const bcf_idpair_t* const ctg = _hdr->id[BCF_DT_CTG];
  for (s32 i = 0; i < _hdr->n[BCF_DT_CTG]; ++i) {
    // Indexing into htslib's contig id-pair array.
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    lengths[ctg[i].key] = ctg[i].val->info[0];
  }
  return lengths;
}

/**
 * @brief Generate a map of sample names and their indexes.
 * @return Map of sample name indexes
 */
StrMap<s32> VcfHeader::GetSampleIndexes() const {
  StrMap<s32> indexes{};
  const bcf_idpair_t* const samples = _hdr->id[BCF_DT_SAMPLE];
  for (s32 i = 0; i < _hdr->n[BCF_DT_SAMPLE]; ++i) {
    // Indexing into htslib's sample id-pair array.
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    indexes[samples[i].key] = samples[i].val->id;
  }
  return indexes;
}

/**
 * @brief Return the value of specified key from the VCF header.
 * @param key The key
 * @return The value
 */
char* VcfHeader::GetValue(const char* const key) const {
  bcf_hrec_t* const rec = BcfHdrGetHrec(_hdr.get(), key);
  if (rec == nullptr) {
    return nullptr;
  }
  return rec->value;
}

bool VcfHeader::HasInfoField(const std::string& name) const {
  return BcfHdrId2Int(_hdr.get(), BCF_DT_ID, name) >= 0;
}

bool VcfHeader::HasFilterField(const std::string& name) const {
  const s32 id = BcfHdrId2Int(_hdr.get(), BCF_DT_ID, name);
  if (id < 0 || id >= _hdr->n[BCF_DT_ID]) {
    return false;
  }

  // Indexing into htslib's INFO/FORMAT/FILTER id-pair array.
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  const bcf_idinfo_t* const id_info = _hdr->id[BCF_DT_ID][id].val;
  return id_info != nullptr && id_info->hrec[BCF_HL_FLT] != nullptr;
}

bool VcfHeader::HasFormatField(const std::string& name) const {
  const s32 id = BcfHdrId2Int(_hdr.get(), BCF_DT_ID, name);
  if (id < 0 || id >= _hdr->n[BCF_DT_ID]) {
    return false;
  }

  // htslib stores INFO, FORMAT, and FILTER IDs in a shared dictionary. Inspect
  // the FORMAT header record so INFO-only IDs are not reported as FORMAT fields.
  // Indexing into htslib's INFO/FORMAT id-pair array.
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  const bcf_idinfo_t* const id_info = _hdr->id[BCF_DT_ID][id].val;
  return id_info != nullptr && id_info->hrec[BCF_HL_FMT] != nullptr;
}

/**
 * @brief Synchronize the header data structures after modifications.
 * This is necessary to ensure that the header is in a consistent state.
 */
void VcfHeader::Sync() const {
  std::unique_lock lock(_state->mutex);
  InvalidateTranslationContexts();
  BcfHdrSync(_hdr.get());
}

namespace {

/**
 * @brief Convert string representation of field type to FieldType enum.
 * @param str String representation of the field type
 * @return FieldType enum
 */
FieldType FieldTypeStringToEnum(const std::string& str) {
  using enum FieldType;
  if (strcasecmp(str.c_str(), "Integer") == 0) {
    return kInteger;
  }
  if (strcasecmp(str.c_str(), "Float") == 0) {
    return kFloat;
  }
  if (strcasecmp(str.c_str(), "Character") == 0) {
    return kCharacter;
  }
  if (strcasecmp(str.c_str(), "Flag") == 0) {
    return kFlag;
  }
  return kString;
}

/**
 * @brief Get the value of a key from a BCF header record.
 * @param hrec BCF header record
 * @param key Key to look for
 * @param default_value Default value if key is not found
 * @return Value associated with the key, or empty string if not found
 */
std::string BcfHrecGetValue(bcf_hrec_t* const hrec, const std::string& key, const std::string& default_value) {
  const s32 idx = BcfHrecFindKey(hrec, key);
  if (idx < 0) {
    return default_value;
  }
  // htslib stores header-record values in a C array; indexing it requires pointer arithmetic.
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  return hrec->vals[idx];
}

/**
 * @brief Get the value of a key from a BCF header record.
 * @param hrec BCF header record
 * @param key Key to look for
 * @return Value associated with the key, or empty string if not found
 */
std::string BcfHrecGetValue(bcf_hrec_t* const hrec, const std::string& key) {
  return BcfHrecGetValue(hrec, key, "");
}

}  // namespace

std::tuple<vec<InfoFieldMetadata>, vec<FormatFieldMetadata>> VcfHeader::GetFieldMetadata() const {
  vec<InfoFieldMetadata> info_metadata;
  vec<FormatFieldMetadata> fmt_metadata;
  // Iterate over all INFO / FORMAT fields in the header
  for (s32 i = 0; i < _hdr->n[BCF_DT_ID]; ++i) {
    // Indexing into htslib's INFO/FORMAT id-pair array.
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    const bcf_idpair_t* const idpair = &_hdr->id[BCF_DT_ID][i];
    const bcf_idinfo_t* const idinfo = idpair->val;
    if (idinfo->hrec[BCF_HL_INFO] != nullptr) {
      // extract metadata for the INFO field
      bcf_hrec_t* const hrec = idinfo->hrec[BCF_HL_INFO];
      const InfoFieldMetadata metadata{
          .id = BcfHrecGetValue(hrec, "ID"),
          .description = BcfHrecGetValue(hrec, "Description"),
          .number = BcfHrecGetValue(hrec, "Number"),
          .type = FieldTypeStringToEnum(BcfHrecGetValue(hrec, "Type")),
      };
      info_metadata.push_back(metadata);
    } else if (idinfo->hrec[BCF_HL_FMT] != nullptr) {
      // extract metadata for the FORMAT field
      bcf_hrec_t* const hrec = idinfo->hrec[BCF_HL_FMT];
      const FormatFieldMetadata metadata{
          .id = BcfHrecGetValue(hrec, "ID"),
          .description = BcfHrecGetValue(hrec, "Description"),
          .number = BcfHrecGetValue(hrec, "Number"),
          .type = FieldTypeStringToEnum(BcfHrecGetValue(hrec, "Type")),
      };
      fmt_metadata.push_back(metadata);
    }
  }
  return std::make_tuple(info_metadata, fmt_metadata);
}

vec<FilterFieldMetadata> VcfHeader::GetFilterMetadata() const {
  vec<FilterFieldMetadata> filter_metadata;
  for (s32 i = 0; i < _hdr->n[BCF_DT_ID]; ++i) {
    // Indexing into htslib's INFO/FORMAT/FILTER id-pair array.
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    const bcf_idpair_t* const id_pair = &_hdr->id[BCF_DT_ID][i];
    const bcf_idinfo_t* const id_info = id_pair->val;
    if (id_info->hrec[BCF_HL_FLT] == nullptr) {
      continue;
    }

    bcf_hrec_t* const hrec = id_info->hrec[BCF_HL_FLT];
    filter_metadata.push_back({
        .id = BcfHrecGetValue(hrec, "ID"),
        .description = BcfHrecGetValue(hrec, "Description"),
    });
  }
  return filter_metadata;
}

}  // namespace xoos::io
