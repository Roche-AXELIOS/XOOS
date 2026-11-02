#include "htslib-bcf-util.h"

#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>

#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>

#include <xoos/error/error.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "xoos/io/htslib-util/htslib-bcf-ptr.h"
#include "xoos/io/htslib-util/htslib-ptr.h"
#include "xoos/io/malloc-ptr.h"

namespace xoos::io {

// --- Header (bcf_hdr_*) ---

BcfHdrPtr BcfHdrRead(htsFile* const fp) {
  BcfHdrPtr hdr{bcf_hdr_read(fp)};
  if (hdr == nullptr) {
    throw error::Error("Failed to read VCF/BCF header");
  }
  return hdr;
}

BcfHdrPtr BcfHdrInit(const std::string& mode) {
  BcfHdrPtr hdr{bcf_hdr_init(mode.c_str())};
  if (hdr == nullptr) {
    throw error::Error("Failed to initialize VCF/BCF header");
  }
  return hdr;
}

BcfHdrPtr BcfHdrDup(const bcf_hdr_t* const hdr) {
  BcfHdrPtr dup{bcf_hdr_dup(hdr)};
  if (dup == nullptr) {
    throw error::Error("Failed to duplicate VCF/BCF header");
  }
  return dup;
}

void BcfHdrWrite(htsFile* const fp, bcf_hdr_t* const hdr) {
  if (bcf_hdr_write(fp, hdr) < 0) {
    throw error::Error("Failed to write VCF/BCF header");
  }
}

void BcfHdrAppend(bcf_hdr_t* const hdr, const std::string& line) {
  if (bcf_hdr_append(hdr, line.c_str()) < 0) {
    throw error::Error("Failed to append line to VCF/BCF header: {}", line);
  }
}

void BcfHdrAddSample(bcf_hdr_t* const hdr, const std::string& sample) {
  if (bcf_hdr_add_sample(hdr, sample.c_str()) < 0) {
    throw error::Error("Failed to add sample to VCF/BCF header: {}", sample);
  }
}

void BcfHdrSync(bcf_hdr_t* const hdr) {
  if (bcf_hdr_sync(hdr) < 0) {
    throw error::Error("Failed to synchronize VCF/BCF header");
  }
}

void BcfHdrRemove(bcf_hdr_t* const hdr, const s32 type, const char* const key) {
  bcf_hdr_remove(hdr, type, key);
}

void BcfHdrRemoveAll(bcf_hdr_t* const hdr, const s32 type) {
  bcf_hdr_remove(hdr, type, nullptr);
}

s32 BcfHdrNumSamples(const bcf_hdr_t* const hdr) {
  return bcf_hdr_nsamples(hdr);
}

s32 BcfHdrName2Id(const bcf_hdr_t* const hdr, const std::string& name) {
  return bcf_hdr_name2id(hdr, name.c_str());
}

const char* BcfHdrId2Name(const bcf_hdr_t* const hdr, const s32 rid) {
  // bcf_hdr_id2name is a macro performing unchecked array indexing, so an out-of-range rid is
  // undefined behavior (segfault) rather than a null return. Bounds-check against the contig
  // dictionary before calling into htslib.
  if (rid < 0 || rid >= hdr->n[BCF_DT_CTG]) {
    throw error::Error("Contig id {} is out of range for VCF/BCF header", rid);
  }
  const auto* const name = bcf_hdr_id2name(hdr, rid);
  if (name == nullptr) {
    throw error::Error("Failed to get contig name for id {} from VCF/BCF header", rid);
  }
  return name;
}

s32 BcfHdrId2Int(const bcf_hdr_t* const hdr, const s32 dict_type, const std::string& id) {
  return bcf_hdr_id2int(hdr, dict_type, id.c_str());
}

const char* BcfHdrInt2Id(const bcf_hdr_t* const hdr, const s32 dict_type, const s32 int_id) {
  // bcf_hdr_int2id is a macro performing unchecked array indexing, so an out-of-range dict_type or
  // int_id is undefined behavior (segfault) rather than a null return. Bounds-check against the
  // selected dictionary before calling into htslib. htslib's dictionary types are BCF_DT_ID,
  // BCF_DT_CTG, BCF_DT_SAMPLE, so the dictionary array size is BCF_DT_SAMPLE + 1.
  const s32 num_dict_types = BCF_DT_SAMPLE + 1;
  if (dict_type < 0 || dict_type >= num_dict_types) {
    throw error::Error("Dictionary type {} is out of range for VCF/BCF header", dict_type);
  }
  // dict_type is bounds-checked against num_dict_types just above, so indexing hdr->n is safe.
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
  if (int_id < 0 || int_id >= hdr->n[dict_type]) {
    throw error::Error("Id {} is out of range for VCF/BCF header dictionary {}", int_id, dict_type);
  }
  const auto* const id = bcf_hdr_int2id(hdr, dict_type, int_id);
  if (id == nullptr) {
    throw error::Error("Failed to get id name for int id {} from VCF/BCF header", int_id);
  }
  return id;
}

bcf_hrec_t* BcfHdrGetHrec(const bcf_hdr_t* const hdr, const std::string& key) {
  return bcf_hdr_get_hrec(hdr, BCF_HL_GEN, key.c_str(), nullptr, nullptr);
}

s32 BcfHrecFindKey(bcf_hrec_t* const hrec, const std::string& key) {
  return bcf_hrec_find_key(hrec, key.c_str());
}

// --- Record lifecycle ---

Bcf1Ptr BcfInit() {
  Bcf1Ptr record{bcf_init()};
  if (record == nullptr) {
    throw error::Error("Failed to initialize VCF/BCF record");
  }
  return record;
}

Bcf1Ptr BcfDup(bcf1_t* const record) {
  Bcf1Ptr dup{bcf_dup(record)};
  if (dup == nullptr) {
    throw error::Error("Failed to duplicate VCF/BCF record");
  }
  return dup;
}

void BcfUnpack(bcf1_t* const record, const s32 which) {
  if (bcf_unpack(record, which) < 0) {
    throw error::Error("Failed to unpack VCF/BCF record");
  }
}

bool BcfRead(htsFile* const fp, const bcf_hdr_t* const hdr, bcf1_t* const record) {
  const auto result = bcf_read(fp, hdr, record);
  if (result == -1) {
    // EOF (normal end of file)
    return false;
  }
  if (result < -1) {
    throw error::Error("VCF/BCF parsing failed");
  }
  return true;
}

void BcfWrite(htsFile* const fp, bcf_hdr_t* const hdr, bcf1_t* const record) {
  if (bcf_write(fp, hdr, record) < 0) {
    throw error::Error("Failed to write VCF/BCF record");
  }
}

void BcfFlush(htsFile* const fp) {
  if (bcf_flush(fp) < 0) {
    throw error::Error("Failed to flush VCF/BCF file");
  }
}

bool VcfParse(kstring_t* const line, const bcf_hdr_t* const hdr, bcf1_t* const record) {
  return vcf_parse(line, hdr, record) >= 0;
}

std::string VcfFormat(const bcf_hdr_t* const hdr, const bcf1_t* const record) {
  kstring_t str = KS_INITIALIZE;
  if (vcf_format(hdr, record, &str) < 0) {
    ks_free(&str);
    throw error::Error("Failed to format VCF record");
  }
  std::string result(str.s, str.l);
  ks_free(&str);
  return result;
}

// --- INFO fields ---

bool BcfRecordHasInfo(bcf1_t* const record, const s32 int_id) {
  return bcf_get_info_id(record, int_id) != nullptr;
}

bool BcfGetInfoFlag(const bcf_hdr_t* const hdr, bcf1_t* const record, const std::string& field) {
  s32 num_entries = 0;
  // Flag retrieval never writes to *dst (see bcf_get_info_flag docs), so pass nullptr.
  return bcf_get_info_values(hdr, record, field.c_str(), nullptr, &num_entries, BCF_HT_FLAG) > 0;
}

namespace {

// Strip a single trailing NUL that htslib appends to string fields read via bcf_read.
void StripTrailingNul(std::string& value) {
  if (!value.empty() && value.back() == '\0') {
    value.pop_back();
  }
}

// Shared helper for typed INFO/FORMAT numeric retrieval.
template <typename T>
std::optional<vec<T>> GetValues(const bcf_hdr_t* const hdr,
                                bcf1_t* const record,
                                const std::string& field,
                                const s32 hts_type,
                                const bool is_info) {
  s32 num_entries = 0;
  // htslib writes a malloc'd buffer pointer into the out-parameter, so it must be the address of a
  // plain raw pointer (not a unique_ptr object). The (void**) cast mirrors htslib's own
  // bcf_get_info_int32/etc. macros. Ownership is transferred to a MallocPtr immediately after the
  // call so the buffer is freed on every return path.
  T* raw = nullptr;
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
  void** const dst = reinterpret_cast<void**>(&raw);  // NOSONAR(S3630,S5356)
  const auto result = is_info ? bcf_get_info_values(hdr, record, field.c_str(), dst, &num_entries, hts_type)
                              : bcf_get_format_values(hdr, record, field.c_str(), dst, &num_entries, hts_type);
  const MallocPtr<T> buffer{raw};
  if (result < 0) {
    return std::nullopt;
  }
  // Some htslib versions leave buffer null when zero values are returned; avoid forming
  // `nullptr + 0`, which is technically undefined behavior.
  if (num_entries <= 0 || buffer == nullptr) {
    return vec<T>{};
  }
  return vec<T>(buffer.get(), buffer.get() + num_entries);
}

// Shared helper for INFO/FORMAT string retrieval.
std::optional<std::string> GetString(const bcf_hdr_t* const hdr,
                                     bcf1_t* const record,
                                     const std::string& field,
                                     const bool is_info) {
  s32 num_entries = 0;
  // htslib writes a malloc'd buffer pointer into the out-parameter, so it must be the address of a
  // plain raw pointer (see comment in GetValues). Ownership is transferred to a MallocPtr
  // immediately after the call so the buffer is freed on every return path.
  char* raw = nullptr;
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
  void** const dst = reinterpret_cast<void**>(&raw);  // NOSONAR(S3630,S5356)
  const auto result = is_info ? bcf_get_info_values(hdr, record, field.c_str(), dst, &num_entries, BCF_HT_STR)
                              : bcf_get_format_values(hdr, record, field.c_str(), dst, &num_entries, BCF_HT_STR);
  const MallocPtr<char> buffer{raw};
  if (result < 0) {
    return std::nullopt;
  }
  if (num_entries <= 0 || buffer == nullptr) {
    return std::string{};
  }
  std::string value(buffer.get(), static_cast<size_t>(num_entries));
  StripTrailingNul(value);
  return value;
}

}  // namespace

template <>
std::optional<vec<s32>> BcfGetInfoValues(const bcf_hdr_t* const hdr, bcf1_t* const record, const std::string& field) {
  return GetValues<s32>(hdr, record, field, BCF_HT_INT, true);
}

template <>
std::optional<vec<f32>> BcfGetInfoValues(const bcf_hdr_t* const hdr, bcf1_t* const record, const std::string& field) {
  return GetValues<f32>(hdr, record, field, BCF_HT_REAL, true);
}

std::optional<std::string> BcfGetInfoString(const bcf_hdr_t* const hdr,
                                            bcf1_t* const record,
                                            const std::string& field) {
  return GetString(hdr, record, field, true);
}

void BcfUpdateInfoInt(const bcf_hdr_t* const hdr,
                      bcf1_t* const record,
                      const std::string& field,
                      const vec<s32>& values) {
  if (bcf_update_info_int32(hdr, record, field.c_str(), values.data(), static_cast<s32>(values.size())) < 0) {
    throw error::Error("Failed to update INFO field {}", field);
  }
}

void BcfUpdateInfoFloat(const bcf_hdr_t* const hdr,
                        bcf1_t* const record,
                        const std::string& field,
                        const vec<f32>& values) {
  if (bcf_update_info_float(hdr, record, field.c_str(), values.data(), static_cast<s32>(values.size())) < 0) {
    throw error::Error("Failed to update INFO field {}", field);
  }
}

void BcfUpdateInfoString(const bcf_hdr_t* const hdr,
                         bcf1_t* const record,
                         const std::string& field,
                         const std::string& value) {
  if (bcf_update_info_string(hdr, record, field.c_str(), value.c_str()) < 0) {
    throw error::Error("Failed to update INFO field {}", field);
  }
}

void BcfUpdateInfoFlag(const bcf_hdr_t* const hdr, bcf1_t* const record, const std::string& field, const bool set) {
  if (bcf_update_info_flag(hdr, record, field.c_str(), nullptr, set ? 1 : 0) < 0) {
    throw error::Error("Failed to update INFO flag {}", field);
  }
}

// --- FORMAT fields ---

bool BcfRecordHasFormat(bcf1_t* const record, const s32 int_id) {
  return bcf_get_fmt_id(record, int_id) != nullptr;
}

template <>
std::optional<vec<s32>> BcfGetFormatValues(const bcf_hdr_t* const hdr, bcf1_t* const record, const std::string& field) {
  return GetValues<s32>(hdr, record, field, BCF_HT_INT, false);
}

template <>
std::optional<vec<f32>> BcfGetFormatValues(const bcf_hdr_t* const hdr, bcf1_t* const record, const std::string& field) {
  return GetValues<f32>(hdr, record, field, BCF_HT_REAL, false);
}

std::optional<std::string> BcfGetFormatString(const bcf_hdr_t* const hdr,
                                              bcf1_t* const record,
                                              const std::string& field) {
  return GetString(hdr, record, field, false);
}

void BcfUpdateFormatInt(const bcf_hdr_t* const hdr,
                        bcf1_t* const record,
                        const std::string& field,
                        const vec<s32>& values) {
  if (bcf_update_format_int32(hdr, record, field.c_str(), values.data(), static_cast<s32>(values.size())) < 0) {
    throw error::Error("Failed to update FORMAT field {}", field);
  }
}

void BcfUpdateFormatFloat(const bcf_hdr_t* const hdr,
                          bcf1_t* const record,
                          const std::string& field,
                          const vec<f32>& values) {
  if (bcf_update_format_float(hdr, record, field.c_str(), values.data(), static_cast<s32>(values.size())) < 0) {
    throw error::Error("Failed to update FORMAT field {}", field);
  }
}

void BcfUpdateFormatString(const bcf_hdr_t* const hdr,
                           bcf1_t* const record,
                           const std::string& field,
                           const vec<std::string>& values) {
  vec<const char*> tmp;
  tmp.reserve(values.size());
  for (const auto& value : values) {
    tmp.emplace_back(value.c_str());
  }
  if (bcf_update_format_string(hdr, record, field.c_str(), tmp.data(), static_cast<s32>(tmp.size())) < 0) {
    throw error::Error("Failed to update FORMAT field {}", field);
  }
}

std::optional<vec<s32>> BcfGetGenotypes(const bcf_hdr_t* const hdr, bcf1_t* const record) {
  s32 num_entries = 0;
  // bcf_get_genotypes writes a malloc'd buffer pointer into the out-parameter, so it must be the
  // address of a plain raw pointer (not a unique_ptr object). The (void**) cast mirrors htslib's own
  // bcf_get_genotypes macro. Ownership is transferred to a MallocPtr immediately after the call so
  // the buffer is freed on every return path.
  s32* raw = nullptr;
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
  void** const dst = reinterpret_cast<void**>(&raw);  // NOSONAR(S3630,S5356)
  const auto result = bcf_get_genotypes(hdr, record, dst, &num_entries);
  const MallocPtr<s32> buffer{raw};
  if (result < 0) {
    return std::nullopt;
  }
  // Guard against `nullptr + 0` when no genotypes are present.
  if (num_entries <= 0 || buffer == nullptr) {
    return vec<s32>{};
  }
  return vec<s32>(buffer.get(), buffer.get() + num_entries);
}

void BcfUpdateGenotypes(const bcf_hdr_t* const hdr, bcf1_t* const record, const vec<s32>& genotypes) {
  if (bcf_update_genotypes(hdr, record, genotypes.data(), static_cast<s32>(genotypes.size())) < 0) {
    throw error::Error("Failed to update FORMAT field GT");
  }
}

// --- Alleles / id / filters ---

void BcfUpdateId(const bcf_hdr_t* const hdr, bcf1_t* const record, const std::string& id) {
  if (bcf_update_id(hdr, record, id.c_str()) < 0) {
    throw error::Error("Failed to update id: {}", id);
  }
}

void BcfUpdateAlleles(const bcf_hdr_t* const hdr, bcf1_t* const record, const vec<std::string>& alleles) {
  vec<const char*> tmp;
  tmp.reserve(alleles.size());
  for (const auto& value : alleles) {
    tmp.emplace_back(value.c_str());
  }
  if (bcf_update_alleles(hdr, record, tmp.data(), static_cast<s32>(tmp.size())) < 0) {
    throw error::Error("Failed to update alleles");
  }
}

void BcfUpdateFilter(const bcf_hdr_t* const hdr, bcf1_t* const record, const vec<s32>& filter_ids) {
  // bcf_update_filter takes a non-const int* (it does not modify the array, but the htslib
  // signature predates const-correctness). Copy into a mutable buffer to bridge the const ref.
  vec<s32> ids = filter_ids;
  s32* const ptr = ids.empty() ? nullptr : ids.data();
  if (bcf_update_filter(hdr, record, ptr, static_cast<s32>(ids.size())) < 0) {
    throw error::Error("Failed to update FILTER");
  }
}

void BcfAddFilter(const bcf_hdr_t* const hdr, bcf1_t* const record, const s32 filter_id) {
  if (bcf_add_filter(hdr, record, filter_id) < 0) {
    throw error::Error("Failed to add FILTER");
  }
}

s32 BcfHasFilter(const bcf_hdr_t* const hdr, bcf1_t* const record, const std::string& filter) {
  // bcf_has_filter takes a non-const char* for the filter name.
  std::string filter_copy = filter;
  return bcf_has_filter(hdr, record, filter_copy.data());
}

bool BcfIsSnp(bcf1_t* const record) {
  return bcf_is_snp(record) != 0;
}

// --- File / index / iterators ---

bool BcfCheckBgzfEof(htsFile* const fp) {
  if (fp->format.compression != bgzf) {
    return true;
  }
  BGZF* const bgzf_fp = hts_get_bgzfp(fp);
  if (bgzf_fp == nullptr) {
    return true;
  }
  return bgzf_check_EOF(bgzf_fp) != 0;
}

TbxPtr TbxIndexLoad(const fs::path& fn) {
  return TbxPtr{tbx_index_load3(fn.c_str(), nullptr, HTS_IDX_SILENT_FAIL)};
}

HtsIdxPtr BcfIndexLoad(const fs::path& fn) {
  return HtsIdxPtr{bcf_index_load3(fn.c_str(), nullptr, HTS_IDX_SILENT_FAIL)};
}

bool BcfIndexBuild(const fs::path& fn, const s32 min_shift) {
  return bcf_index_build(fn.c_str(), min_shift) == 0;
}

s32 TbxName2Id(tbx_t* const tbx, const std::string& name) {
  return tbx_name2id(tbx, name.c_str());
}

HtsItrPtr TbxItrQueryI(const tbx_t* const tbx, const s32 tid, const s64 beg, const s64 end) {
  return HtsItrPtr{tbx_itr_queryi(tbx, tid, beg, end)};
}

HtsItrPtr BcfItrQueryI(const hts_idx_t* const idx, const s32 tid, const s64 beg, const s64 end) {
  return HtsItrPtr{bcf_itr_queryi(idx, tid, beg, end)};
}

bool TbxItrNext(htsFile* const fp, tbx_t* const tbx, hts_itr_t* const itr, kstring_t* const line) {
  const auto ret = tbx_itr_next(fp, tbx, itr, line);
  if (ret < -1) {
    throw error::Error("Failed to get next record from tabix-indexed file '{}'", fp->fn);
  }
  // Return true if a record was read, false at end-of-iterator (-1).
  return ret >= 0;
}

bool BcfItrNext(htsFile* const fp, hts_itr_t* const itr, bcf1_t* const record) {
  // bcf_itr_next takes a void* record out-parameter; the bcf1_t*->void* conversion is the
  // documented htslib calling convention for this macro.
  const auto ret = bcf_itr_next(fp, itr, record);  // NOSONAR(S5356)
  if (ret < -1) {
    throw error::Error("Failed to get next record from VCF/BCF '{}'", fp->fn);
  }
  // Return true if a record was read, false at end-of-iterator (-1).
  return ret >= 0;
}

}  // namespace xoos::io
