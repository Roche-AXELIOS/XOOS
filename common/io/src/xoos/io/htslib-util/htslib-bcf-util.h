/**
 * @file htslib-bcf-util.h
 * @brief Utility functions for working with HTSlib BCF/VCF/tabix handling.
 *
 * This header provides an extremely light wrapper around HTSlib BCF/VCF functions, mirroring the
 * approach taken by htslib-util.h for SAM/BAM/CRAM. The primary intentions are to:
 *   1. Work with native C++ types like std::string, std::vector, and fs::path
 *   2. Add exceptions for error handling, due to the inconsistent error handling of HTSlib
 *   3. Be very lightweight and not attempt to abstract the HTSlib API
 *
 * Functions that wrap a call with a meaningful error code throw xoos::error::Error on failure.
 * Functions that wrap a "silent fail" lookup (e.g. loading an optional index, querying an
 * iterator that may legitimately be empty) instead return an empty/owning pointer or a boolean so
 * the caller can decide how to react, matching the underlying htslib semantics.
 */
#pragma once

#include <optional>
#include <string>
#include <type_traits>

#include <htslib/tbx.h>
#include <htslib/vcf.h>

#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "htslib-bcf-ptr.h"
#include "htslib-ptr.h"

namespace xoos::io {

// --- Header (bcf_hdr_*) ---

/// Read a BCF/VCF header from an open file. Throws on failure.
BcfHdrPtr BcfHdrRead(htsFile* fp);

/// Initialize a new BCF/VCF header in the given mode (e.g. "w"). Throws on failure.
BcfHdrPtr BcfHdrInit(const std::string& mode);

/// Duplicate a BCF/VCF header. Throws on failure.
BcfHdrPtr BcfHdrDup(const bcf_hdr_t* hdr);

/// Write a BCF/VCF header to an open file. Throws on failure.
void BcfHdrWrite(htsFile* fp, bcf_hdr_t* hdr);

/// Append a single meta-information line to a header. Throws on failure.
void BcfHdrAppend(bcf_hdr_t* hdr, const std::string& line);

/// Add a sample to a header. Throws on failure.
void BcfHdrAddSample(bcf_hdr_t* hdr, const std::string& sample);

/// Synchronize a header's internal dictionaries after modification. Throws on failure.
void BcfHdrSync(bcf_hdr_t* hdr);

/// Remove a single header line of the given type (one of BCF_HL_*) identified by `key`.
/// `key` must be a valid null-terminated C string.
void BcfHdrRemove(bcf_hdr_t* hdr, s32 type, const char* key);

/// Remove all header lines of the given type (one of BCF_HL_*).
void BcfHdrRemoveAll(bcf_hdr_t* hdr, s32 type);

/// Number of samples in the header.
s32 BcfHdrNumSamples(const bcf_hdr_t* hdr);

/// Convert a contig name to its numeric id. Returns a negative value if not found (no throw).
s32 BcfHdrName2Id(const bcf_hdr_t* hdr, const std::string& name);

/// Convert a contig numeric id to its name. Bounds-checks rid against the contig dictionary and
/// throws if it is out of range or not present.
const char* BcfHdrId2Name(const bcf_hdr_t* hdr, s32 rid);

/// Look up the integer id of a dictionary entry (e.g. BCF_DT_ID). Returns negative if absent (no throw).
s32 BcfHdrId2Int(const bcf_hdr_t* hdr, s32 dict_type, const std::string& id);

/// Look up the name of a dictionary entry by its integer id. Bounds-checks dict_type and int_id
/// against the selected dictionary and throws if either is out of range or not present.
const char* BcfHdrInt2Id(const bcf_hdr_t* hdr, s32 dict_type, s32 int_id);

/// Get a generic header record (BCF_HL_GEN) by key. Returns nullptr if not found (no throw).
bcf_hrec_t* BcfHdrGetHrec(const bcf_hdr_t* hdr, const std::string& key);

/// Find the index of a key within a header record. Returns negative if absent (no throw).
s32 BcfHrecFindKey(bcf_hrec_t* hrec, const std::string& key);

// --- Record lifecycle ---

/// Allocate a new BCF record. Throws on failure.
Bcf1Ptr BcfInit();

/// Duplicate a BCF record. Throws on failure.
Bcf1Ptr BcfDup(bcf1_t* record);

/// Unpack a BCF record up to the given level. Throws on failure.
void BcfUnpack(bcf1_t* record, s32 which);

/// Read the next BCF/VCF record. Returns true if a record was read, false on EOF.
/// Throws on a parsing error.
bool BcfRead(htsFile* fp, const bcf_hdr_t* hdr, bcf1_t* record);

/// Write a BCF/VCF record. Throws on failure.
void BcfWrite(htsFile* fp, bcf_hdr_t* hdr, bcf1_t* record);

/// Flush buffered output for a BCF/VCF file. Throws on failure.
void BcfFlush(htsFile* fp);

/// Parse a single VCF text line into a record. Returns true on success, false on parse failure.
bool VcfParse(kstring_t* line, const bcf_hdr_t* hdr, bcf1_t* record);

/// Format a BCF record as a VCF text line. Throws on failure.
std::string VcfFormat(const bcf_hdr_t* hdr, const bcf1_t* record);

// --- INFO fields ---

/// Whether the record carries the INFO field with the given integer id (no header check).
bool BcfRecordHasInfo(bcf1_t* record, s32 int_id);

/// Whether the INFO flag is set on the record.
bool BcfGetInfoFlag(const bcf_hdr_t* hdr, bcf1_t* record, const std::string& field);

/// Numeric value types supported by the typed INFO/FORMAT getters. The htslib value APIs only
/// expose 32-bit integer (BCF_HT_INT) and 32-bit float (BCF_HT_REAL) storage, so only these two
/// instantiations are defined.
template <typename T>
inline constexpr bool kIsSupportedBcfValueType = std::is_same_v<T, s32> || std::is_same_v<T, f32>;

/// Retrieve typed INFO values. Returns nullopt if the field is absent or on error.
/// Only `s32` and `f32` are supported; any other type fails to compile.
// The primary template is intentionally never instantiated: it exists only to emit a clear
// static_assert for unsupported types. Only the s32/f32 specializations below are used.
template <typename T>
std::optional<vec<T>> BcfGetInfoValues([[maybe_unused]] const bcf_hdr_t* hdr,
                                       [[maybe_unused]] bcf1_t* record,
                                       [[maybe_unused]] const std::string& field) {  // NOSONAR(S6009)
  static_assert(kIsSupportedBcfValueType<T>, "BcfGetInfoValues supports only s32 (BCF_HT_INT) and f32 (BCF_HT_REAL)");
  return std::nullopt;
}

template <>
std::optional<vec<s32>> BcfGetInfoValues(const bcf_hdr_t* hdr, bcf1_t* record, const std::string& field);
template <>
std::optional<vec<f32>> BcfGetInfoValues(const bcf_hdr_t* hdr, bcf1_t* record, const std::string& field);

/// Retrieve a string INFO value. Returns nullopt if the field is absent or on error.
std::optional<std::string> BcfGetInfoString(const bcf_hdr_t* hdr, bcf1_t* record, const std::string& field);

void BcfUpdateInfoInt(const bcf_hdr_t* hdr, bcf1_t* record, const std::string& field, const vec<s32>& values);
void BcfUpdateInfoFloat(const bcf_hdr_t* hdr, bcf1_t* record, const std::string& field, const vec<f32>& values);
void BcfUpdateInfoString(const bcf_hdr_t* hdr, bcf1_t* record, const std::string& field, const std::string& value);
void BcfUpdateInfoFlag(const bcf_hdr_t* hdr, bcf1_t* record, const std::string& field, bool set);

// --- FORMAT fields ---

/// Whether the record carries the FORMAT field with the given integer id (no header check).
bool BcfRecordHasFormat(bcf1_t* record, s32 int_id);

/// Retrieve typed FORMAT values. Returns nullopt if the field is absent or on error.
/// Only `s32` and `f32` are supported; any other type fails to compile.
// See BcfGetInfoValues: this primary template exists only for the static_assert and is never used.
template <typename T>
std::optional<vec<T>> BcfGetFormatValues([[maybe_unused]] const bcf_hdr_t* hdr,
                                         [[maybe_unused]] bcf1_t* record,
                                         [[maybe_unused]] const std::string& field) {  // NOSONAR(S6009)
  static_assert(kIsSupportedBcfValueType<T>, "BcfGetFormatValues supports only s32 (BCF_HT_INT) and f32 (BCF_HT_REAL)");
  return std::nullopt;
}

template <>
std::optional<vec<s32>> BcfGetFormatValues(const bcf_hdr_t* hdr, bcf1_t* record, const std::string& field);
template <>
std::optional<vec<f32>> BcfGetFormatValues(const bcf_hdr_t* hdr, bcf1_t* record, const std::string& field);

/// Retrieve a string FORMAT value. Returns nullopt if the field is absent or on error.
std::optional<std::string> BcfGetFormatString(const bcf_hdr_t* hdr, bcf1_t* record, const std::string& field);

void BcfUpdateFormatInt(const bcf_hdr_t* hdr, bcf1_t* record, const std::string& field, const vec<s32>& values);
void BcfUpdateFormatFloat(const bcf_hdr_t* hdr, bcf1_t* record, const std::string& field, const vec<f32>& values);
void BcfUpdateFormatString(const bcf_hdr_t* hdr,
                           bcf1_t* record,
                           const std::string& field,
                           const vec<std::string>& values);

/// Retrieve the raw genotype integers. Returns nullopt on error (negative return from htslib).
std::optional<vec<s32>> BcfGetGenotypes(const bcf_hdr_t* hdr, bcf1_t* record);

/// Update the genotype (GT) field from a flat genotype buffer. Throws on failure.
void BcfUpdateGenotypes(const bcf_hdr_t* hdr, bcf1_t* record, const vec<s32>& genotypes);

// --- Alleles / id / filters ---

void BcfUpdateId(const bcf_hdr_t* hdr, bcf1_t* record, const std::string& id);
void BcfUpdateAlleles(const bcf_hdr_t* hdr, bcf1_t* record, const vec<std::string>& alleles);

/// Replace the record's filters with the given filter ids (empty clears all). Throws on failure.
void BcfUpdateFilter(const bcf_hdr_t* hdr, bcf1_t* record, const vec<s32>& filter_ids);

/// Append a single filter id to the record. Throws on failure.
void BcfAddFilter(const bcf_hdr_t* hdr, bcf1_t* record, s32 filter_id);

/// Whether the record has the named filter. Returns the raw htslib code (1 = yes, 0 = no, -1 = undefined).
s32 BcfHasFilter(const bcf_hdr_t* hdr, bcf1_t* record, const std::string& filter);

/// Whether the record is a SNP.
bool BcfIsSnp(bcf1_t* record);

// --- File / index / iterators ---

/// Check that a BGZF-compressed file has a valid EOF marker. Returns true when valid or when the
/// file is not BGZF-compressed; returns false when the EOF marker is missing (file may be truncated).
bool BcfCheckBgzfEof(htsFile* fp);

/// Load a tabix (.tbi/.csi) index. Returns an empty pointer if no index is available (silent fail).
TbxPtr TbxIndexLoad(const fs::path& fn);

/// Load a BCF (.csi) index. Returns an empty pointer if no index is available (silent fail).
HtsIdxPtr BcfIndexLoad(const fs::path& fn);

/// Build a BCF/VCF index. Returns true on success, false on failure (does not throw).
bool BcfIndexBuild(const fs::path& fn, s32 min_shift);

/// Convert a contig name to its tabix id. Returns negative if not found (no throw).
s32 TbxName2Id(tbx_t* tbx, const std::string& name);

/// Query a tabix iterator for a region. Returns an empty pointer if the region cannot be queried.
HtsItrPtr TbxItrQueryI(const tbx_t* tbx, s32 tid, s64 beg, s64 end);

/// Query a BCF iterator for a region. Returns an empty pointer if the region cannot be queried.
HtsItrPtr BcfItrQueryI(const hts_idx_t* idx, s32 tid, s64 beg, s64 end);

/// Advance a tabix iterator into the given kstring. Returns true if a line was read, false at end
/// of iterator. Throws on a read error.
bool TbxItrNext(htsFile* fp, tbx_t* tbx, hts_itr_t* itr, kstring_t* line);

/// Advance a BCF iterator into the given record. Returns true if a record was read, false at end of
/// iterator. Throws on a read error.
bool BcfItrNext(htsFile* fp, hts_itr_t* itr, bcf1_t* record);

}  // namespace xoos::io
