/**
 * @file htslib-util.h
 * @brief Utility functions for working with HTSlib, including SAM/BAM/CRAM file handling.
 *
 * This header provides an extremely light wrapper around HTSlib functions, the primary intentions are to:
 *   1. Work with native C++ types like std::string and std::vector
 *   2. Add exceptions for error handling, due to the inconsistent error handling of HTSlib
 *   3. Be very lightweight and not attempt to abstract the HTSlib API
 */
#pragma once

#include <filesystem>
#include <optional>

#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "htslib-ptr.h"
#include "xoos/types/vec.h"

namespace xoos::io {

HtsIdxPtr HtsIdxLoad(const fs::path& fn, int fmt);

HtsIdxPtr SamIndexLoad(htsFile* fp, const fs::path& fn);

void HtsFlush(htsFile* fp);

HtsFilePtr HtsOpenFormat(const fs::path& fn, const std::string& mode, const htsFormat* fmt);

HtsItrPtr SamItrRegion(const hts_idx_t* idx, sam_hdr_t* hdr, hts_reglist_t* reglist, u32 reg_count);

int SamItrMultiNext(htsFile* fp, hts_itr_t* itr, bam1_t* record);

HtsFilePtr HtsOpen(const fs::path& fn, const std::string& mode);

SamHdrPtr SamHdrInit();

// this returns a null-terminated char pointer to the name of the target with the given tid so we leave
// it as a char*
const char* SamHdrTid2Name(const sam_hdr_t* hdr, s32 tid);

s64 SamHdrTid2Length(const sam_hdr_t* hdr, s32 tid);

s32 SamHdrNRef(const sam_hdr_t* hdr);

SamFilePtr SamOpen(const fs::path& fn, const std::string& mode);

void SamWrite1(samFile* fp, const sam_hdr_t* hdr, const bam1_t* record);

SamHdrPtr SamHdrRead(samFile* fp);

void SamHdrWrite(samFile* fp, const sam_hdr_t* hdr);

vec<u32> SamParseCigar(const std::string& in);

struct SqHdrLine {
  std::string sn{};
  u32 ln{};
};

struct PgHdrLine {
  // ID is required
  std::string id{};
  // Program name (PN) is required
  std::string program_name{};
  // Version number (VN) is optional but recommended
  std::optional<std::string> version_number{};
  // Command line (CL) is optional but recommended
  std::optional<std::string> command_line{};
  // Previous program ID (PP) is optional
  std::optional<std::string> previous_program_id{};
};

void SamHdrAddSqLine(sam_hdr_t* hdr, const SqHdrLine& line);

// Add a PG line to the SAM header
// ID and PN are required, VN and CL are optional but recommended, PP is optional
void SamHdrAddPgLine(sam_hdr_t* hdr, const PgHdrLine& line);

void SamIdxSave(htsFile* fp);

void SamIdxInit(htsFile* fp, sam_hdr_t* hdr, int min_shift, const std::string& idx_fn);

bool SamItrNext(htsFile* fp, hts_itr_t* itr, bam1_t* record);

HtsItrPtr SamItrQueryI(hts_idx_t* idx, s32 tid, s64 beg, s64 end);

HtsItrPtr SamItrQueryS(const hts_idx_t* idx, sam_hdr_t* hdr, const std::string& region);

HtsItrPtr SamItrQuerySNoThrow(const hts_idx_t* idx, sam_hdr_t* hdr, const std::string& region);

s32 SamHdrName2Tid(sam_hdr_t* hdr, const std::string& target);

void BamAuxAppend(bam1_t* record, const std::string& name, u32 value);

void BamAuxAppend(bam1_t* record, const std::string& name, const std::string& value);

void BamAuxUpdateStr(bam1_t* record, const std::string_view& tag, const std::string_view& value);

template <typename T>
std::optional<T> BamAuxGet(const bam1_t* record, const std::string& tag);

void BamSet1(bam1_t* bam,
             const std::string_view& qname,
             u16 flag,
             s32 tid,
             s64 pos,
             u8 mapq,
             const vec<u32>& cigar,
             s32 mtid,
             s64 mpos,
             s32 isize,
             const std::string_view& seq,
             const std::string_view& qual,
             size_t l_aux);

}  // namespace xoos::io
