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

#include <optional>

#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "htslib-ptr.h"
#include "xoos/types/vec.h"

namespace xoos::io {

HtsIdxPtr HtsIdxLoad(const fs::path& fn, s32 fmt);

HtsIdxPtr SamIndexLoad(htsFile* fp, const fs::path& fn);

// Report whether an HTS file is BGZF-compressed and carries an intact EOF marker.
// Returns false if the file is not BGZF-compressed or the EOF marker is missing
// (a sign the file may be truncated). Callers own the user-facing diagnostic, so
// this is a predicate rather than a throwing check. Must be called before loading
// an index.
bool HtsHasBgzfEof(htsFile* fp);

// Read the next alignment record (sam_read1). Returns true if a record was read,
// false at end of file. Throws on read/parse error.
bool SamRead1(samFile* fp, sam_hdr_t* hdr, bam1_t* record);

// Load an existing FASTA index (fai_load3 with index creation disabled).
// Throws if the index is missing or cannot be opened.
FaIdxPtr FaiLoad(const fs::path& fn);

void HtsFlush(htsFile* fp);

HtsFilePtr HtsOpenFormat(const fs::path& fn, const std::string& mode, const htsFormat* fmt);

HtsItrPtr SamItrRegion(const hts_idx_t* idx, sam_hdr_t* hdr, hts_reglist_t* reglist, u32 reg_count);

s32 SamItrMultiNext(htsFile* fp, hts_itr_t* itr, bam1_t* record);

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
  std::string sn;
  u32 ln{};
};

struct PgHdrLine {
  // ID is required
  std::string id;
  // Program name (PN) is required
  std::string program_name;
  // Version number (VN) is optional but recommended
  std::optional<std::string> version_number;
  // Command line (CL) is optional but recommended
  std::optional<std::string> command_line;
  // Previous program ID (PP) is optional
  std::optional<std::string> previous_program_id;
};

constexpr auto kPgTag = "PG";
constexpr auto kPgProgramIdTag = "ID";
constexpr auto kPgProgramNameTag = "PN";
constexpr auto kPgVersionNumberTag = "VN";
constexpr auto kPgCommandLineTag = "CL";
constexpr auto kPgPreviousProgramIdTag = "PP";

void SamHdrAddSqLine(sam_hdr_t* hdr, const SqHdrLine& line);

// Add a PG line to the SAM header
// ID and PN are required, VN and CL are optional but recommended, PP is optional.
// If a PG line with the same ID already exists, a unique suffix will be appended to the new ID to avoid collision
// If PP is not specified, the last PG line (if any) will be used as the previous program.
void SamHdrAddPgLine(sam_hdr_t* hdr, const PgHdrLine& line);

void SamIdxSave(htsFile* fp);

void SamIdxInit(htsFile* fp, sam_hdr_t* hdr, s32 min_shift, const std::string& idx_fn);

bool SamItrNext(htsFile* fp, hts_itr_t* itr, bam1_t* record);

/// Query a BAM index by numeric tid and 0-based half-open coordinates [beg, end).
/// Wraps sam_itr_queryi() which uses the same coordinate convention.
HtsItrPtr SamItrQueryI(const hts_idx_t* idx, s32 tid, s64 beg, s64 end);

/// Query a BAM index using a samtools-style region string "chr:start-end" (1-based, closed).
/// Wraps sam_itr_querys(). Callers must format the string in 1-based closed coordinates;
/// formatting a xoos::genome::Region produces that form from 0-based half-open coordinates.
HtsItrPtr SamItrQueryS(const hts_idx_t* idx, sam_hdr_t* hdr, const std::string& region);

/// Same as SamItrQueryS but returns nullptr instead of throwing on parse failure.
HtsItrPtr SamItrQuerySNoThrow(const hts_idx_t* idx, sam_hdr_t* hdr, const std::string& region);

s32 SamHdrName2Tid(sam_hdr_t* hdr, const std::string& target);

void BamAuxAppend(bam1_t* record, const std::string& name, u32 value);

void BamAuxAppend(bam1_t* record, const std::string& name, s32 value);

void BamAuxAppend(bam1_t* record, const std::string& name, const std::string& value);

// Append an aux tag of an arbitrary htslib type code. Used when copying tags
// verbatim from another record, where the type and raw bytes are already known.
void BamAuxAppend(bam1_t* record, const std::string& name, char type, u32 len, const u8* data);

void BamAuxUpdateStr(bam1_t* record, const std::string_view& tag, const std::string_view& value);

// Return the rightmost reference coordinate (exclusive) covered by the record's
// alignment. Pure computation; htslib reports no error state here, so unlike the
// fallible wrappers above this does not throw.
s64 BamEndPos(const bam1_t* record);

// Compute the htslib bin number for a [beg, end) reference interval.
// Mirrors hts_reg2bin: `min_shift` is the bits for the smallest bin and
// `n_lvls` the number of binning levels (BAM defaults are 14 and 5).
// Pure computation; htslib reports no error state here, so this does not throw.
u16 HtsReg2Bin(s64 beg, s64 end, s32 min_shift, s32 n_lvls);

template <typename T>
std::optional<T> BamAuxGet(const bam1_t* record, const std::string& tag);

/// Wrapper around bam_set1. @p pos and @p mpos are 0-based mapping positions.
void BamSet1(bam1_t* bam_record,
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

Bam1Ptr BamInit1();

Bam1Ptr BamDup1(const bam1_t* src);

HtsItrMultiPtr SamItrMultiQueryS(const hts_idx_t* idx, sam_hdr_t* hdr, const std::string& region);

HtsItrMultiPtr SamItrMultiQueryI(const hts_idx_t* idx, s32 tid, s64 beg, s64 end);

void SamIndexBuild3(const fs::path& bam_filename, const fs::path& index_filename, s32 min_shift, s32 threads);

}  // namespace xoos::io
