#include "htslib-util.h"

#include <array>
#include <cerrno>
#include <cstddef>
#include <filesystem>
#include <gsl/gsl>
#include <optional>
#include <string>
#include <string_view>
#include <system_error>

#include <fmt/format.h>
#include <htslib/bgzf.h>
#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/tbx.h>

#include <xoos/error/error.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/util/adopt.h>

#include "xoos/io/htslib-util/htslib-ptr.h"
#include "xoos/io/malloc-ptr.h"

namespace xoos::io {

namespace {

/// htslib sam_hdr_name2tid return value indicating that the BAM header could not be parsed.
constexpr s32 kSamHdrName2TidParseError = -2;
constexpr size_t kPgProgramNameTagIndex = 0;
constexpr size_t kPgProgramNameValueIndex = 1;
constexpr size_t kPgFirstOptionalTagIndex = 2;
constexpr size_t kPgFirstOptionalValueIndex = 3;
constexpr size_t kPgSecondOptionalTagIndex = 4;
constexpr size_t kPgSecondOptionalValueIndex = 5;
constexpr size_t kPgThirdOptionalTagIndex = 6;
constexpr size_t kPgThirdOptionalValueIndex = 7;
constexpr size_t kPgFieldCountOneTag = kPgProgramNameValueIndex + 1;
constexpr size_t kPgFieldCountTwoTags = kPgFirstOptionalValueIndex + 1;
constexpr size_t kPgFieldCountThreeTags = kPgSecondOptionalValueIndex + 1;
constexpr size_t kPgFieldCountFourTags = kPgThirdOptionalValueIndex + 1;
constexpr size_t kPgMaxFieldCount = kPgFieldCountFourTags;

/// Returns the system error string for the current errno value (thread-safe).
std::string ErrnoMessage() {
  return std::generic_category().message(errno);  // NOSONAR(S985) — S985 flags any errno read; here errno is captured
                                                  // immediately by std::generic_category().message() before any other
                                                  // syscall can change it, which is the standard thread-safe idiom
}

/// Renders a region list to a human-readable string for error messages.
std::string ToString(const hts_reglist_t* const reglist) {
  if (reglist != nullptr && reglist->count > 0) {
    return reglist->reg;
  }
  return "No regions";
}

}  // namespace

HtsIdxPtr HtsIdxLoad(const fs::path& fn, const s32 fmt) {
  HtsIdxPtr idx{hts_idx_load(fn.c_str(), fmt)};
  if (idx == nullptr) {
    throw error::Error("Failed to read index for HTS file: '{}'", fn);
  }
  return idx;
}

HtsIdxPtr SamIndexLoad(htsFile* const fp, const fs::path& fn) {
  HtsIdxPtr idx{sam_index_load(fp, fn.c_str())};
  if (idx == nullptr) {
    throw error::Error("Failed to read index for BAM/CRAM file: '{}'", fn);
  }
  return idx;
}

bool HtsHasBgzfEof(htsFile* const fp) {
  BGZF* const bgzf_file = hts_get_bgzfp(fp);
  if (bgzf_file == nullptr) {
    return false;
  }
  // M23_103: extract call result before comparison to avoid side effect in RHS of != operator.
  const auto eof_status = bgzf_check_EOF(bgzf_file);
  return eof_status != 0;
}

bool SamRead1(samFile* const fp, sam_hdr_t* const hdr, bam1_t* const record) {
  const auto ret = sam_read1(fp, hdr, record);
  if (ret < -1) {
    throw error::Error("Error reading BAM record from '{}'", fp->fn);
  }
  // ret == -1 signals end of file; anything >= 0 is the record length.
  return ret >= 0;
}

FaIdxPtr FaiLoad(const fs::path& fn) {
  // Pass index creation flag 0 so a missing index is reported rather than built.
  FaIdxPtr fai{fai_load3(fn.c_str(), nullptr, nullptr, 0)};
  if (fai == nullptr) {
    throw error::Error("Failed to open index for FASTA file: '{}'", fn);
  }
  return fai;
}

void HtsFlush(htsFile* const fp) {
  if (hts_flush(fp) < 0) {
    throw error::Error("Failed to flush HTS file: '{}'", fp->fn);
  }
}

HtsFilePtr HtsOpenFormat(const fs::path& fn, const std::string& mode, const htsFormat* const fmt) {
  HtsFilePtr fp{hts_open_format(fn.c_str(), mode.c_str(), fmt)};
  if (fp == nullptr) {
    throw error::Error("Failed to open HTS file '{}': {}", fn, ErrnoMessage());
  }
  return fp;
}

HtsItrPtr SamItrRegion(const hts_idx_t* const idx,
                       sam_hdr_t* const hdr,
                       hts_reglist_t* const reglist,
                       const u32 reg_count) {
  HtsItrPtr itr_ptr(sam_itr_regions(idx, hdr, reglist, reg_count));
  // Check if the iterator was created successfully if not, throw an error with the appropriate message
  if (itr_ptr == nullptr) {
    throw error::Error(
        "Failed to create iterator for regions in BAM file with reglist '{}': {}", ToString(reglist), ErrnoMessage());
  }
  return itr_ptr;
}

// Reads the next record from a multi-region iterator. Returns the record length (>= 0) on success,
// -1 at end of iterator, or throws on error (< -1).
s32 SamItrMultiNext(htsFile* const fp, hts_itr_t* const itr, bam1_t* const record) {
  const s32 return_value = sam_itr_multi_next(fp, itr, record);
  if (return_value < -1) {
    throw error::Error("Error reading alignment '{}': {}", fp->fn, ErrnoMessage());
  }
  return return_value;
}

HtsFilePtr HtsOpen(const fs::path& fn, const std::string& mode) {
  HtsFilePtr fp{hts_open(fn.c_str(), mode.c_str())};
  if (fp == nullptr) {
    throw error::Error("Failed to open HTS file '{}': {}", fn, ErrnoMessage());
  }
  return fp;
}

SamHdrPtr SamHdrInit() {
  SamHdrPtr hdr{sam_hdr_init()};
  if (hdr == nullptr) {
    throw error::Error("Failed to initialize SAM header");
  }
  return hdr;
}

// returning a null-terminated char pointer to the name of the target with the given tid
const char* SamHdrTid2Name(const sam_hdr_t* const hdr, const s32 tid) {
  const auto* const name = sam_hdr_tid2name(hdr, tid);
  if (name == nullptr) {
    throw error::Error("Failed to get name for tid {} from BAM header", tid);
  }
  return name;
}

s64 SamHdrTid2Length(const sam_hdr_t* const hdr, const s32 tid) {
  const auto len = sam_hdr_tid2len(hdr, tid);
  if (len <= 0) {
    throw error::Error("Failed to get length for tid {} from BAM header", tid);
  }
  return len;
}

s32 SamHdrNRef(const sam_hdr_t* const hdr) {
  const auto nref = sam_hdr_nref(hdr);
  if (nref < 0) {
    throw error::Error("Failed to get number of references from BAM header");
  }
  return nref;
}

SamFilePtr SamOpen(const fs::path& fn, const std::string& mode) {
  SamFilePtr in{sam_open(fn.c_str(), mode.c_str())};
  if (in == nullptr) {
    throw error::Error("Failed to open BAM file: {}", fn);
  }
  return in;
}

void SamWrite1(samFile* const fp, const sam_hdr_t* const hdr, const bam1_t* const record) {
  if (sam_write1(fp, hdr, record) < 0) {
    throw error::Error("Failed to write BAM record '{}' to '{}': {}", bam_get_qname(record), fp->fn, ErrnoMessage());
  }
}

SamHdrPtr SamHdrRead(samFile* const fp) {
  SamHdrPtr hdr{sam_hdr_read(fp)};
  if (hdr == nullptr) {
    throw error::Error("Failed to read header from BAM file: {}", fp->fn);
  }
  return hdr;
}

void SamHdrWrite(samFile* const fp, const sam_hdr_t* const hdr) {
  if (sam_hdr_write(fp, hdr) == -1) {
    throw error::Error("Failed to write BAM header to '{}'", fp->fn);
  }
}

vec<u32> SamParseCigar(const std::string& in) {
  u32* cigar_buf = nullptr;
  size_t cigar_len = 0;
  if (sam_parse_cigar(in.c_str(), nullptr, &cigar_buf, &cigar_len) == -1) {
    throw error::Error("Failed to parse CIGAR string '{}'", in);
  }
  const auto cigar = MallocPtr<u32>(cigar_buf);
  return vec<u32>{cigar.get(), cigar.get() + cigar_len};
}

void SamHdrAddSqLine(sam_hdr_t* const hdr, const SqHdrLine& line) {
  const auto ln = std::to_string(line.ln);
  // sam_hdr_add_line is a C variadic API (NULL-terminated key/value list).
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-vararg,hicpp-vararg)
  if (sam_hdr_add_line(hdr, "SQ", "SN", line.sn.c_str(), "LN", ln.c_str(), nullptr) == -1) {
    throw error::Error("Failed to add SQ header line for '{}'", line.sn);
  }
}

void SamHdrAddPgLine(sam_hdr_t* const hdr, const PgHdrLine& line) {
  if (line.id.empty()) {
    throw error::Error("Cannot add PG header line with empty ID");
  }
  if (line.program_name.empty()) {
    throw error::Error("Cannot add PG header line with empty program name");
  }

  std::array<std::optional<std::string>, kPgMaxFieldCount> tags = {kPgProgramNameTag,
                                                                   line.program_name,
                                                                   std::nullopt,
                                                                   std::nullopt,
                                                                   std::nullopt,
                                                                   std::nullopt,
                                                                   std::nullopt,
                                                                   std::nullopt};
  size_t pg_field_count = kPgFieldCountOneTag;
  // Need to silence SonarQube S6009 (Replace this const reference to "std::string" by a "std::string_view")
  // std::string is used here intentionally because sam_hdr_add_pg expects C-style strings (const char*) and
  // std::string_view cannot be converted to const char* directly.
  const auto add_optional_tag = [&tags, &pg_field_count](const std::string& key,  // NOSONAR cpp:S6009
                                                         const std::optional<std::string>& value) {
    if (key != "" && value.has_value() && pg_field_count + 1 < tags.size()) {
      tags.at(pg_field_count) = key;
      ++pg_field_count;
      tags.at(pg_field_count) = *value;
      ++pg_field_count;
    }
  };
  add_optional_tag(kPgVersionNumberTag, line.version_number);
  add_optional_tag(kPgCommandLineTag, line.command_line);
  add_optional_tag(kPgPreviousProgramIdTag, line.previous_program_id);

  s32 result{};
  // sam_hdr_add_pg handles @PG ID collisions and links new records to existing program chains.
  // NOLINTBEGIN(cppcoreguidelines-pro-type-vararg,hicpp-vararg) - sam_hdr_add_pg is a C variadic htslib API
  switch (pg_field_count) {
    case kPgFieldCountOneTag:  // NOSONAR cpp:S1151 - high line count due to variadic arguments
      result = sam_hdr_add_pg(hdr,
                              line.id.c_str(),
                              tags[kPgProgramNameTagIndex].value().c_str(),
                              tags[kPgProgramNameValueIndex].value().c_str(),
                              nullptr);
      break;
    case kPgFieldCountTwoTags:  // NOSONAR cpp:S1151 - high line count due to variadic arguments
      result = sam_hdr_add_pg(hdr,
                              line.id.c_str(),
                              tags[kPgProgramNameTagIndex].value().c_str(),
                              tags[kPgProgramNameValueIndex].value().c_str(),
                              tags[kPgFirstOptionalTagIndex].value().c_str(),
                              tags[kPgFirstOptionalValueIndex].value().c_str(),
                              nullptr);
      break;
    case kPgFieldCountThreeTags:  // NOSONAR cpp:S1151 - high line count due to variadic arguments
      result = sam_hdr_add_pg(hdr,
                              line.id.c_str(),
                              tags[kPgProgramNameTagIndex].value().c_str(),
                              tags[kPgProgramNameValueIndex].value().c_str(),
                              tags[kPgFirstOptionalTagIndex].value().c_str(),
                              tags[kPgFirstOptionalValueIndex].value().c_str(),
                              tags[kPgSecondOptionalTagIndex].value().c_str(),
                              tags[kPgSecondOptionalValueIndex].value().c_str(),
                              nullptr);
      break;
    case kPgFieldCountFourTags:  // NOSONAR cpp:S1151 - high line count due to variadic arguments
      result = sam_hdr_add_pg(hdr,
                              line.id.c_str(),
                              tags[kPgProgramNameTagIndex].value().c_str(),
                              tags[kPgProgramNameValueIndex].value().c_str(),
                              tags[kPgFirstOptionalTagIndex].value().c_str(),
                              tags[kPgFirstOptionalValueIndex].value().c_str(),
                              tags[kPgSecondOptionalTagIndex].value().c_str(),
                              tags[kPgSecondOptionalValueIndex].value().c_str(),
                              tags[kPgThirdOptionalTagIndex].value().c_str(),
                              tags[kPgThirdOptionalValueIndex].value().c_str(),
                              nullptr);
      break;
    default:
      throw error::Error("Invalid number of @PG header tags");
  }
  // NOLINTEND(cppcoreguidelines-pro-type-vararg,hicpp-vararg)

  if (result != 0) {
    throw error::Error("Failed to add PG header line for ID '{}'", line.id);
  }
}

void SamIdxSave(htsFile* const fp) {
  if (sam_idx_save(fp) < 0) {
    throw error::Error("Failed to save BAM index for '{}'", fp->fn);
  }
}

void SamIdxInit(htsFile* const fp, sam_hdr_t* const hdr, const s32 min_shift, const std::string& idx_fn) {
  if (sam_idx_init(fp, hdr, min_shift, idx_fn.c_str()) < 0) {
    throw error::Error("Failed to initialize BAM index for '{}'", fp->fn);
  }
}

bool SamItrNext(htsFile* const fp, hts_itr_t* const itr, bam1_t* const record) {
  const auto ret = sam_itr_next(fp, itr, record);
  if (ret < -1) {
    throw error::Error("Failed to get next record from BAM '{}': {}", fp->fn, ErrnoMessage());
  }
  // Return true if a record was read, false if no more records are available.
  return ret >= 0;
}

HtsItrPtr SamItrQueryI(const hts_idx_t* const idx, s32 tid, s64 beg, s64 end) {
  auto itr = HtsItrPtr{sam_itr_queryi(idx, tid, beg, end)};
  if (itr == nullptr) {
    throw error::Error("Failed to create iterator {}-{}:{}: {}", tid, beg, end, ErrnoMessage());
  }
  return itr;
}

HtsItrPtr SamItrQueryS(const hts_idx_t* const idx, sam_hdr_t* const hdr, const std::string& region) {
  auto itr = HtsItrPtr{sam_itr_querys(idx, hdr, region.c_str())};
  if (itr == nullptr) {
    throw error::Error("Failed to create iterator {} : {}", region, ErrnoMessage());
  }
  return itr;
}

HtsItrPtr SamItrQuerySNoThrow(const hts_idx_t* const idx, sam_hdr_t* const hdr, const std::string& region) {
  return HtsItrPtr{sam_itr_querys(idx, hdr, region.c_str())};
}

s32 SamHdrName2Tid(sam_hdr_t* const hdr, const std::string& target) {
  const auto tid = sam_hdr_name2tid(hdr, target.c_str());
  if (tid == -1) {
    throw error::Error("Failed to find target '{}' in BAM header", target);
  }
  if (tid == kSamHdrName2TidParseError) {
    throw error::Error("Failed to parse BAM header");
  }
  return tid;
}

void BamAuxAppend(bam1_t* const record, const std::string& name, const char type, const u32 len, const u8* const data) {
  // htslib reads exactly two tag characters from `name`; a shorter name would
  // make bam_aux_append read past the end of the buffer.
  constexpr size_t kBamAuxTagNameLength = 2;
  if (name.size() != kBamAuxTagNameLength) {
    throw error::Error("BAM aux tag name must be exactly 2 characters, got '{}'", name);
  }
  // gsl::narrow throws if `len` does not fit in s32 rather than silently wrapping.
  if (bam_aux_append(record, name.c_str(), type, gsl::narrow<s32>(len), data) == -1) {
    throw error::Error("Failed to append BAM aux tag '{}' of type '{}': {}", name, type, ErrnoMessage());
  }
}

void BamAuxUpdateStr(bam1_t* const record, const std::string_view& tag, const std::string_view& value) {
  // Length is passed explicitly to bam_aux_update_str, so the non-null-terminated
  // data() pointer is safe here.
  // NOLINTNEXTLINE(bugprone-suspicious-stringview-data-usage)
  if (bam_aux_update_str(record, tag.data(), static_cast<s32>(value.length() + 1), value.data()) != 0) {
    throw error::Error("Failed to update BAM aux tag '{}'", tag);
  }
}

s64 BamEndPos(const bam1_t* const record) {
  return bam_endpos(record);
}

u16 HtsReg2Bin(const s64 beg, const s64 end, const s32 min_shift, const s32 n_lvls) {
  return static_cast<u16>(hts_reg2bin(beg, end, min_shift, n_lvls));
}

// S5356/S6022: bam_aux_append requires a const uint8_t* byte buffer. Casting to
// unsigned char* (u8 = uint8_t = unsigned char) is explicitly permitted by the C++
// aliasing rules ([basic.lval]): any object may be read through an unsigned char
// pointer. The sizeof() length argument ensures only the correct number of bytes are read.
void BamAuxAppend(bam1_t* const record, const std::string& name, const u32 value) {  // NOSONAR(S5356,S6022)
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
  if (bam_aux_append(record, name.c_str(), 'I', sizeof(u32), reinterpret_cast<const u8*>(&value)) == -1) {
    throw error::Error("Failed to append BAM aux tag '{}' with value '{}': {}", name, value, ErrnoMessage());
  }
}

void BamAuxAppend(bam1_t* const record, const std::string& name, const s32 value) {  // NOSONAR(S5356,S6022)
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
  if (bam_aux_append(record, name.c_str(), 'i', sizeof(s32), reinterpret_cast<const u8*>(&value)) == -1) {
    throw error::Error("Failed to append BAM aux tag '{}' with value '{}': {}", name, value, ErrnoMessage());
  }
}

void BamAuxAppend(bam1_t* const record, const std::string& name, const std::string& value) {  // NOSONAR(S5356)
  const auto len = static_cast<s32>(value.size()) + 1;
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
  if (bam_aux_append(record, name.c_str(), 'Z', len, reinterpret_cast<const u8*>(value.c_str())) == -1) {
    throw error::Error("Failed to append BAM aux tag '{}' with value '{}': {}", name, value, ErrnoMessage());
  }
}

// S5536: explicit template specialization; must live in one .cpp to avoid multiple
// definitions — called by downstream modules that link this library.
template <>  // NOSONAR(S5536)
std::optional<u32> BamAuxGet(const bam1_t* record, const std::string& tag) {
  const auto* const tag_ptr = bam_aux_get(record, tag.c_str());
  if (tag_ptr == nullptr) {
    return std::nullopt;
  }
  return bam_aux2i(tag_ptr);
}

// S5536: explicit template specialization; must live in one .cpp to avoid multiple
// definitions — called by downstream modules that link this library.
template <>  // NOSONAR(S5536)
std::optional<s64> BamAuxGet(const bam1_t* record, const std::string& tag) {
  const auto* const tag_ptr = bam_aux_get(record, tag.c_str());
  if (tag_ptr == nullptr) {
    return std::nullopt;
  }
  return bam_aux2i(tag_ptr);
}

// S5536: explicit template specialization; must live in one .cpp to avoid multiple
// definitions — called by downstream modules that link this library.
template <>  // NOSONAR(S5536)
std::optional<std::string> BamAuxGet(const bam1_t* record, const std::string& tag) {
  const auto* const tag_ptr = bam_aux_get(record, tag.c_str());
  if (tag_ptr == nullptr) {
    return std::nullopt;
  }
  const auto* str = bam_aux2Z(tag_ptr);
  if (str == nullptr) {
    return std::nullopt;
  }
  return std::string(str);
}

void BamSet1(bam1_t* const bam_record,  // NOSONAR(S107) — S107 fires because bam_set1() requires all 13 positional
                                        // parameters; introducing a parameter struct would add complexity without
                                        // benefit since BamSet1 is a thin wrapper over the C API
             const std::string_view& qname,
             const u16 flag,
             const s32 tid,
             const s64 pos,
             const u8 mapq,
             const vec<u32>& cigar,
             const s32 mtid,
             const s64 mpos,
             const s32 isize,
             const std::string_view& seq,
             const std::string_view& qual,
             const size_t l_aux) {
  // bam_set1 receives explicit lengths (qname.length(), seq.length()) alongside the
  // string_view data() pointers, so the non-null-terminated buffers are safe here.
  // NOLINTBEGIN(bugprone-suspicious-stringview-data-usage)
  if (bam_set1(bam_record,
               qname.length(),
               qname.data(),
               flag,
               tid,
               pos,
               mapq,
               cigar.size(),
               cigar.data(),
               mtid,
               mpos,
               isize,
               seq.length(),
               seq.data(),
               qual.data(),
               l_aux) < 0) {
    throw error::Error("Failed to set BAM record for read '{}': {}", qname, ErrnoMessage());
  }
  // NOLINTEND(bugprone-suspicious-stringview-data-usage)
}

Bam1Ptr BamInit1() {
  Bam1Ptr bam_record = AdoptUnique<bam1_t, Bam1Deleter>(bam_init1());
  if (bam_record == nullptr) {
    throw error::Error("Failed to allocate BAM record");
  }
  return bam_record;
}

Bam1Ptr BamDup1(const bam1_t* const src) {
  Bam1Ptr bam_record = AdoptUnique<bam1_t, Bam1Deleter>(bam_dup1(src));
  if (bam_record == nullptr) {
    throw error::Error("Failed to duplicate BAM record '{}'", bam_get_qname(src));
  }
  return bam_record;
}

HtsItrMultiPtr SamItrMultiQueryS(const hts_idx_t* const idx, sam_hdr_t* const hdr, const std::string& region) {
  HtsItrMultiPtr itr{sam_itr_querys(idx, hdr, region.c_str())};
  if (itr == nullptr) {
    throw error::Error("Failed to create multi-region iterator for '{}': {}", region, ErrnoMessage());
  }
  return itr;
}

HtsItrMultiPtr SamItrMultiQueryI(const hts_idx_t* const idx, const s32 tid, const s64 beg, const s64 end) {
  HtsItrMultiPtr itr{sam_itr_queryi(idx, tid, beg, end)};
  if (itr == nullptr) {
    throw error::Error("Failed to create multi-region iterator for {}-{}:{}: {}", tid, beg, end, ErrnoMessage());
  }
  return itr;
}

void SamIndexBuild3(const fs::path& bam_filename,
                    const fs::path& index_filename,
                    const s32 min_shift,
                    const s32 threads) {
  if (sam_index_build3(bam_filename.c_str(), index_filename.c_str(), min_shift, threads) != 0) {
    throw error::Error("Failed to build index for BAM file '{}'", bam_filename);
  }
}

}  // namespace xoos::io
