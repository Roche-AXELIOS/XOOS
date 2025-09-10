#include "htslib-util.h"

#include <xoos/error/error.h>

#include "xoos/io/htslib-util/htslib-ptr.h"
#include "xoos/io/malloc-ptr.h"

namespace xoos::io {

HtsIdxPtr HtsIdxLoad(const fs::path& fn, int fmt) {
  HtsIdxPtr idx{hts_idx_load(fn.c_str(), fmt)};
  if (idx == nullptr) {
    throw error::Error("Failed to read index for HTS file: '{}'", fn);
  }
  return idx;
}

HtsIdxPtr SamIndexLoad(htsFile* fp, const fs::path& fn) {
  HtsIdxPtr idx{sam_index_load(fp, fn.c_str())};
  if (idx == nullptr) {
    throw error::Error("Failed to read index for BAM/CRAM file: '{}'", fn);
  }
  return idx;
}

void HtsFlush(htsFile* fp) {
  if (hts_flush(fp) < 0) {
    throw error::Error("Failed to flush HTS file: '{}'", fp->fn);
  }
}

HtsFilePtr HtsOpenFormat(const fs::path& fn, const std::string& mode, const htsFormat* fmt) {
  HtsFilePtr fp{hts_open_format(fn.c_str(), mode.c_str(), fmt)};
  if (fp == nullptr) {
    throw error::Error("Failed to open HTS file '{}': {}", fn, strerror(errno));
  }
  return fp;
}

static std::string ToString(const hts_reglist_t* reglist) {
  if (reglist != nullptr && reglist->count > 0) {
    return reglist->reg;
  }
  return "No regions";
}

HtsItrPtr SamItrRegion(const hts_idx_t* idx, sam_hdr_t* hdr, hts_reglist_t* reglist, const u32 reg_count) {
  HtsItrPtr itr_ptr(sam_itr_regions(idx, hdr, reglist, reg_count));
  // Check if the iterator was created successfully if not, throw an error with the appropriate message
  if (itr_ptr == nullptr) {
    throw error::Error(
        "Failed to create iterator for regions in BAM file with reglist '{}': {}", ToString(reglist), strerror(errno));
  }
  return itr_ptr;
}

// This function reads the next record from a multi-region iterator, returning -1 if there are no more records to read
// or 0 if there is more. If an error occurs (<-1), it throws an error with the appropriate message.
int SamItrMultiNext(htsFile* fp, hts_itr_t* itr, bam1_t* record) {
  const int return_value = sam_itr_multi_next(fp, itr, record);
  if (return_value < -1) {
    throw error::Error("Error reading alignment '{}': {}", fp->fn, strerror(errno));
  }
  return return_value;
}

HtsFilePtr HtsOpen(const fs::path& fn, const std::string& mode) {
  HtsFilePtr fp{hts_open(fn.c_str(), mode.c_str())};
  if (fp == nullptr) {
    throw error::Error("Failed to open HTS file '{}': {}", fn, strerror(errno));
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
const char* SamHdrTid2Name(const sam_hdr_t* hdr, const s32 tid) {
  const auto* name = sam_hdr_tid2name(hdr, tid);
  if (name == nullptr) {
    throw error::Error("Failed to get name for tid {} from BAM header", tid);
  }
  return name;
}

s64 SamHdrTid2Length(const sam_hdr_t* hdr, const s32 tid) {
  const auto len = sam_hdr_tid2len(hdr, tid);
  if (len <= 0) {
    throw error::Error("Failed to get length for tid {} from BAM header", tid);
  }
  return len;
}

s32 SamHdrNRef(const sam_hdr_t* hdr) {
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

void SamWrite1(samFile* fp, const sam_hdr_t* hdr, const bam1_t* record) {
  if (sam_write1(fp, hdr, record) < 0) {
    throw error::Error("Failed to write BAM record '{}' to '{}': {}", bam_get_qname(record), fp->fn, strerror(errno));
  }
}

SamHdrPtr SamHdrRead(samFile* fp) {
  SamHdrPtr hdr{sam_hdr_read(fp)};
  if (hdr == nullptr) {
    throw error::Error("Failed to read header from BAM file: {}", fp->fn);
  }
  return hdr;
}

void SamHdrWrite(samFile* fp, const sam_hdr_t* hdr) {
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

void SamHdrAddSqLine(sam_hdr_t* hdr, const SqHdrLine& line) {
  const auto ln = std::to_string(line.ln);
  if (sam_hdr_add_line(hdr, "SQ", "SN", line.sn.c_str(), "LN", ln.c_str(), nullptr) == -1) {
    throw error::Error("Failed to add SQ header line for '{}'", line.sn);
  }
}

void SamHdrAddPgLine(sam_hdr_t* hdr, const PgHdrLine& line) {
  if (line.id.empty()) {
    throw error::Error("Cannot add PG header line with empty ID");
  }
  if (line.program_name.empty()) {
    throw error::Error("Cannot add PG header line with empty program name");
  }
  // Build the @PG header line as a single string
  std::string pg_string = fmt::format("@PG\tID:{}\tPN:{}", line.id, line.program_name);
  if (line.version_number) {
    pg_string += fmt::format("\tVN:{}", line.version_number.value());
  }
  if (line.command_line) {
    pg_string += fmt::format("\tCL:{}", line.command_line.value());
  }
  if (line.previous_program_id) {
    pg_string += fmt::format("\tPP:{}", line.previous_program_id.value());
  }

  if (sam_hdr_add_lines(hdr, pg_string.c_str(), pg_string.length()) == -1) {
    throw error::Error("Failed to add PG header line for ID '{}'", line.id);
  }
}

void SamIdxSave(htsFile* fp) {
  if (sam_idx_save(fp) < 0) {
    throw error::Error("Failed to save BAM index for '{}'", fp->fn);
  }
}

void SamIdxInit(htsFile* fp, sam_hdr_t* hdr, const int min_shift, const std::string& idx_fn) {
  if (sam_idx_init(fp, hdr, min_shift, idx_fn.c_str()) < 0) {
    throw error::Error("Failed to initialize BAM index for '{}'", fp->fn);
  }
}

bool SamItrNext(htsFile* fp, hts_itr_t* itr, bam1_t* record) {
  const auto ret = sam_itr_next(fp, itr, record);
  if (ret < -1) {
    throw error::Error("Failed to get next record from BAM '{}': {}", fp->fn, strerror(errno));
  }
  // Return true if a record was read, false if no more records are available.
  return ret >= 0;
}

HtsItrPtr SamItrQueryI(hts_idx_t* idx, s32 tid, s64 beg, s64 end) {
  auto itr = HtsItrPtr{sam_itr_queryi(idx, tid, beg, end)};
  if (itr == nullptr) {
    throw error::Error("Failed to create iterator {}-{}:{}: {}", tid, beg, end, strerror(errno));
  }
  return itr;
}

HtsItrPtr SamItrQueryS(const hts_idx_t* idx, sam_hdr_t* hdr, const std::string& region) {
  auto itr = HtsItrPtr{sam_itr_querys(idx, hdr, region.c_str())};
  if (itr == nullptr) {
    throw error::Error("Failed to create iterator {} : {}", region, strerror(errno));
  }
  return itr;
}

HtsItrPtr SamItrQuerySNoThrow(const hts_idx_t* idx, sam_hdr_t* hdr, const std::string& region) {
  return HtsItrPtr{sam_itr_querys(idx, hdr, region.c_str())};
}

s32 SamHdrName2Tid(sam_hdr_t* hdr, const std::string& target) {
  const auto tid = sam_hdr_name2tid(hdr, target.c_str());
  if (tid == -1) {
    throw error::Error("Failed to find target '{}' in BAM header", target);
  }
  if (tid == -2) {
    throw error::Error("Failed to parse BAM header");
  }
  return tid;
}

void BamAuxUpdateStr(bam1_t* record, const std::string_view& tag, const std::string_view& value) {
  if (bam_aux_update_str(record, tag.data(), static_cast<s32>(value.length() + 1), value.data()) != 0) {
    throw error::Error("Failed to update BAM aux tag '{}'", tag);
  }
}

void BamAuxAppend(bam1_t* record, const std::string& name, const u32 value) {
  if (bam_aux_append(record, name.c_str(), 'I', sizeof(u32), reinterpret_cast<const u8*>(&value)) == -1) {
    throw error::Error("Failed to append BAM aux tag '{}' with value '{}': {}", name, value, strerror(errno));
  }
}

void BamAuxAppend(bam1_t* record, const std::string& name, const std::string& value) {
  const auto len = static_cast<int>(value.size()) + 1;
  if (bam_aux_append(record, name.c_str(), 'Z', len, reinterpret_cast<const u8*>(value.c_str())) == -1) {
    throw error::Error("Failed to append BAM aux tag '{}' with value '{}': {}", name, value, strerror(errno));
  }
}

template <>
std::optional<u32> BamAuxGet(const bam1_t* record, const std::string& tag) {
  auto* const tag_ptr = bam_aux_get(record, tag.c_str());
  if (tag_ptr == nullptr) {
    return std::nullopt;
  }
  return bam_aux2i(tag_ptr);
}

template <>
std::optional<std::string> BamAuxGet(const bam1_t* record, const std::string& tag) {
  auto* const tag_ptr = bam_aux_get(record, tag.c_str());
  if (tag_ptr == nullptr) {
    return std::nullopt;
  }
  const auto* str = bam_aux2Z(tag_ptr);
  if (str == nullptr) {
    return std::nullopt;
  }
  return std::string(str);
}

void BamSet1(bam1_t* bam,
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
  if (bam_set1(bam,
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
    throw error::Error("Failed to set BAM record for read '{}': {}", qname, strerror(errno));
  }
}

}  // namespace xoos::io
