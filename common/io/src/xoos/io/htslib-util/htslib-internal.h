#pragma once

// NOLINTBEGIN

#include <htslib/khash.h>

#include <xoos/error/error.h>
#include <xoos/io/htslib-util/htslib-util.h>

namespace xoos::io::htslib {

// The following is copied from htslib 1.21, these structs are not
// part of the public API, but are used internally by htslib. But
// they are useful to access internal details of the index. These structures
// are unlikely to ever need updating, but they may at some point be modified by
// the upstream htslib project.
//
// https://github.com/samtools/htslib/blob/1.21/hts.c#L2207
// The MIT/Expat License

// Finds the special meta bin
//  ((1<<(3 * n_lvls + 3)) - 1) / 7 + 1
#define META_BIN(idx) ((idx)->n_bins + 1)

typedef struct {
  int32_t m, n;
  uint64_t loff;
  hts_pair64_t* list;
} bins_t;

// Forward declaration of the internal struct fields for clarity.
// This avoids needing to include private htslib headers.
// These definitions are based on the htslib source code.
KHASH_MAP_INIT_INT(bin, bins_t)
typedef khash_t(bin) bidx_t;

typedef struct {
  hts_pos_t n, m;
  uint64_t* offset;
} lidx_t;

struct hts_idx_t {
  int fmt, min_shift, n_lvls, n_bins;
  uint32_t l_meta;
  int32_t n, m;
  uint64_t n_no_coor;
  bidx_t** bidx;
  lidx_t* lidx;
  uint8_t* meta;  // MUST have a terminating NUL on the end
  int tbi_n, last_tbi_tid;

  struct {
    uint32_t last_bin, save_bin;
    hts_pos_t last_coor;
    int last_tid, save_tid, finished;
    uint64_t last_off, save_off;
    uint64_t off_beg, off_end;
    uint64_t n_mapped, n_unmapped;
  } z;  // keep internal states

  BGZF* otf_fp;  // Index on-the-fly output file
};

struct HtsIdxDeleter {
  void operator()(hts_idx_t* idx) {
    if (idx != nullptr) {
      hts_idx_destroy(reinterpret_cast<::hts_idx_t*>(idx));
    }
  }
};

using HtsIdxPtr = std::unique_ptr<hts_idx_t, HtsIdxDeleter>;

HtsIdxPtr HtsIdxLoad(const fs::path& bam, int fmt) {
  HtsIdxPtr idx{reinterpret_cast<hts_idx_t*>(hts_idx_load(bam.c_str(), fmt))};
  if (idx == nullptr) {
    throw error::Error("Failed to read index from BAM file: {}", bam);
  }
  return idx;
}

// NOLINTEND

}  // namespace xoos::io::htslib
