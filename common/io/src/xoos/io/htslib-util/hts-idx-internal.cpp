#include "xoos/io/htslib-util/hts-idx-internal.h"

#include <algorithm>
#include <gsl/gsl>  // NOLINT(misc-include-cleaner)
#include <limits>

#include <htslib/hts.h>
#include <htslib/khash.h>

#include <xoos/types/int.h>
#include <xoos/types/vec.h>

// htslib does not expose the layout of `hts_idx_t`, the binning index, or the linear index in its
// public headers. The definitions below are copied verbatim from htslib's hts.c so we can read the
// virtual file offsets needed for workload estimation. They are private implementation details and
// are deliberately confined to this translation unit.
//
// Source: htslib hts.c, v1.23.1 (MIT/Expat License).
//
// HTS_VERSION encodes the htslib release as MAJOR*100000 + MINOR*100 + PATCH (e.g. 1.23.1 -> 102301).
// The struct layout below is stable across 1.19 through 1.23.1; the only change in that range is the
// trailing `otf_fp` field added in 1.20, which precedes nothing we read. The guard fails the build on
// an untested htslib so the layout can be re-verified before bumping the upper bound, rather than
// silently misreading memory.
#if !defined(HTS_VERSION)
#error "HTS_VERSION is not defined; cannot verify htslib internal struct layout."
#endif
#if HTS_VERSION < 101900 || HTS_VERSION > 102301
#error \
    "Unverified htslib version for BAM index internals. Re-check the hts_idx_t/bins_t/lidx_t " \
    "layout in htslib's hts.c against hts-idx-internal.cpp, then widen this version guard."
#endif

// NOLINTBEGIN

namespace {

// Finds the special meta bin: ((1 << (3 * n_lvls + 3)) - 1) / 7 + 1
#define META_BIN(idx) ((idx)->n_bins + 1)

typedef struct {
  int32_t m, n;
  uint64_t loff;
  hts_pair64_t* list;
} bins_t;

KHASH_MAP_INIT_INT(bin, bins_t)
typedef khash_t(bin) bidx_t;

typedef struct {
  hts_pos_t n, m;
  uint64_t* offset;
} lidx_t;

// Mirror of htslib's private `struct hts_idx_t`. Named distinctly so it does not collide with the
// public opaque `hts_idx_t` forward declaration in <htslib/hts.h>.
struct HtsIdxLayout {
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

  BGZF* otf_fp;
};

}  // namespace

// NOLINTEND

namespace xoos::io {

namespace {

// Reinterpret the public opaque index handle as the private layout copied above. All internal field
// access goes through this single cast so the unsafe reinterpretation lives in one place.
const HtsIdxLayout* Internal(const hts_idx_t* const idx) {
  // NOSONAR, NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast) - sole entry point to htslib internals
  return reinterpret_cast<const HtsIdxLayout*>(idx);
}

// Span over the per-reference linear index. Empty when the linear index is absent.
gsl::span<const lidx_t> LinearIndices(const HtsIdxLayout* const internal) {
  if (internal->lidx == nullptr) {
    return {};
  }
  return {internal->lidx, ToUnsigned(std::max(internal->n, 0))};
}

// Span over the per-reference binning index pointers.
gsl::span<bidx_t* const> BinIndices(const HtsIdxLayout* const internal) {
  return {internal->bidx, ToUnsigned(std::max(internal->n, 0))};
}

// Span over the linear-index offsets for one reference sequence. Empty when there are none.
gsl::span<const u64> LinearOffsetSpan(const lidx_t& lidx) {
  if (lidx.offset == nullptr || lidx.n <= 0) {
    return {};
  }
  return {lidx.offset, ToUnsigned(lidx.n)};
}

// Span over one binning-index chunk list. Empty when there are none.
gsl::span<const hts_pair64_t> ChunkSpan(const bins_t& bin) {
  if (bin.list == nullptr || bin.n <= 0) {
    return {};
  }
  return {bin.list, ToUnsigned(bin.n)};
}

// Read the metadata-bin chunk for a reference, if present. The meta bin stores the reference's
// overall begin/end virtual offsets in a single chunk (.u = begin, .v = end).
const hts_pair64_t* MetaChunk(const HtsIdxLayout* const internal, const bidx_t* const bidx) {
  const auto meta = kh_get(bin, bidx, META_BIN(internal));
  if (meta != kh_end(bidx) && kh_val(bidx, meta).n > 0) {
    return gsl::span(kh_val(bidx, meta).list, ToUnsigned(kh_val(bidx, meta).n)).data();
  }
  return nullptr;
}

// Whether to take the minimum or maximum mapped virtual offset for a reference sequence.
enum class OffsetBound : u8 {
  kMin,
  kMax
};

// Walk a reference sequence's binning and linear index to find its minimum or maximum virtual file
// offset. Returns 0 when the reference has no data. Shared by the public min/max accessors so the
// internal traversal exists in exactly one place.
u64 BoundMappedOffset(const HtsIdxLayout* const internal, const s32 tid, const OffsetBound bound) {
  if (tid < 0 || tid >= internal->n) {
    return 0;
  }
  const auto* const bidx = BinIndices(internal)[tid];
  if (bidx == nullptr) {
    return 0;
  }

  // The metadata bin, when present, records the reference's begin (.u) and end (.v) offsets directly.
  if (const auto* const meta = MetaChunk(internal, bidx); meta != nullptr) {
    return bound == OffsetBound::kMin ? meta->u : meta->v;
  }

  // Otherwise scan every chunk and every linear-index window for the extreme offset.
  const auto combine = [bound](const u64 acc, const u64 value) {
    return bound == OffsetBound::kMin ? std::min(acc, value) : std::max(acc, value);
  };
  auto result = bound == OffsetBound::kMin ? std::numeric_limits<u64>::max() : u64{0};
  bool seen = false;
  for (khint_t k = 0; k < kh_end(bidx); ++k) {
    if (kh_exist(bidx, k) == 0) {
      continue;
    }
    for (const auto& chunk : ChunkSpan(kh_val(bidx, k))) {
      result = combine(result, bound == OffsetBound::kMin ? chunk.u : chunk.v);
      seen = true;
    }
  }
  if (const auto lidx = LinearIndices(internal); !lidx.empty()) {
    const auto offsets = LinearOffsetSpan(lidx[tid]);
    for (const auto offset : offsets) {
      result = combine(result, offset);
      seen = true;
    }
  }
  return seen ? result : 0;
}

}  // namespace

s32 BamIndexNumRef(const hts_idx_t* const idx) {
  return Internal(idx)->n;
}

u64 BamIndexNumNoCoor(const hts_idx_t* const idx) {
  return Internal(idx)->n_no_coor;
}

vec<u64> BamIndexLinearOffsets(const hts_idx_t* const idx, const s32 tid) {
  const auto* const internal = Internal(idx);
  const auto lidx = LinearIndices(internal);
  if (tid < 0 || tid >= internal->n || lidx.empty()) {
    return {};
  }
  const auto offsets = LinearOffsetSpan(lidx[tid]);
  return {offsets.begin(), offsets.end()};
}

u64 BamIndexMaxMappedOffset(const hts_idx_t* const idx, const s32 tid) {
  return BoundMappedOffset(Internal(idx), tid, OffsetBound::kMax);
}

u64 BamIndexMinMappedOffset(const hts_idx_t* const idx, const s32 tid) {
  return BoundMappedOffset(Internal(idx), tid, OffsetBound::kMin);
}

}  // namespace xoos::io
