/**
 * @file hts-idx-internal.h
 * @brief Accessors for htslib index (`hts_idx_t`) internals that htslib does not expose publicly.
 *
 * htslib's `hts_idx_t` is an opaque type. The binning and linear index data needed to estimate
 * per-region workload (virtual file offsets, 16kb linear-index windows) is only reachable through
 * htslib's private struct layout and khash macros. Those are not covered by htslib's API stability
 * guarantees, so all such access is confined to `hts-idx-internal.cpp` behind the functions declared
 * here. Prefer the public wrappers in `htslib-util.h` for anything that does not require these
 * internals.
 *
 * The implementation pins the supported htslib version range with a compile-time check, so an
 * upgrade that could change the struct layout fails to build rather than silently misreading memory.
 */
#pragma once

#include <htslib/hts.h>

#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::io {

/**
 * Number of reference sequences covered by the index.
 *
 * @param idx Loaded BAM index. Must not be null.
 */
s32 BamIndexNumRef(const hts_idx_t* idx);

/**
 * Number of unplaced (no-coordinate / unmapped) reads recorded in the index.
 *
 * @param idx Loaded BAM index. Must not be null.
 */
u64 BamIndexNumNoCoor(const hts_idx_t* idx);

/**
 * Linear-index virtual file offsets for a reference sequence.
 *
 * The linear index records the virtual file offset at the start of each 16kb genomic window. The
 * returned offsets are ordered by window. Returns an empty vector when `tid` is out of range or the
 * reference has no linear-index entries.
 *
 * @param idx Loaded BAM index. Must not be null.
 * @param tid Reference sequence ID.
 */
vec<u64> BamIndexLinearOffsets(const hts_idx_t* idx, s32 tid);

/**
 * Maximum virtual file offset for a reference sequence: the end of reads mapped to it.
 *
 * @param idx Loaded BAM index. Must not be null.
 * @param tid Reference sequence ID.
 * @return Maximum virtual file offset, or 0 if `tid` is invalid or has no data.
 */
u64 BamIndexMaxMappedOffset(const hts_idx_t* idx, s32 tid);

/**
 * Minimum virtual file offset for a reference sequence: the start of reads mapped to it.
 *
 * @param idx Loaded BAM index. Must not be null.
 * @param tid Reference sequence ID.
 * @return Minimum virtual file offset, or 0 if `tid` is invalid or has no data.
 */
u64 BamIndexMinMappedOffset(const hts_idx_t* idx, s32 tid);

}  // namespace xoos::io
