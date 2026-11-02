#pragma once

#include <atomic>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <xoos/types/int.h>

namespace xoos::sex_predict {

/// One chunk of one chromosome to be counted by a single scan task. `start`/
/// `end` are 0-based, half-open [start, end), the XOOS convention.
/// `region_index` routes the chunk's count back to the right logical region
/// (e.g. 0=chrX/1=chrY for the ratio scan, or an autosome's position in the
/// sampled list for the depth-profile scan) -- callers assign the meaning.
/// Shared by the default method's ratio scan and the cascade's depth-profile
/// phase, which both tile chromosomes into chunks and count filtered reads
/// per chunk via the same taskflow pattern.
struct ScanChunk {
  s32 tid;
  s64 start;
  s64 end;
  u32 region_index;
};

/// Count primary, mapped, non-duplicate, non-QC-fail reads with MAPQ >=
/// `min_mapq` whose start falls in the chunk's [start, end) 0-based
/// half-open range. If `min_reads > 0`, stops early once the shared
/// `total_counted` reaches it.
u64 CountFilteredReads(htsFile* bam_file,
                       const hts_idx_t* idx,
                       const ScanChunk& chunk,
                       u8 min_mapq,
                       u64 min_reads,
                       std::atomic<u64>& total_counted);

/// Length of contig `tid` in `header`.
s64 TargetLen(const sam_hdr_t* header, s32 tid);

}  // namespace xoos::sex_predict
