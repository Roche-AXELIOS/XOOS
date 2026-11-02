#include "scan-chunk.h"

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/io/htslib-util/read-util.h>

namespace xoos::sex_predict {

namespace {

// SAM spec sentinel for "mapping quality not available"; excluded
// explicitly since it's never `< min_mapq`.
constexpr u8 kMapQNotAvailable = 255;

}  // namespace

u64 CountFilteredReads(htsFile* const bam_file,
                       const hts_idx_t* const idx,
                       const ScanChunk& chunk,
                       const u8 min_mapq,
                       const u64 min_reads,
                       std::atomic<u64>& total_counted) {
  // Unsigned to satisfy clang-tidy's hicpp-signed-bitwise (rec->core.flag is unsigned).
  constexpr u32 kExcludeFlags = ToUnsigned(BAM_FUNMAP) | ToUnsigned(BAM_FSECONDARY) | ToUnsigned(BAM_FQCFAIL) |
                                ToUnsigned(BAM_FDUP) | ToUnsigned(BAM_FSUPPLEMENTARY);
  u64 count = 0;
  const io::Bam1Ptr read = io::BamInit1();
  const io::HtsItrPtr itr = io::SamItrQueryI(idx, chunk.tid, chunk.start, chunk.end);
  while (io::SamItrNext(bam_file, itr.get(), read.get())) {
    const bam1_t* const rec = read.get();
    if ((rec->core.flag & kExcludeFlags) != 0 || rec->core.qual < min_mapq || rec->core.qual == kMapQNotAvailable) {
      continue;
    }
    // Count a read only in the chunk that contains its start position, so a
    // read spanning a chunk boundary is not double-counted by both chunks.
    if (rec->core.pos < chunk.start || rec->core.pos >= chunk.end) {
      continue;
    }
    ++count;
    if (min_reads > 0 && total_counted.fetch_add(1, std::memory_order_relaxed) + 1 >= min_reads) {
      break;
    }
  }
  return count;
}

s64 TargetLen(const sam_hdr_t* const header, const s32 tid) {
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic) — htslib C array access
  return ToSigned(header->target_len[tid]);
}

}  // namespace xoos::sex_predict
