#include "ref-region.h"

#include <xoos/io/htslib-util/htslib-util.h>

namespace xoos::svc {

StrMap<RefRegion> GetRefRegions(const AlignmentReader& reader) {
  // Given a BAM file, header and index this function reads the header and computes a
  // map of RefRegion structs that store information relating to the contigs
  // contained within the BAM header.
  StrMap<RefRegion> contig_map;
  io::Bam1Ptr b(bam_init1());
  // Iterate through each contig in the header one by one. The start position recorded for the contig is the start
  // position of the first alignment found for the contig in the BAM. If there are no alignments found for the contig in
  // the BAM we skip and don't record an entry in the map. This way only contigs that have alignments are recorded.
  for (int i = 0; i < reader.hdr->n_targets; ++i) {
    const std::string chrom = io::SamHdrTid2Name(reader.hdr.get(), i);
    const io::HtsItrPtr tmp_itr(io::SamItrQuerySNoThrow(reader.idx.get(), reader.hdr.get(), chrom));
    if (tmp_itr == nullptr) {
      continue;
    }
    const auto has_next = io::SamItrNext(reader.fp.get(), tmp_itr.get(), b.get());
    if (!has_next) {
      continue;
    }
    const u64 start = b->core.pos - 1 < 0 ? 0UL : static_cast<u64>(b->core.pos - 1);
    const u64 len = io::SamHdrTid2Length(reader.hdr.get(), i);
    const u64 end = len;
    const RefRegion region{
        .name = chrom,
        .ref_id = static_cast<u32>(i),
        .length = len,
        .start_position = start,
        .end_position = end,
    };
    contig_map.try_emplace(chrom, region);
  }
  return contig_map;
}

}  // namespace xoos::svc
