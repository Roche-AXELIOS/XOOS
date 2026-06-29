#include "mark-duplicate/mean-base-quality.h"

#include <numeric>

#include <xoos/error/error.h>
#include <xoos/io/htslib-util/htslib-util.h>

namespace xoos::read_collapser {

f64 MeanBaseQ(const AlignmentPtr& alignment) {
  const u8* qual = bam_get_qual(alignment->record.get());
  const s32 l_qseq = alignment->record->core.l_qseq;
  if (l_qseq == 0) {
    return 0.0;
  }
  return std::accumulate(qual, qual + l_qseq, 0.0) / l_qseq;
}

AlignmentPtr FindAlignmentWithMaxMeanBaseQ(const vec<AlignmentPtr>& alignments) {
  if (alignments.empty()) {
    return nullptr;
  }

  auto max_alignment = alignments.front();
  auto max_mean_baseq = MeanBaseQ(max_alignment);
  for (size_t i = 1; i < alignments.size(); ++i) {
    const auto& i_alignment = alignments.at(i);
    if (const f64 i_mean_baseq = MeanBaseQ(i_alignment); i_mean_baseq > max_mean_baseq) {
      max_alignment = i_alignment;
      max_mean_baseq = i_mean_baseq;
    }
  }

  return max_alignment;
}

}  // namespace xoos::read_collapser
