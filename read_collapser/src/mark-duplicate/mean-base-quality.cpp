#include "mark-duplicate/mean-base-quality.h"

#include <numeric>
#include <span>

#include <xoos/error/error.h>
#include <xoos/io/htslib-util/htslib-util.h>

namespace xoos::read_collapser {

f64 MeanBaseQ(const AlignmentPtr& alignment) {
  const s32 l_qseq = alignment->record->core.l_qseq;
  if (l_qseq == 0) {
    return 0.0;
  }
  const std::span<const u8> qual{bam_get_qual(alignment->record.get()), static_cast<size_t>(l_qseq)};
  return std::accumulate(qual.begin(), qual.end(), 0.0) / l_qseq;
}

AlignmentPtr FindAlignmentWithMaxMeanBaseQ(const vec<AlignmentPtr>& alignments) {
  if (alignments.empty()) {
    return nullptr;
  }

  // Prefer full-length reads over partial reads.
  // Record the read with the highest mean base quality in each category.
  AlignmentPtr best_full = nullptr;
  f64 best_full_baseq = -1.0;
  AlignmentPtr best_partial = nullptr;
  f64 best_partial_baseq = -1.0;

  for (const auto& alignment : alignments) {
    const f64 baseq = MeanBaseQ(alignment);
    if (alignment->IsPartial()) {
      if (baseq > best_partial_baseq) {
        best_partial = alignment;
        best_partial_baseq = baseq;
      }
    } else {
      if (baseq > best_full_baseq) {
        best_full = alignment;
        best_full_baseq = baseq;
      }
    }
  }
  // If there are any full-length reads, return the one with the highest mean base quality. Otherwise, return the
  // partial read with the highest mean base quality.
  return best_full != nullptr ? best_full : best_partial;
}

}  // namespace xoos::read_collapser
