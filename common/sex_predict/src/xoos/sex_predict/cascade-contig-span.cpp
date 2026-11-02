#include "cascade-contig-span.h"

namespace xoos::sex_predict {

namespace {

// Bases trimmed from the start and end of every contig before any depth or
// heterozygosity sampling. Two purposes, both reference-build agnostic:
//
//   * PAR removal on the allosomes WITHOUT build detection. PAR1 lies within
//     the first ~2.8 Mb and PAR2 within the last ~330 kb of chrX and chrY alike,
//     on GRCh37 and GRCh38 both, so this trim covers every combination with no
//     coordinate table. That matters because the build may be unrecognized and
//     chrY PAR may or may not have been hard-masked at alignment -- we cannot
//     assume either.
//   * Subtelomeric repeat avoidance, which is good practice on autosomes too.
//
// Applied uniformly to every contig so that all groups get identical treatment
// and the ratios between them stay comparable.
constexpr s64 kContigStartTrim = 3'000'000;
constexpr s64 kContigEndTrim = 500'000;

}  // namespace

std::optional<std::pair<s64, s64>> TrimmedSpan(const s64 chrom_len) {
  const s64 start = kContigStartTrim;
  const s64 end = chrom_len - kContigEndTrim;
  if (end <= start) {
    return std::nullopt;
  }
  return std::make_pair(start, end);
}

}  // namespace xoos::sex_predict
