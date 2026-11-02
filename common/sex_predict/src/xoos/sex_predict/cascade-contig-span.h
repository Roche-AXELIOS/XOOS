#pragma once

#include <optional>
#include <utility>

#include <xoos/types/int.h>

namespace xoos::sex_predict {

/// The scannable span of a contig after end-trimming (PAR removal plus
/// subtelomeric repeat avoidance, both reference-build agnostic -- see the
/// README's "Pseudoautosomal regions" section), or `std::nullopt` when the
/// contig is too short to survive it. Shared by the depth-profile and
/// heterozygosity phases, which both window-sample within a contig's
/// trimmed interior.
std::optional<std::pair<s64, s64>> TrimmedSpan(s64 chrom_len);

}  // namespace xoos::sex_predict
