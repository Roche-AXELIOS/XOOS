#pragma once

#include <xoos/types/vec.h>

#include "variant.h"

namespace xoos::svc::phaser {

/**
 * @brief Build connected phase blocks and assign haplotypes from pairwise phasing results.
 *
 * Phased pair results are treated as graph edges. Connected components become phase blocks,
 * and cis/trans evidence is used to assign each variant to one of two haplotypes.
 *
 * @param results Pairwise phasing results with variant indices into `variants`.
 * @param variants Variants used to order the resulting phase blocks by genomic position.
 * @return Phase blocks containing at least two connected variants.
 */
vec<PhaseBlock> BuildPhaseBlocks(const vec<PhasingResult>& results, const vec<Variant>& variants);

}  // namespace xoos::svc::phaser
