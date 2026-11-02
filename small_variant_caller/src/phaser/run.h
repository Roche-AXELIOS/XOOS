#pragma once

#include "params.h"

namespace xoos::svc::phaser {

/**
 * @brief Run read-based phasing and write the annotated output VCF.
 *
 * The input VCF is read, candidate phase groups are built, reads are examined, and the
 * resulting phase annotations and merged records are written to the requested output.
 *
 * @param params Input paths, phasing thresholds, and output options.
 */
// Public entry point called from apps/phase-variants.cpp and the phaser tests; Sonar analyses only the src tree.
void Run(const Params& params);  // NOSONAR(cpp:S5536)

}  // namespace xoos::svc::phaser
