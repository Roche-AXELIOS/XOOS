#pragma once

#include "cli/copy-number-caller-options.h"
#include "observations/observations.h"

namespace xoos::cnc {

/**
 * @brief Extract tumor-only B-allele frequency (BAF) observations in memory.
 *
 * Scans the single-sample tumor VCF, keeps likely germline heterozygous SNP
 * sites confirmed against the population VCF, applies depth/BAF filters, and
 * excludes any site overlapping the optional SNP blocklist BED. Returns the
 * ref and alt allele-depth observations without touching the filesystem, so the
 * result can be consumed directly by downstream in-memory stages or written out
 * by the caller.
 *
 * @param options Copy-number-caller options (tumor/population VCFs, filters,
 *                optional snp_blocklist_bed_fname, sample metadata).
 * @return RefAltObservations holding the per-site ref and alt depths.
 */
RefAltObservations ExtractTumorOnlyBafs(const CopyNumberCallerOptions& options);

/// Extract tumor-only BAFs and write them to the module's default BAF output file.
void TumorOnlyBafsMain(const CopyNumberCallerOptions& options);

}  // namespace xoos::cnc
