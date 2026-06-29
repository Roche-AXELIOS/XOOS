#pragma once

#include <map>

#include <xoos/types/fs.h>
#include <xoos/types/vec.h>

#include "core/interval.h"

namespace xoos::alignment_metrics {

/**
 * A mapping from chromosome name to a vector of intervals on that chromosome.
 */
using RegionsByChromosome = std::map<std::string, vec<Interval>, std::less<>>;

/**
 * Read regions from a BED file, skipping contigs not present in the BAM file header if a BAM file is provided.
 */
RegionsByChromosome GetRegionsFromBed(const fs::path& bed_filename, const fs::path& bam_filename);

/**
 * Read regions from a BAM file header. The resulting regions will cover every reference contig
 * present in the BAM file.
 */
RegionsByChromosome GetRegionsFromBam(const fs::path& bam_filename);

}  // namespace xoos::alignment_metrics
