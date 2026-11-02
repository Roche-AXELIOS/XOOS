#pragma once

#include <xoos/types/fs.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>

#include "core/interval.h"
#include "core/region-labels-set.h"

namespace xoos::alignment_metrics {

/**
 * A mapping from chromosome name to a vector of intervals on that chromosome.
 */
using RegionsByChromosome = StrMap<vec<Interval>>;

/**
 * Read regions from a BED file, skipping contigs not present in the BAM file header if a BAM file is provided.
 * Extra BED columns are ignored; the returned intervals are unlabeled (labels == 0).
 */
RegionsByChromosome GetRegionsFromBed(const fs::path& bed_filename, const fs::path& bam_filename);

/**
 * Read regions from a BED file, optionally interpreting the 4th column as a comma-separated list of
 * region labels used to stratify metrics.
 *
 * When `stratify_by_label` is true, every region must have a non-empty 4th column; each label in it is
 * registered in `registry` and its bit is set in the interval's `labels` bit vector. When
 * `stratify_by_label` is false, extra columns are ignored and `registry` is left untouched (all returned
 * intervals are unlabeled).
 *
 * Contigs not present in the BAM header are skipped with a warning.
 *
 * @param bed_filename Path to the BED file containing regions of interest.
 * @param bam_filename Path to the BAM file used for chromosome validation.
 * @param registry Region labels set populated with any labels encountered (only when `stratify_by_label`).
 * @param stratify_by_label Whether to parse the 4th column as region labels.
 * @throws error::Error if the BED file is missing, a labeled region lacks a 4th column, or a label
 *         name is invalid / exceeds the supported label count.
 */
RegionsByChromosome GetRegionsFromBed(const fs::path& bed_filename,
                                      const fs::path& bam_filename,
                                      RegionLabelsSet& registry,
                                      bool stratify_by_label);

/**
 * Read regions from a BAM file header. The resulting regions will cover every reference contig
 * present in the BAM file.
 */
RegionsByChromosome GetRegionsFromBam(const fs::path& bam_filename);

}  // namespace xoos::alignment_metrics
