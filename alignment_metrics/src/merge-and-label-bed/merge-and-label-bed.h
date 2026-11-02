#pragma once

#include "core/region-labels-set.h"
#include "io/bed-reader.h"
#include "merge-and-label-bed/merge-and-label-bed-options.h"

namespace xoos::alignment_metrics::merge_and_label_bed {

/**
 * Write merged, labeled intervals to a 4-column BED file (`chrom  start  end  label1,label2,...`),
 * ordered by chromosome then position. Intervals with no labels are written with an empty 4th
 * column omitted (3 columns). A `##RocheCommandLine` metadata comment is written first when
 * `command_line_info` is populated.
 */
void WriteLabeledBed(const RegionsByChromosome& regions,
                     const RegionLabelsSet& labels,
                     const fs::path& output,
                     const io::CommandLineInfo& command_line_info);

/**
 * Entry point for the merge_and_label_bed workflow: merge/label all inputs and write the result.
 */
void MergeAndLabelBed(const MergeAndLabelBedOptions& options);

}  // namespace xoos::alignment_metrics::merge_and_label_bed
