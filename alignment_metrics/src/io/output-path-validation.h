#pragma once

#include <string>

#include <xoos/types/fs.h>
#include <xoos/types/vec.h>

#include "alignment-metrics-options.h"
#include "core/region-labels-set.h"
#include "metadata/dataset-metadata.h"

namespace xoos::alignment_metrics {

vec<fs::path> CollectOutputPaths(const AlignmentMetricsOptions& options,
                                 const DatasetMetadata& metadata,
                                 const fs::path& out_dir,
                                 const std::string& output_prefix);

void ValidateOutputFilesDoNotExist(const AlignmentMetricsOptions& options,
                                   const DatasetMetadata& metadata,
                                   const RegionLabelsSet& labels,
                                   const std::string& output_prefix);

}  // namespace xoos::alignment_metrics
