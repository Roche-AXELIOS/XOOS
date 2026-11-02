#pragma once

#include <xoos/io/metadata-util.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "cli/copy-number-caller-options.h"
#include "segmentation/genomic-segments.h"
#include "segmentation/segmentation-options.h"

namespace xoos::cnc {

using std::size_t;
void ParentSpecificBinarySegmentationMain(const CopyNumberCallerOptions& options);
vec<GenomicSegment> ParentSpecificBinarySegmentation(const Observations& logrs,
                                                     const Observations& dh_vals,
                                                     const vec<GenomicSegment>& seed_segments,
                                                     const SegmentationOptions& segmentation_options,
                                                     const SampleMetadataOptions& sample_metadata_options,
                                                     size_t n_threads,
                                                     const fs::path& logr_segments_out,
                                                     const fs::path& baf_segments_out,
                                                     const io::CommandLineInfo& command_line_info);
}  // namespace xoos::cnc
