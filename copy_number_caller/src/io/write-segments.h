#pragma once

#include <optional>

#include <xoos/io/metadata-util.h>
#include <xoos/sex_predict/sex.h>
#include <xoos/types/fs.h>
#include <xoos/types/vec.h>

#include "observations/observations.h"
#include "segmentation/genomic-segments.h"
#include "segmentation/segment-type.h"

namespace xoos::cnc {

void WriteSegments(const fs::path& out_fname,
                   const vec<GenomicSegment>& segments,
                   const SegmentType& mode,
                   const std::optional<io::CommandLineInfo>& command_line_info);

void WriteLogRSegments(const vec<GenomicSegment>& segments,
                       const fs::path& logr_segments_out,
                       const Observations& logrs,
                       std::string_view sample_id,
                       sex_predict::Sex sex,
                       const io::CommandLineInfo& command_line_info);

void WriteBafSegments(const vec<GenomicSegment>& segments,
                      const fs::path& baf_segments_out,
                      const io::CommandLineInfo& command_line_info);

}  // namespace xoos::cnc
