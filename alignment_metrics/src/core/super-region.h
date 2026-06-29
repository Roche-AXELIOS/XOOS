#pragma once

#include <string>

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "core/interval.h"

namespace xoos::alignment_metrics {

/**
 * A super region is a contiguous genomic interval that contains one or more subregions (regions of interest).
 * Super regions are used to group together nearby regions for more efficient processing of alignments.
 */
struct SuperRegion {
  std::string chromosome;
  s64 start;
  s64 end;
  vec<Interval> subregions;
};

// Creates a SuperRegion from a chromosome name and a vector of subregions.
SuperRegion CreateSuperRegion(const std::string& chromosome, const vec<Interval>& sub_regions);

// Creates an htslib region list from a SuperRegion for BAM file querying.
io::HtsRegListPtr CreateHtsRegList(const SuperRegion& super_region);

}  // namespace xoos::alignment_metrics
