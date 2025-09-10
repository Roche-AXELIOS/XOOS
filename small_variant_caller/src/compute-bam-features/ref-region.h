#pragma once

#include <string>

#include <xoos/error/error.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>

#include "compute-bam-features/alignment-reader.h"

namespace xoos::svc {

struct RefRegion {
  // name of reference region in the BAM header. Example chr1, chr2, chr3
  std::string name;
  // ref_id assigned to the reference region in the BAM header
  u32 ref_id{};
  // length of the reference region in the BAM header
  u64 length{};
  // 0-based start position of the first alignment of the reference
  u64 start_position{0};
  // 0-based end position of the last alignment of the reference
  u64 end_position{0};
};

StrMap<RefRegion> GetRefRegions(const AlignmentReader& reader);

}  // namespace xoos::svc
