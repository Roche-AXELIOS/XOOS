#pragma once

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::yc_decode {

struct HtslibCigarOp {
  u8 op{};     // htslib CIGAR op, e.g. `BAM_CMATCH`
  u32 len{0};  // length of the CIGAR op
  auto operator<=>(const HtslibCigarOp&) const = default;
};

struct HtslibCigarAlignInfo {
  u8 op{};           // htslib CIGAR op, e.g. `BAM_CMATCH`
  u64 read_pos{0};   // 0-based start position in the read
  u64 ref_pos{0};    // 0-based start position in the reference
  u64 read_span{0};  // number of read bases consumed by the CIGAR op
  u64 ref_span{0};   // number of reference bases consumed by the CIGAR op
  auto operator<=>(const HtslibCigarAlignInfo&) const = default;
};

// utility functions to handle htslib CIGAR operators
bool ConsumesRead(u32 op);
bool ConsumesReference(u32 op);
vec<HtslibCigarAlignInfo> GetCigarAlignInfos(const bam1_t* record);
std::string GetCigarString(const vec<HtslibCigarOp>& cigar_ops);
std::string GetCigarString(const vec<HtslibCigarAlignInfo>& cigar_infos);
}  // namespace xoos::yc_decode
