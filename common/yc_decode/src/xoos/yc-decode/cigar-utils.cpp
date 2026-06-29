#include "xoos/yc-decode/cigar-utils.h"

#include <fmt/format.h>

#include <xoos/util/string-functions.h>

namespace xoos::yc_decode {

bool ConsumesRead(const u32 op) {
  static constexpr u32 kQueryBitMask = 0b01;
  // CIGAR ops: M, I, S, =, X
  return (bam_cigar_type(op) & kQueryBitMask) != 0;
}

bool ConsumesReference(const u32 op) {
  static constexpr u32 kReferenceBitMask = 0b10;
  // CIGAR ops: M, D, N, =, X
  return (bam_cigar_type(op) & kReferenceBitMask) != 0;
}

u32 GetReadLength(const vec<HtslibCigarOp>& op) {
  u32 len = 0;
  for (const auto& [op_type, op_len] : op) {
    if (ConsumesRead(op_type)) {
      len += op_len;
    }
  }
  return len;
}

vec<HtslibCigarAlignInfo> GetCigarAlignInfos(const bam1_t* const record) {
  vec<HtslibCigarAlignInfo> infos;
  const u32* const cigar = bam_get_cigar(record);
  if (const auto num_ops = record->core.n_cigar; num_ops > 0) {
    u64 read_pos = 0;
    auto ref_pos = ToUnsigned(record->core.pos);
    for (u32 i = 0; i < num_ops; ++i) {
      auto op = bam_cigar_op(cigar[i]);
      if (op != BAM_CHARD_CLIP) {
        const auto op_len = bam_cigar_oplen(cigar[i]);
        const u64 read_span = ConsumesRead(op) ? op_len : 0;
        const u64 ref_span = ConsumesReference(op) ? op_len : 0;
        infos.emplace_back(op, read_pos, ref_pos, read_span, ref_span);
        read_pos += read_span;
        ref_pos += ref_span;
      }
    }
  }
  return infos;
}

std::string GetCigarString(const vec<HtslibCigarOp>& cigar_ops) {
  vec<std::string> cigar_str;
  cigar_str.reserve(cigar_ops.size());
  for (const auto& [op, len] : cigar_ops) {
    cigar_str.emplace_back(fmt::format("{}{}", len, BAM_CIGAR_STR[op]));
  }
  return string::Join(cigar_str, "");
}

std::string GetCigarString(const vec<HtslibCigarAlignInfo>& cigar_infos) {
  vec<std::string> cigar_str;
  cigar_str.reserve(cigar_infos.size());
  for (const auto& info : cigar_infos) {
    cigar_str.emplace_back(fmt::format("{}{}", std::max(info.read_span, info.ref_span), BAM_CIGAR_STR[info.op]));
  }
  return string::Join(cigar_str, "");
}

}  // namespace xoos::yc_decode
