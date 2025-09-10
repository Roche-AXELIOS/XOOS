#include "xoos/yc-decode/cigar-utils.h"

#include <fmt/format.h>

#include <xoos/util/string-functions.h>

namespace xoos::yc_decode {

/**
 * @brief Check whether the htslib CIGAR operator consumes read base(s)
 * @param op htslib CIGAR operator
 * @return true if the operator consumes read base(s), false otherwise
 */
bool ConsumesRead(const u32 op) {
  return op == BAM_CMATCH || op == BAM_CINS || op == BAM_CSOFT_CLIP || op == BAM_CEQUAL || op == BAM_CDIFF;
}

/**
 * @brief Check whether the htslib CIGAR operator consumes reference base(s)
 * @param op htslib CIGAR operator
 * @return true if the operator consumes reference base(s), false otherwise
 */
bool ConsumesReference(const u32 op) {
  return op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF;
}

/**
 * @brief Extract alignment information for each CIGAR operator from the BAM record.
 * @param record BAM record
 * @return Vector of CIGAR operator alignment information; hard-clips are not included,
 */
vec<HtslibCigarAlignInfo> GetCigarAlignInfos(const bam1_t* record) {
  vec<HtslibCigarAlignInfo> infos;
  const u32* cigar = bam_get_cigar(record);
  if (const auto num_ops = record->core.n_cigar; num_ops > 0) {
    u64 read_pos = 0;
    u64 ref_pos = record->core.pos;
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

/**
 * @brief Extract the string representation for htslib CIGAR operators.
 * @param cigar_ops Vector of CIGAR operators
 * @return CIGAR string
 */
std::string GetCigarString(const vec<HtslibCigarOp>& cigar_ops) {
  vec<std::string> cigar_str;
  cigar_str.reserve(cigar_ops.size());
  for (const auto& [op, len] : cigar_ops) {
    cigar_str.emplace_back(fmt::format("{}{}", len, BAM_CIGAR_STR[op]));
  }
  return string::Join(cigar_str, "");
}

/**
 * @brief Extract the string representation for htslib CIGAR operators.
 * @param cigar_infos Vector of CIGAR operator alignment information
 * @return CIGAR string
 */
std::string GetCigarString(const vec<HtslibCigarAlignInfo>& cigar_infos) {
  vec<std::string> cigar_str;
  cigar_str.reserve(cigar_infos.size());
  for (const auto& info : cigar_infos) {
    cigar_str.emplace_back(fmt::format("{}{}", std::max(info.read_span, info.ref_span), BAM_CIGAR_STR[info.op]));
  }
  return string::Join(cigar_str, "");
}

}  // namespace xoos::yc_decode
