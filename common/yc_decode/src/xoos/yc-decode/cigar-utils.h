#pragma once

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::yc_decode {

/**
 * @brief Basic information of a single htslib CIGAR operation.
 */
struct HtslibCigarOp {
  // htslib CIGAR op, e.g. `BAM_CMATCH`
  u8 op{};
  // length of the CIGAR op
  u32 len{0};
  auto operator<=>(const HtslibCigarOp&) const = default;
};

/**
 * @brief Alignment information of a single htslib CIGAR operation.
 */
struct HtslibCigarAlignInfo {
  // htslib CIGAR op, e.g. `BAM_CMATCH`
  u8 op{};
  // 0-based start position in the read
  u64 read_pos{0};
  // 0-based start position in the reference
  u64 ref_pos{0};
  // number of read bases consumed by the CIGAR op
  u64 read_span{0};
  // number of reference bases consumed by the CIGAR op
  u64 ref_span{0};
  auto operator<=>(const HtslibCigarAlignInfo&) const = default;
};

/**
 * @brief Check whether the htslib CIGAR operator consumes read base(s)
 * @param op htslib CIGAR operator
 * @return true if the operator consumes read base(s), false otherwise
 */
bool ConsumesRead(u32 op);

/**
 * @brief Check whether the htslib CIGAR operator consumes reference base(s)
 * @param op htslib CIGAR operator
 * @return true if the operator consumes reference base(s), false otherwise
 */
bool ConsumesReference(u32 op);

/**
 * @brief Get the read length from the vector of htslib CIGAR operators.
 * @param op CIGAR operators
 * @return read length
 */
u32 GetReadLength(const vec<HtslibCigarOp>& op);

/**
 * @brief Extract alignment information for each CIGAR operator from the BAM record.
 * @param record BAM record
 * @return Vector of CIGAR operator alignment information; hard-clips are not included,
 */
vec<HtslibCigarAlignInfo> GetCigarAlignInfos(const bam1_t* record);

/**
 * @brief Extract the string representation for htslib CIGAR operators.
 * @param cigar_ops Vector of CIGAR operators
 * @return CIGAR string
 */
std::string GetCigarString(const vec<HtslibCigarOp>& cigar_ops);

/**
 * @brief Extract the string representation for htslib CIGAR operators.
 * @param cigar_infos Vector of CIGAR operator alignment information
 * @return CIGAR string
 */
std::string GetCigarString(const vec<HtslibCigarAlignInfo>& cigar_infos);
}  // namespace xoos::yc_decode
