#pragma once
#include <optional>

#include "xoos/io/vcf/vcf-header.h"

namespace xoos::svc {

/**
 * @brief Structure to hold the indexes of tumor and normal samples in a VCF file.
 */
struct TumorNormalSampleIndexes {
  s32 tumor_sample_idx{-1};
  s32 normal_sample_idx{-1};
};

/**
 * @brief Extract the tumor and normal sample names from the VCF header and return their corresponding sample IDs
 * as a pair.
 *
 * @return An optional pair of integers representing the tumor and normal sample IDs in the VCF file. If either
 * sample name is not found, or their IDs are the same or invalid, then std::nullopt is returned.
 * @note The function looks for the "tumor_sample" and "normal_sample" keys in the VCF header.
 * @note Sample IDs are 1-based
 */
std::optional<TumorNormalSampleIndexes> GetTumorNormalSampleIndexes(const io::VcfHeaderPtr& hdr);

}  // namespace xoos::svc
