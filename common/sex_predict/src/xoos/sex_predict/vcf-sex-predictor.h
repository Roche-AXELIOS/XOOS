#pragma once

#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>

#include "sex.h"

namespace xoos::sex_predict {

/**
 * @brief Extract chromosome-level median FORMAT field DP from VCF file.
 * @param vcf_path Path to a VCF file
 * @return Map of chromosome to median DP
 */
StrUnorderedMap<u32> GetChromosomeMedianDP(const fs::path& vcf_path);

/**
 * @brief Extract chromosome-level median FORMAT field DP from VCF file for specified autosome and chrX only.
 * @param vcf_path Path to a VCF file
 * @param autosome_name Autosomal chromosome name
 * @param chr_x_name Chromosome X name
 * @return Map of chromosome to median DP
 */
StrUnorderedMap<u32> GetChromosomeMedianDP(const fs::path& vcf_path,
                                           const std::string& autosome_name,
                                           const std::string& chr_x_name);

/**
 * @brief Predict sex based on the ratio of median DP between an autosome and chromosome X.
 * @details A normal human diploid genome has 2 copies of autosomes (e.g. chr1).
 * A biological female has 2 copies of chrX; her expected median DP ratio of chr1 to chrX is ~1.0.
 * A biological male has 1 copy of chrX and 1 copy of chrY; his expected median DP ratio of chr1 to chrX is ~2.0.
 * @param autosome_dp Median DP of an autosomal chromosome
 * @param chr_x_dp Median DP of chromosome X
 * @return Predicted sex
 */
Sex PredictSex(u32 autosome_dp, u32 chr_x_dp);

/**
 * @brief Predict sex from chromosome-level median DP.
 * @param chr_to_dp Map of chromosome to median DP
 * @param autosome_name Autosomal chromosome name
 * @param chr_x_name Chromosome X name
 * @return Predicted sex
 */
Sex PredictSex(const StrUnorderedMap<u32>& chr_to_dp, const std::string& autosome_name, const std::string& chr_x_name);

/**
 * @brief Predict sex from VCF file.
 * @param vcf_path Path to a VCF file
 * @param autosome_name Autosomal chromosome name
 * @param chr_x_name Chromosome X name
 * @return Predicted sex
 */
Sex PredictSex(const fs::path& vcf_path, const std::string& autosome_name, const std::string& chr_x_name);

}  // namespace xoos::sex_predict
