#pragma once

#include <xoos/types/fs.h>
#include <xoos/types/str-container.h>

#include "variant-info.h"

namespace xoos::svc {

/**
 * @brief Protocols for normalizing feature values.
 * @details Supported protocols:
 * 1. "none": no normalization is applied to feature values.
 * 2. "median-dp": feature values are normalized by the median depth in VCF file.
 */
enum class FeatureNormalization {
  kNone,
  kMedianDp
};

/**
 * @brief Struct to hold depth of coverage values for normal and tumor samples.
 * @details The struct contains three fields:
 * 1. `normal`: median depth for the normal sample.
 * 2. `tumor`: median depth for the tumor sample.
 * 3. `total`: total median depth; sum of `normal` and `tumor`.
 * In the case of single-sample VCFs, `tumor` will be zero, but `normal` will still hold the median value.
 * Therefore, `total` will be equal to `normal` in single-sample scenarios.
 */
struct DepthTuple {
  u32 normal;
  u32 tumor;
  u32 total;
};

/**
 * @brief Struct to hold chromosome-level median depth for tumor and normal samples.
 * @details The struct contains two maps:
 * 1. `normal`: maps chromosome names to median depth values for the normal sample.
 * 2. `tumor`: maps chromosome names to median depth values for the tumor sample.
 * In the case of single-sample VCFs, only `normal` will be populated.
 */
struct ChromMedianDepth {
  StrUnorderedMap<u32> normal{};
  StrUnorderedMap<u32> tumor{};

  /**
   * @brief Get the DepthTuple for a given chromosome.
   * @post If the chromosome is not found, the depth values in the returned DpTuple will be zero.
   * @param chrom Chromosome name
   * @return DepthTuple containing normal, tumor, and total median depth for the chromosome
   */
  DepthTuple GetValue(const std::string& chrom) const;
};

/**
 * @brief Extract chromosome-level median depth from extracted VCF features.
 * @pre Assumes that VCF features were not filtered; there is one VCF feature for each variant in the VCF.
 * @post Both tumor and normal depth values are extracted from the VCF features. However, tumor depth values will be
 * zero if features are extracted from a single-sample VCF.
 * @param vcf_features Map of variant ID to VCF feature
 * @return Struct containing maps of chromosome to median depth
 */
ChromMedianDepth GetChromosomeMedianDepth(const VarIdToVcfFeatures& vcf_features);

/**
 * @brief Extract chromosome-level median depth from a VCF file.
 * @post For single-sample VCFs, depth values are extracted from the only sample as the normal sample.
 * @post For tumor-normal VCFs, depth values are extracted from both tumor and normal samples.
 * @param vcf_path Path of input VCF file
 * @return Struct containing maps of chromosome to median depth
 */
ChromMedianDepth GetChromosomeMedianDepth(const fs::path& vcf_path);

/**
 * @brief Get the sum of chromosome-level median depth values across all chromosomes.
 * @param chrom_median_depth Struct containing maps of chromosome to median depth
 * @return DepthTuple containing the sum of normal, tumor, and total depth values across all chromosomes
 */
DepthTuple GetSumOfChromosomeMedianDepth(const ChromMedianDepth& chrom_median_depth);

/**
 * @brief Get the median of chromosome-level median depth values across all chromosomes.
 * @param chrom_median_depth Struct containing maps of chromosome to median depth
 * @return DepthTuple containing the median normal, tumor, and total depth values across all chromosomes
 */
DepthTuple GetMedianOfChromosomeMedianDepth(const ChromMedianDepth& chrom_median_depth);

/**
 * @brief Validate the chromosome-level median depth values for feature normalization.
 * @details Checks that the normal sample has valid median depth values, and if there is sample context in the feature
 * columns, also checks that the tumor sample has valid median depth values. If any of the required median depth
 * values are missing or zero, an error is thrown because median depth normalization cannot be performed.
 * @param chr_med_dp Struct containing maps of chromosome to median depth
 * @param has_sample_context Boolean indicating whether there is sample context in the feature columns
 */
void ValidateChromMedianDepthForNormalization(const ChromMedianDepth& chr_med_dp, bool has_sample_context);

}  // namespace xoos::svc
