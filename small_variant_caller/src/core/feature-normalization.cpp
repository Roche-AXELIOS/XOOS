#include "feature-normalization.h"

#include <filesystem>
#include <ranges>

#include <xoos/types/vec.h>
#include <xoos/util/math.h>

#include "column-names.h"
#include "compute-vcf-features/vcf-header-util.h"
#include "vcf-fields.h"
#include "xoos/error/error.h"
#include "xoos/io/vcf/vcf-reader.h"
#include "xoos/log/logging.h"

namespace xoos::svc {

DepthTuple ChromMedianDepth::GetValue(const std::string& chrom) const {
  const auto tumor_itr = tumor.find(chrom);
  u32 tumor_val = 0;
  if (tumor_itr != tumor.end()) {
    tumor_val = tumor_itr->second;
  }
  const auto normal_itr = normal.find(chrom);
  u32 normal_val = 0;
  if (normal_itr != normal.end()) {
    normal_val = normal_itr->second;
  }
  return DepthTuple{.normal = normal_val, .tumor = tumor_val, .total = tumor_val + normal_val};
}

ChromMedianDepth GetChromosomeMedianDepth(const VarIdToVcfFeatures& vcf_features) {
  ChromMedianDepth result{};
  StrUnorderedMap<vec<u32>> tumor_chrom_to_dp_vals{};
  StrUnorderedMap<vec<u32>> normal_chrom_to_dp_vals{};
  StrUnorderedMap<std::unordered_set<u64>> seen_positions{};
  for (const auto& [vid, feat] : vcf_features) {
    if (!seen_positions[vid.chrom].insert(vid.pos).second) {
      // variant at the same position shares the same DP value, skip it
      continue;
    }
    tumor_chrom_to_dp_vals[vid.chrom].emplace_back(feat.tumor_dp);
    normal_chrom_to_dp_vals[vid.chrom].emplace_back(feat.normal_dp);
  }

  // Calculate median DP for each chromosome for the tumor sample
  for (auto& [chrom, dp_vals] : tumor_chrom_to_dp_vals) {
    if (!dp_vals.empty()) {
      result.tumor[chrom] = math::Median(dp_vals);
    }
  }

  // Calculate median DP for each chromosome for the normal sample
  for (auto& [chrom, dp_vals] : normal_chrom_to_dp_vals) {
    if (!dp_vals.empty()) {
      result.normal[chrom] = math::Median(dp_vals);
    }
  }

  return result;
}

ChromMedianDepth GetChromosomeMedianDepth(const fs::path& vcf_path) {
  io::VcfReader vcf_reader(vcf_path);
  const auto hdr = vcf_reader.GetHeader();
  if (!hdr->HasInfoField(kFieldDp)) {
    // DP field not found in the VCF header, return empty result
    return {};
  }

  s32 tumor_idx = -1;
  s32 normal_idx = -1;
  // Get tumor_sample and normal_sample indexes from the VCF header
  const auto tn_sample_idx = GetTumorNormalSampleIndexes(hdr);
  if (tn_sample_idx.has_value()) {
    // tumor/normal sample indexes not found in the VCF header, return empty result
    tumor_idx = tn_sample_idx->tumor_sample_idx;
    normal_idx = tn_sample_idx->normal_sample_idx;
  } else {
    // assume single-sample VCF, use sample index 0 as normal sample
    normal_idx = 0;
  }

  // Iterate through the VCF records and collect DP values for each chromosome and sample
  ChromMedianDepth results;
  StrUnorderedMap<vec<u32>> tumor_chr_to_dp_vec{};
  StrUnorderedMap<vec<u32>> normal_chr_to_dp_vec{};
  StrUnorderedMap<std::unordered_set<hts_pos_t>> chr_to_pos_set{};
  while (const auto& record = vcf_reader.GetNextRecord()) {
    const auto chrom = record->Chromosome();
    if (!chr_to_pos_set[chrom].insert(record->Position()).second) {
      // record at the same position shares the same DP value, skip it
      continue;
    }
    const auto& values = record->GetFormatFieldNoCheck<s32>(kFieldDp);
    if (0 <= tumor_idx && std::cmp_less(tumor_idx, values.size())) {
      tumor_chr_to_dp_vec[chrom].emplace_back(values.at(tumor_idx));
    }
    if (0 <= normal_idx && std::cmp_less(normal_idx, values.size())) {
      normal_chr_to_dp_vec[chrom].emplace_back(values.at(normal_idx));
    }
  }

  // Calculate the median DP for each chromosome in the tumor sample
  for (auto& [chrom, dp_vals] : tumor_chr_to_dp_vec) {
    if (!dp_vals.empty()) {
      results.tumor[chrom] = math::Median(dp_vals);
    }
  }

  // Calculate the median DP for each chromosome in the normal sample
  for (auto& [chrom, dp_vals] : normal_chr_to_dp_vec) {
    if (!dp_vals.empty()) {
      results.normal[chrom] = math::Median(dp_vals);
    }
  }

  return results;
}

DepthTuple GetSumOfChromosomeMedianDepth(const ChromMedianDepth& chrom_median_depth) {
  u32 tumor_dp = 0;
  u32 normal_dp = 0;
  for (const auto& dp : chrom_median_depth.tumor | std::views::values) {
    if (dp > 0) {
      tumor_dp += dp;
    }
  }
  for (const auto& dp : chrom_median_depth.normal | std::views::values) {
    if (dp > 0) {
      normal_dp += dp;
    }
  }
  return DepthTuple{.normal = normal_dp, .tumor = tumor_dp, .total = tumor_dp + normal_dp};
}

DepthTuple GetMedianOfChromosomeMedianDepth(const ChromMedianDepth& chrom_median_depth) {
  vec<u32> tumor_dps{};
  vec<u32> normal_dps{};
  for (const auto& dp : chrom_median_depth.tumor | std::views::values) {
    if (dp > 0) {
      tumor_dps.emplace_back(dp);
    }
  }
  for (const auto& dp : chrom_median_depth.normal | std::views::values) {
    if (dp > 0) {
      normal_dps.emplace_back(dp);
    }
  }
  DepthTuple result{};
  if (!tumor_dps.empty()) {
    result.tumor = math::Median(tumor_dps);
  }
  if (!normal_dps.empty()) {
    result.normal = math::Median(normal_dps);
  }
  result.total = result.tumor + result.normal;
  return result;
}

void ValidateChromMedianDepthForNormalization(const ChromMedianDepth& chr_med_dp, const bool has_sample_context) {
  // Calculate the sum for all chromosomes.
  // If the sum of median depth values across all chromosomes is zero for the normal sample, then it means that all
  // median depth values for the normal sample are zero, which is invalid for normalization. An error is thrown in this
  // case because median depth normalization cannot be performed with invalid median depth values.
  // If there is sample context in the feature columns, then the same check is also performed for the tumor sample
  // because tumor sample median depth values are also required for normalization in this case.
  const auto& [normal_sum, tumor_sum, total_sum] = GetSumOfChromosomeMedianDepth(chr_med_dp);
  if (normal_sum == 0) {
    // "normal_dp" is always required, regardless of sample context of feature columns.
    const std::string sample = has_sample_context ? "normal sample" : "sample";
    throw error::Error(
        "The chromosome median depth values for the {} are all zero. Feature normalization cannot be "
        "performed. Please check your VCF files for non-zero FORMAT field '{}' and the configuration file for '{}' in "
        "the VCF features list.",
        sample,
        kFieldDp,
        kNameNormalDP);
  }
  if (has_sample_context && tumor_sum == 0) {
    // If there is sample context in the feature columns, then "tumor_dp" is required for median depth
    // calculation for the tumor sample
    throw error::Error(
        "The chromosome median depth values for the tumor sample are all zero. Feature normalization cannot be "
        "performed. Please check your VCF files for non-zero FORMAT field '{}' and the configuration file for '{}' in "
        "the VCF features list.",
        kFieldDp,
        kNameTumorDP);
  }
  // log the median of chromosome median depth values
  const auto& [normal_m, tumor_m, total_m] = GetMedianOfChromosomeMedianDepth(chr_med_dp);
  if (has_sample_context) {
    Logging::Info("Median of chromosome median depths: {} (Normal: {}, Tumor: {})", total_m, normal_m, tumor_m);
  } else {
    Logging::Info("Median of chromosome median depths: {}", total_m);
  }
}

}  // namespace xoos::svc
