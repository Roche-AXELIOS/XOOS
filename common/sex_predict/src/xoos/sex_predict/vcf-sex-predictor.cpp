#include "vcf-sex-predictor.h"

#include <xoos/error/error.h>
#include <xoos/io/vcf/vcf-reader.h>

#include "xoos/types/vec.h"
#include "xoos/util/math.h"

namespace xoos::sex_predict {

using enum Sex;

// VCF FORMAT field name for DP
static const std::string kFieldDp = "DP";

StrUnorderedMap<u32> GetChromosomeMedianDP(const fs::path& vcf_path) {
  const io::VcfReader vcf_reader(vcf_path);
  if (!vcf_reader.GetHeader()->HasInfoField(kFieldDp)) {
    // If the VCF does not have DP field in the header, return an empty map
    return {};
  }

  // Iterate through the VCF records and collect DP values for each chromosome
  StrUnorderedMap<vec<u32>> chr_to_dp_vec{};
  StrUnorderedMap<std::set<hts_pos_t>> chr_to_pos_set{};
  while (const auto& record = vcf_reader.GetNextRecord()) {
    const auto chrom = record->Chromosome();
    if (!chr_to_pos_set[chrom].insert(record->Position()).second) {
      // record at the same position shares the same DP value, skip it
      continue;
    }
    const auto& values = record->GetFormatFieldNoCheck<s32>(kFieldDp);
    if (!values.empty()) {
      chr_to_dp_vec[chrom].emplace_back(values.front());
    }
  }

  // Calculate the median DP for each chromosome
  StrUnorderedMap<u32> chr_to_dp{};
  for (auto& [chrom, dp_vals] : chr_to_dp_vec) {
    if (!dp_vals.empty()) {
      chr_to_dp[chrom] = math::Median(dp_vals);
    }
  }
  return chr_to_dp;
}

/**
 * @brief Helper function to extract median DP from VCF reader for the set region.
 * @param reader VCF reader with region set
 * @return Median DP value for the set region
 */
static u32 GetMedianDP(io::VcfReader& reader) {
  vec<u32> dps{};
  std::set<hts_pos_t> pos_set{};
  // Note that the region must be already set in the VCF reader
  // Use `GetNextRegionRecord` (instead of `GetNextRecord`) to iterate through records in the set region
  while (const auto& record = reader.GetNextRegionRecord(BCF_UN_FMT)) {
    if (!pos_set.insert(record->Position()).second) {
      // record at the same position shares the same DP value, skip it
      continue;
    }
    const auto& values = record->GetFormatFieldNoCheck<s32>(kFieldDp);
    if (!values.empty()) {
      dps.emplace_back(values.front());
    }
  }
  return math::Median(dps);
}

StrUnorderedMap<u32> GetChromosomeMedianDP(const fs::path& vcf_path,
                                           const std::string& autosome_name,
                                           const std::string& chr_x_name) {
  if (autosome_name.empty() || chr_x_name.empty()) {
    return {};
  }
  io::VcfReader vcf_reader(vcf_path);
  if (!vcf_reader.GetHeader()->HasInfoField(kFieldDp)) {
    return {};
  }

  StrUnorderedMap<u32> chr_to_dp{};
  if (vcf_reader.HasIndex()) {
    const auto ctg_lengths = vcf_reader.GetHeader()->GetContigLengths();
    const auto autosome_len_itr = ctg_lengths.find(autosome_name);
    const auto chr_x_len_itr = ctg_lengths.find(chr_x_name);
    if (autosome_len_itr != ctg_lengths.end() && chr_x_len_itr != ctg_lengths.end()) {
      // If chromosome lengths are found, extract DP values for autosomes and chrX only.
      if (vcf_reader.SetRegion(autosome_name, 0, autosome_len_itr->second)) {
        chr_to_dp[autosome_name] = GetMedianDP(vcf_reader);
      }
      if (vcf_reader.SetRegion(chr_x_name, 0, chr_x_len_itr->second)) {
        chr_to_dp[chr_x_name] = GetMedianDP(vcf_reader);
      }
      return chr_to_dp;
    }
  }

  // Otherwise, extract median DP values for all chromosomes and return the values for autosome and chrX only.
  const auto all_chrom_to_dp = GetChromosomeMedianDP(vcf_path);
  const auto autosome_itr = all_chrom_to_dp.find(autosome_name);
  if (autosome_itr != all_chrom_to_dp.end()) {
    chr_to_dp[autosome_name] = autosome_itr->second;
  }
  const auto chr_x_itr = all_chrom_to_dp.find(chr_x_name);
  if (chr_x_itr != all_chrom_to_dp.end()) {
    chr_to_dp[chr_x_name] = chr_x_itr->second;
  }
  return chr_to_dp;
}

Sex PredictSex(const u32 autosome_dp, const u32 chr_x_dp) {
  if (chr_x_dp == 0) {
    // Avoid division by zero
    return kUnknown;
  }
  // Ratio <= 0.5: chrX DP is much higher than autosome DP
  static constexpr f64 kLowRatioThreshold = 0.5;
  // Ratio >= 2.5: chrX DP is much lower than autosome DP
  static constexpr f64 kHighRatioThreshold = 2.5;
  const f64 ratio = autosome_dp / static_cast<f64>(chr_x_dp);
  if (ratio <= kLowRatioThreshold || ratio >= kHighRatioThreshold) {
    return kUnknown;
  }
  static constexpr s64 kFemaleRatio = 1;
  static constexpr s64 kMaleRatio = 2;
  // Round the ratio to the nearest integer
  const auto rounded_ratio = lround(ratio);
  if (rounded_ratio == kFemaleRatio) {
    // Autosome and chrX have similar DP
    return kFemale;
  }
  if (rounded_ratio == kMaleRatio) {
    // Autosome DP is about double of chrX DP
    return kMale;
  }
  return kUnknown;
}

Sex PredictSex(const StrUnorderedMap<u32>& chr_to_dp, const std::string& autosome_name, const std::string& chr_x_name) {
  const auto autosome_itr = chr_to_dp.find(autosome_name);
  const auto chr_x_itr = chr_to_dp.find(chr_x_name);
  if (autosome_itr == chr_to_dp.end() || chr_x_itr == chr_to_dp.end()) {
    return kUnknown;
  }
  return PredictSex(autosome_itr->second, chr_x_itr->second);
}

Sex PredictSex(const fs::path& vcf_path, const std::string& autosome_name, const std::string& chr_x_name) {
  const auto chr_to_dp = GetChromosomeMedianDP(vcf_path, autosome_name, chr_x_name);
  if (chr_to_dp.empty()) {
    return kUnknown;
  }
  return PredictSex(chr_to_dp, autosome_name, chr_x_name);
}

}  // namespace xoos::sex_predict
