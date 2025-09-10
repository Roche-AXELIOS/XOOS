#include "filtering.h"

#include <algorithm>
#include <ranges>

#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>

#include "util/seq-util.h"

namespace xoos::svc {
constexpr int kMaxSubstId = 13;
constexpr int kIndelSubstId = 11;
constexpr double kFpRateIndel = 0.005;
constexpr double kFpRateSnv = 0.0005;

/**
 * @brief Filters a given variant based on the provided filters and thresholds
 * @param vid The VariantInfo struct containing info about the variant
 * @param uvf The UnifiedVariantFeature struct containing features of the variant
 * @param urf The UnifiedReferenceFeature struct containing features of the reference
 * @param settings The FilterSettings struct containing relevant filters and thresholds to filter the variant against
 * @return A string representing the failure reason or PASS
 */
std::vector<std::string> FilterVariant(const VariantId& vid,
                                       const UnifiedVariantFeature& uvf,
                                       const UnifiedReferenceFeature& urf,
                                       const FilterSettings& settings) {
  std::vector<std::string> fail_reasons;
  // TODO: Update the way filtering is done based on workflow/User profiles in future change
  // TODO: value == threshold should not fail? Use `<` instead of `<=`?
  if (uvf.mapq_mean < settings.min_mapq) {
    fail_reasons.emplace_back(kFilteringMapQualityId);
  }
  if (uvf.baseq_mean < settings.min_baseq) {
    fail_reasons.emplace_back(kFilteringBaseQualityId);
  }
  const std::string vardesc = GetVariantCorrelationKey(vid, false);

  auto is_hotspot = settings.hotspots.contains(vardesc);

  if (is_hotspot) {
    if (uvf.weighted_score < settings.hotspot_weighted_counts_threshold) {
      fail_reasons.emplace_back(kFilteringCountsId);
    }
    if (uvf.ml_score < settings.hotspot_ml_threshold) {
      fail_reasons.emplace_back(kFilteringMLScoreId);
    }
  } else {
    auto threshold = settings.weighted_counts_thresholds.at(SubstIndex(vid.ref, vid.alt));
    if (uvf.weighted_score < threshold) {
      fail_reasons.emplace_back(kFilteringCountsId);
    }
    if (uvf.ml_score < settings.ml_threshold) {
      fail_reasons.emplace_back(kFilteringMLScoreId);
    }
  }
  if (settings.blocklist.contains(vardesc)) {
    // Variant is blocklisted
    fail_reasons.emplace_back(kFilteringBlocklistedId);
  }
  if (static_cast<float>(uvf.duplex + uvf.nonduplex) / static_cast<float>(urf.support + uvf.duplex + uvf.nonduplex) <
      settings.min_af_threshold) {
    // TODO: check whether the AF denominator should be tallied for all ALT alleles at this position?
    // AF is too low
    fail_reasons.emplace_back(kFilteringAFId);
  }
  if (fail_reasons.empty()) {
    fail_reasons.emplace_back(kFilteringPassId);
  }
  return fail_reasons;
}

/**
 * Filters a phased variant call based on a set list of checks
 * @param key A left-padded variant correlation key used to query the blocklist
 * @param settings a PhasedFilterSettings struct that contains relevant filters and thesholds to filter against
 * @param allele_depth The allele depth of the variant
 * @param allele_freq The allele frequency of the variant
 * @param mapping_qual The mapping quality of the variant
 * @param base_qual The base quality of the variant
 * @return A string representing the failure reason or PASS
 */
std::vector<std::string> FilterPhasedVariant(const std::string& key,
                                             const PhasedFilterSettings& settings,
                                             const u32 allele_depth,
                                             const float allele_freq,
                                             const u8 mapping_qual,
                                             const u8 base_qual) {
  std::vector<std::string> fail_reasons;
  // TODO: Confirm if blocklist filtering needed here
  if (settings.blocklist.contains(key)) {
    fail_reasons.emplace_back(kFilteringBlocklistedId);
  }
  if (allele_depth < settings.min_alt_counts) {
    fail_reasons.emplace_back(kFilteringMinAltCountsId);
  }
  if (base_qual < settings.min_baseq) {
    fail_reasons.emplace_back(kFilteringBaseQualityId);
  }
  if (mapping_qual < settings.min_mapq) {
    fail_reasons.emplace_back(kFilteringMapQualityId);
  }
  if (allele_freq < settings.min_allele_frequency || allele_freq > settings.max_allele_frequency) {
    fail_reasons.emplace_back(kFilteringAFId);
  }
  if (fail_reasons.empty()) {
    fail_reasons.emplace_back(kFilteringPassId);
  }
  return fail_reasons;
}

// Helper function to convert an integer to a string and pads it with leading zeros to a specified total width
template <typename T>
static std::string ConvertIntToStringLeftPadWith0(T value, int total_width = 10) {
  std::string result = std::to_string(value);
  const int zeros_to_add = total_width - static_cast<int>(result.size());
  if (zeros_to_add > 0) {
    return std::string(zeros_to_add, '0') + result;
  }
  return result;
}

/**
 * @brief Generate a string for a variant.
 * @param chrom Chromosome name
 * @param position Position on the chromosome
 * @param ref Reference allele
 * @param alt Alternate allele
 * @param pad_left Flag to indicate whether to pad the position with leading zeros
 * @return String representation for the variant
 */
std::string GetVariantCorrelationKey(
    const std::string& chrom, u64 position, const std::string& ref, const std::string& alt, const bool pad_left) {
  return fmt::format(
      "{}_{}_{}_{}", chrom, (pad_left ? ConvertIntToStringLeftPadWith0(position) : std::to_string(position)), ref, alt);
}

/**
 * @brief Generate a string for a variant.
 * @param vid VariantId struct for the variant
 * @param pad_left Flag to indicate whether to pad the position with leading zeros
 * @return String representation for the variant
 */
std::string GetVariantCorrelationKey(const VariantId& vid, const bool pad_left) {
  return GetVariantCorrelationKey(vid.chrom, vid.pos, vid.ref, vid.alt, pad_left);
}

/**
 * @brief Extract the string set of hotspot variants from a VCF file.
 * @param vcf_path Path to the VCF file containing hotspot variants
 * @return A set of strings representing the hotspot variants
 */
StrUnorderedSet LoadHotspotVariants(const fs::path& vcf_path) {
  const io::VcfReader reader(vcf_path);
  io::VcfRecordPtr record;
  StrUnorderedSet hotspots;
  while ((record = reader.GetNextRecord())) {
    auto prefix = fmt::format("{}_{}_{}_", record->Chromosome(), record->Position(), record->Allele(0));
    for (auto allele = 1; allele < record->NumAlleles(); allele++) {
      hotspots.insert(prefix + record->Allele(allele));
    }
  }
  return hotspots;
}

/**
 * @brief Retain features for SNVs in the specified VCF.
 * @param vcf Path to the VCF file containing SNVs to be retained
 * @param features Features to be filtered, keyed by chromosome and position
 * @return A map of chromosome to position to variant info containing only the truth variants
 */
ChromToVariantInfoMap FilterVariantsByVcf(const fs::path& vcf, const ChromToVariantInfoMap& features) {
  // Extract SNVs from the VCF file and store them in a set for quick lookup
  StrUnorderedSet snv_set;
  // TODO: support indels?
  const io::VcfReader reader(vcf);
  while (const auto& vcf_record = reader.GetNextRecord()) {
    const auto& ref = vcf_record->Allele(0);
    if (ref.size() == 1 && !IsAnyNotACTG(ref)) {
      for (int i = 1; i < vcf_record->NumAlleles(); ++i) {
        const auto& alt = vcf_record->Allele(i);
        if (alt.size() == 1 && !IsAnyNotACTG(alt)) {
          snv_set.insert(GetVariantCorrelationKey(vcf_record->Chromosome(), vcf_record->Position(), ref, alt, false));
        }
      }
    }
  }
  // Filter the features map to retain only those variants that are in the SNV set
  ChromToVariantInfoMap positives;
  for (const auto& [chr, pos_data] : features) {
    for (const auto& [pos, vid_data] : pos_data) {
      for (const auto& [vid, var_feat] : vid_data) {
        if (snv_set.contains(GetVariantCorrelationKey(vid, false))) {
          positives[chr][pos][vid] = var_feat;
        }
      }
    }
  }
  return positives;
}

/**
 * @brief Calculate weighted count thresholds for each substitution type based on the provided features and panel size.
 * @param features Map of chromosome to position to variant info
 * @param panel_size Size of the panel used for calculating thresholds
 * @param default_threshold Default threshold value to use if no variants are found
 * @return A map of substitution type ID to calculated threshold value
 */
std::unordered_map<int, double> CalculateWeightedCountThresholdsPerSubstitutionType(
    const ChromToVariantInfoMap& features, const u32 panel_size, const double default_threshold) {
  std::unordered_map<int, double> thresholds;
  if (0 == panel_size) {
    for (auto subst_id = 1; subst_id <= kMaxSubstId; subst_id++) {
      thresholds[subst_id] = default_threshold;
    }
    return thresholds;
  }

  // Calculate weighted count thresholds
  for (auto subst_id = 1; subst_id <= kMaxSubstId; subst_id++) {
    thresholds[subst_id] = default_threshold;
    // extract the weighted scores for the current substitution type
    std::vector<double> weighted_scores;
    for (const auto& pos_data : std::views::values(features)) {
      for (const auto& vid_data : std::views::values(pos_data)) {
        for (const auto& [vid, feat] : vid_data) {
          if (SubstIndex(vid.ref, vid.alt) == subst_id) {
            weighted_scores.push_back(feat.weighted_score);
          }
        }
      }
    }
    if (weighted_scores.empty()) {
      // use the default threshold
      continue;
    }

    // Count the number of variants with weighted scores greater than or equal to each possible integer score
    // Normalize the counts by the panel size
    std::vector<double> scaled_score_counts;
    const int max_score = static_cast<int>(*std::ranges::max_element(weighted_scores));
    for (int i = 1; i <= max_score + 1; i++) {
      scaled_score_counts.push_back(
          static_cast<double>(std::ranges::count_if(weighted_scores, [i](double v) { return v >= i; })) /
          static_cast<double>(panel_size));
    }

    // Find the index of the scaled score count that is closest to the false positive rate
    const double diff_amount = subst_id == kIndelSubstId ? kFpRateIndel : kFpRateSnv;
    double min_diff = diff_amount;
    size_t min_diff_index = 0;
    for (size_t i = 0; i < scaled_score_counts.size(); i++) {
      const auto projectedfps = scaled_score_counts[i];
      auto diff = std::abs(projectedfps - diff_amount);
      if (diff < min_diff) {
        min_diff = diff;
        min_diff_index = i;
      }
    }
    // min_diff_index is the index of scaled_score_counts that is closest to diff_amount
    const double false_positive_threshold = std::max(default_threshold, static_cast<double>(min_diff_index + 1));
    thresholds[subst_id] = false_positive_threshold;
  }

  return thresholds;
}

}  // namespace xoos::svc
