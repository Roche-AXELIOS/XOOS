#pragma once

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>

#include "core/variant-info.h"

namespace xoos::svc {

// VCF outout FILTER field IDs
const std::string kFilteringMapQualityId = "mapquality";
const std::string kFilteringBaseQualityId = "basequality";
const std::string kFilteringAFId = "Allele_Freq";
const std::string kFilteringMinAltCountsId = "Min_Alt_Counts";
const std::string kFilteringCountsId = "Counts";
const std::string kFilteringMLScoreId = "ML_Score";
const std::string kFilteringBlocklistedId = "BLOCKLISTED";
const std::string kFilteringForcedId = "FORCED";
const std::string kFilteringPassId = "PASS";
const std::string kFilteringFailId = "FAIL";
const std::string kFilteringFailSomaticTNId = "FAIL_Somatic_TN_ML";
const std::string kFilteringFailSomaticTNGermlineId = "FAIL_Somatic_TN_Germline";
const std::string kFilteringMissingFeatureId = "missing_feature";
const std::string kFilteringFalsePositiveId = "false_positive";
const std::string kFilteringMultialleleFormatId = "multiallele_format";
const std::string kFilteringMultiallelePartnerId = "multiallele_partner";
const std::string kFilteringMultialleleConflictId = "multiallele_conflict";
const std::string kFilteringNonAcgtRefAltId = "non_acgt_ref_alt";

// VCF output descriptions
const std::string kFilteringPassDesc = "All filters passed";
const std::string kFilteringGermlinePassDesc = "Site contains at least one allele that passes filters";
const std::string kFilteringGermlineFailDesc = "Variant is filtered out by ML";
const std::string kFilteringMapQualityDesc = "Filtered due to mapping quality";
const std::string kFilteringAFDesc = "Filtered due to allele frequency";
const std::string kFilteringMinAltCountsDesc = "Filtered due to alt counts";
const std::string kFilteringBaseQualityDesc = "Filtered due to base quality";
const std::string kFilteringCountsDesc = "Filtered due to low counts";
const std::string kFilteringMLScoreDesc = "Filtered due to low ML Score";
const std::string kFilteringForcedDesc = "Filtered due to being a forced variant call";
const std::string kFilteringBlocklistedDesc = "Variant blocklisted";
const std::string kFilteringFailSomaticTNDesc = "Filtered due to low Somatic TN ML Score";
const std::string kFilteringFailSomaticTNGermlineDesc =
    "Filtered due to being classified as a possible germline variant";
const std::string kFilteringMissingFeatureDesc = "Filtered due to missing ML feature";
const std::string kFilteringFalsePositiveDesc = "Filtered due to ML model classification as false positive";
const std::string kFilteringMultialleleFormatDesc = "Filtered due to multi-allele record format failure";
const std::string kFilteringMultiallelePartnerDesc = "Filtered due to missing multi-allele partner";
const std::string kFilteringMultialleleConflictDesc =
    "Filtered due to conflicting multi-allele ML model classification";
const std::string kFilteringNonAcgtRefAltDesc = "Filtered due to non-ACGT reference or alternate allele(s)";

// Settings for `somatic` workflow variants filtering
struct FilterSettings {
  FilterSettings(u8 min_mapq,
                 u8 min_baseq,
                 const std::unordered_map<int, double>& weighted_counts_thresholds,
                 float ml_threshold,
                 StrUnorderedSet blocklist,
                 float hotspot_weighted_counts_threshold,
                 float hotspot_ml_threshold,
                 float min_af_threshold,
                 StrUnorderedSet hotspots)
      : min_mapq(min_mapq),
        min_baseq(min_baseq),
        weighted_counts_thresholds(weighted_counts_thresholds),
        ml_threshold(ml_threshold),
        blocklist(std::move(blocklist)),
        hotspot_weighted_counts_threshold(hotspot_weighted_counts_threshold),
        hotspot_ml_threshold(hotspot_ml_threshold),
        min_af_threshold(min_af_threshold),
        hotspots(std::move(hotspots)) {
  }

  const u8 min_mapq{0};
  const u8 min_baseq{0};
  const std::unordered_map<int, double> weighted_counts_thresholds{};
  const float ml_threshold{0};
  const StrUnorderedSet blocklist{};
  const float hotspot_weighted_counts_threshold{0};
  const float hotspot_ml_threshold{0};
  const float min_af_threshold{0};
  const StrUnorderedSet hotspots{};
};

// Settings for `somatic` workflow phased variants filtering
struct PhasedFilterSettings {
  PhasedFilterSettings(
      u8 min_mapq, u8 min_baseq, StrUnorderedSet blocklist, float min_af, float max_af, u32 min_alt_counts)
      : min_mapq(min_mapq),
        min_baseq(min_baseq),
        blocklist(std::move(blocklist)),
        min_allele_frequency(min_af),
        max_allele_frequency(max_af),
        min_alt_counts(min_alt_counts) {
  }

  const u8 min_mapq{0};
  const u8 min_baseq{0};
  const StrUnorderedSet blocklist{};
  const float min_allele_frequency{0};
  const float max_allele_frequency{0};
  const u32 min_alt_counts{0};
};

// TODO: move these two functions to a more appropriate place, e.g. `variant-info`?
std::string GetVariantCorrelationKey(const VariantId& vi, bool pad_left = true);
std::string GetVariantCorrelationKey(
    const std::string& chrom, u64 position, const std::string& ref, const std::string& alt, bool pad_left = true);

// the following functions are only used in the `somatic` workflow
std::vector<std::string> FilterVariant(const VariantId& vid,
                                       const UnifiedVariantFeature& uvf,
                                       const UnifiedReferenceFeature& urf,
                                       const FilterSettings& settings);
std::vector<std::string> FilterPhasedVariant(const std::string& key,
                                             const PhasedFilterSettings& settings,
                                             u32 allele_depth,
                                             float allele_freq,
                                             u8 mapping_qual,
                                             u8 base_qual);
StrUnorderedSet LoadHotspotVariants(const fs::path& vcf_path);
ChromToVariantInfoMap FilterVariantsByVcf(const fs::path& vcf, const ChromToVariantInfoMap& features);
using FilteredVcfFeatures = std::pair<ChromToVcfFeaturesMap, StrUnorderedMap<std::string>>;
std::unordered_map<int, double> CalculateWeightedCountThresholdsPerSubstitutionType(
    const ChromToVariantInfoMap& features, u32 panel_size, double default_threshold);
}  // namespace xoos::svc
