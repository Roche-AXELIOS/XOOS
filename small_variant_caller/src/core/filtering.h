#pragma once

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>

#include "bam-feature-collection.h"
#include "core/variant-info.h"

namespace xoos::svc {

// VCF outout FILTER field IDs
constexpr auto* kFilteringMapQualityId = "mapquality";
constexpr auto* kFilteringBaseQualityId = "basequality";
constexpr auto* kFilteringAFId = "Allele_Freq";
constexpr auto* kFilteringMinAltCountsId = "Min_Alt_Counts";
constexpr auto* kFilteringCountsId = "Counts";
constexpr auto* kFilteringMLScoreId = "ML_Score";
constexpr auto* kFilteringBlocklistedId = "BLOCKLISTED";
constexpr auto* kFilteringForcedId = "FORCED";
constexpr auto* kFilteringPassId = "PASS";
constexpr auto* kFilteringFailId = "FAIL";
constexpr auto* kFilteringFailMinTumorSupportId = "min_tumor_support";
constexpr auto* kFilteringFailMinSupportId = "min_support";
constexpr auto* kFilteringFailMaxNormalSupportId = "max_normal_support";
constexpr auto* kFilteringFailMinAfId = "min_af";
constexpr auto* kFilteringFailMinTumorAfId = "min_tumor_af";
constexpr auto* kFilteringFailMinDpRatioId = "min_dp_ratio";
constexpr auto* kFilteringFailSomaticTNGermlineId = "GermlineTagging";
constexpr auto* kFilteringFailMaxIndelSizeId = "max_indel_size";
constexpr auto* kFilteringMissingFeatureId = "missing_feature";
constexpr auto* kFilteringFalsePositiveId = "false_positive";
constexpr auto* kFilteringMultialleleFormatId = "multiallele_format";
constexpr auto* kFilteringMultiallelePartnerId = "multiallele_partner";
constexpr auto* kFilteringMultialleleConflictId = "multiallele_conflict";
constexpr auto* kFilteringNonAcgtRefAltId = "non_acgt_ref_alt";

// FILTER tags applied by the germline-mitochondria AF/depth filtering.
constexpr auto* kFilterMitoLowAf = "low_af";
constexpr auto* kFilterMitoLowDp = "low_dp";
constexpr auto* kFilterMitoBlocklist = "blocklist";

// VCF output descriptions
constexpr auto* kFilteringPassDesc = "All filters passed";
constexpr auto* kFilteringFailDesc = "One or more filters failed";
constexpr auto* kFilteringGermlinePassDesc = "Site contains at least one allele that passes filters";
constexpr auto* kFilteringGermlineFailDesc = "Variant is filtered out by ML";
constexpr auto* kFilteringMapQualityDesc = "Filtered due to mapping quality";
constexpr auto* kFilteringAFDesc = "Filtered due to allele frequency";
constexpr auto* kFilteringMinAltCountsDesc = "Filtered due to alt counts";
constexpr auto* kFilteringBaseQualityDesc = "Filtered due to base quality";
constexpr auto* kFilteringCountsDesc = "Filtered due to low counts";
constexpr auto* kFilteringMLScoreDesc = "Filtered due to low ML Score";
constexpr auto* kFilteringForcedDesc = "Filtered due to being a forced variant call";
constexpr auto* kFilteringBlocklistedDesc = "Variant blocklisted";
constexpr auto* kFilteringFailMinTumorSupportDesc = "Filtered due to low supporting reads in tumor sample";
constexpr auto* kFilteringFailMinSupportDesc = "Filtered due to low supporting reads";
constexpr auto* kFilteringFailMaxNormalSupportDesc = "Filtered due to high supporting reads in normal sample";
constexpr auto* kFilteringFailMinAfDesc = "Filtered due to low allele frequency";
constexpr auto* kFilteringFailMinTumorAfDesc = "Filtered due to low allele frequency in tumor sample";
constexpr auto* kFilteringMinDpRatioDesc = "Filtered due to low DP relative to chromosome median DP";
constexpr auto* kFilteringFailSomaticTNGermlineDesc = "Filtered due to being classified as a possible germline variant";
constexpr auto* kFilteringFailMaxIndelSizeDesc = "Filtered due to large indel size";
constexpr auto* kFilteringMissingFeatureDesc = "Filtered due to missing ML feature";
constexpr auto* kFilteringFalsePositiveDesc = "Filtered due to ML model classification as false positive";
constexpr auto* kFilteringMultialleleFormatDesc = "Filtered due to multi-allele record format failure";
constexpr auto* kFilteringMultiallelePartnerDesc = "Filtered due to missing multi-allele partner";
constexpr auto* kFilteringMultialleleConflictDesc = "Filtered due to conflicting multi-allele ML model classification";
constexpr auto* kFilteringNonAcgtRefAltDesc = "Filtered due to non-ACGT reference or alternate allele(s)";
constexpr auto* kFilterMitoBlocklistDesc = "Variant start is in a blocklist BED region";

// Settings for `somatic` workflow variants filtering
struct FilterSettings {
  FilterSettings(const u8 min_mapq,
                 const u8 min_baseq,
                 const std::unordered_map<s32, f64>& weighted_counts_thresholds,
                 const f32 ml_threshold,
                 StrUnorderedSet blocklist,
                 const f32 hotspot_weighted_counts_threshold,
                 const f32 hotspot_ml_threshold,
                 const f32 min_af_threshold,
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

  u8 min_mapq{0};
  u8 min_baseq{0};
  std::unordered_map<s32, f64> weighted_counts_thresholds;
  f32 ml_threshold{0};
  StrUnorderedSet blocklist;
  f32 hotspot_weighted_counts_threshold{0};
  f32 hotspot_ml_threshold{0};
  f32 min_af_threshold{0};
  StrUnorderedSet hotspots;
};

// Settings for `somatic` workflow phased variants filtering
struct PhasedFilterSettings {
  PhasedFilterSettings(const u8 min_mapq,
                       const u8 min_baseq,
                       StrUnorderedSet blocklist,
                       const f32 min_af,
                       const f32 max_af,
                       const u32 min_alt_counts)
      : min_mapq(min_mapq),
        min_baseq(min_baseq),
        blocklist(std::move(blocklist)),
        min_allele_frequency(min_af),
        max_allele_frequency(max_af),
        min_alt_counts(min_alt_counts) {
  }

  u8 min_mapq{0};
  u8 min_baseq{0};
  StrUnorderedSet blocklist;
  f32 min_allele_frequency{0};
  f32 max_allele_frequency{0};
  u32 min_alt_counts{0};
};

std::string GetVariantCorrelationKey(const VariantId& vid, bool pad_left);
std::string GetVariantCorrelationKey(
    const std::string& chrom, u64 position, const std::string& ref, const std::string& alt, bool pad_left);
std::string GetVariantCorrelationKey(const std::string& chrom,
                                     u64 position,
                                     const std::string& ref,
                                     const std::string& alt);

// the following functions are only used in the `somatic` workflow
std::vector<std::string> FilterVariant(const VariantId& vid,
                                       const UnifiedVariantFeature& uvf,
                                       const UnifiedReferenceFeature& urf,
                                       const FilterSettings& settings);
std::vector<std::string> FilterPhasedVariant(const std::string& key,
                                             const PhasedFilterSettings& settings,
                                             u32 allele_depth,
                                             f32 allele_freq,
                                             u8 mapping_qual,
                                             u8 base_qual);
StrUnorderedSet LoadHotspotVariants(const fs::path& vcf_path);

std::unordered_map<s32, f64> CalculateWeightedCountThresholdsPerSubstitutionType(
    const BamRegionFeatureCollection& features, u32 panel_size, f64 default_threshold);
}  // namespace xoos::svc
