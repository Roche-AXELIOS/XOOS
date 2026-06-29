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
const std::string kFilteringFailSomaticTNMLId = "Somatic_TN_ML";
const std::string kFilteringFailSomaticTNNormalADId = "NormalAD";
const std::string kFilteringFailSomaticTNGermlineId = "GermlineTagging";
const std::string kFilteringMissingFeatureId = "missing_feature";
const std::string kFilteringFalsePositiveId = "false_positive";
const std::string kFilteringMultialleleFormatId = "multiallele_format";
const std::string kFilteringMultiallelePartnerId = "multiallele_partner";
const std::string kFilteringMultialleleConflictId = "multiallele_conflict";
const std::string kFilteringNonAcgtRefAltId = "non_acgt_ref_alt";

// VCF output descriptions
const std::string kFilteringPassDesc = "All filters passed";
const std::string kFilteringFailDesc = "One or more filters failed";
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
const std::string kFilteringFailSomaticTNMLDesc = "Filtered due to low Somatic TN ML Score";
const std::string kFilteringFailSomaticTNNormalADDesc = "Filtered due to high Normal AD counts in somatic variant call";
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

  const u8 min_mapq{0};
  const u8 min_baseq{0};
  const std::unordered_map<s32, f64> weighted_counts_thresholds{};
  const f32 ml_threshold{0};
  const StrUnorderedSet blocklist{};
  const f32 hotspot_weighted_counts_threshold{0};
  const f32 hotspot_ml_threshold{0};
  const f32 min_af_threshold{0};
  const StrUnorderedSet hotspots{};
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

  const u8 min_mapq{0};
  const u8 min_baseq{0};
  const StrUnorderedSet blocklist{};
  const f32 min_allele_frequency{0};
  const f32 max_allele_frequency{0};
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
                                             f32 allele_freq,
                                             u8 mapping_qual,
                                             u8 base_qual);
StrUnorderedSet LoadHotspotVariants(const fs::path& vcf_path);

std::unordered_map<s32, f64> CalculateWeightedCountThresholdsPerSubstitutionType(
    const BamRegionFeatureCollection& features, u32 panel_size, f64 default_threshold);
}  // namespace xoos::svc
