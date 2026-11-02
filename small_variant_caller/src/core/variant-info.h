#pragma once

#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <xoos/types/int.h>
#include <xoos/types/str-container.h>

#include "genotype.h"
#include "variant-id.h"
#include "xoos/types/float.h"

namespace xoos::svc {

/**
 * Synopsis:
 * This file contains data structures and utility functions for handling variants and their feature column enums.
 */

s32 SubstIndex(const std::string& ref, const std::string& alt);

constexpr s32 kIndelSubstIndex{10};  // `ID` has substitution index of 10

// Enum equivalent of all possible feature column names
enum class UnifiedFeatureCols : u8 {
  kChrom,
  kPos,
  kRef,
  kAlt,
  kContext,
  kContextIndex,
  kSubTypeIndex,
  kWeightedDepth,
  kSupport,
  kMapqMin,
  kMapqMax,
  kMapqSum,
  kMapqSumLowbq,
  kMapqSumSimplex,
  kMapqMean,
  kMapqMeanLowbq,
  kMapqMeanSimplex,
  kMapqLT60Count,
  kMapqLT40Count,
  kMapqLT30Count,
  kMapqLT20Count,
  kMapqLT60Ratio,
  kMapqLT40Ratio,
  kMapqLT30Ratio,
  kMapqLT20Ratio,
  kBaseqMin,
  kBaseqMax,
  kBaseqSum,
  kBaseqMean,
  kBaseqLT20Count,
  kDistanceMin,
  kDistanceMax,
  kDistanceSum,
  kDistanceSumLowbq,
  kDistanceSumSimplex,
  kDistanceMean,
  kDistanceMeanLowbq,
  kDistanceMeanSimplex,
  kFamilysizeSum,
  kFamilysizeLT3Count,
  kFamilysizeLT5Count,
  kFamilysizeMean,
  kNonDuplex,
  kDuplex,
  kDuplexLowbq,
  kDuplexConcordant,
  kDuplexSimplex,
  kDuplexDiscordant,
  kSimplex,
  kDuplexAF,
  kDuplexOther,
  kDuplexOtherAF,
  kDuplexGQ,
  kAfqGT01,
  kAfqGT11,
  kAfqGT12,
  kAfqGT,
  kRefDuplexAF,
  kPlusOnly,
  kMinusOnly,
  kWeightedScore,
  kStrandBias,
  kBamTnAfRatio,
  kSupportReverse,
  kAlignmentBias,
  kMLScore,
  kFilterStatus,
  kRefWeightedDepth,
  kRefNonhomopolymerWeightedDepth,
  kRefSupport,
  kRefNonhomopolymerSupport,
  kRefMapqMin,
  kRefMapqMax,
  kRefMapqSum,
  kRefMapqSumLowbq,
  kRefMapqSumSimplex,
  kRefMapqMean,
  kRefMapqMeanLowbq,
  kRefMapqMeanSimplex,
  kRefMapqLT60Count,
  kRefMapqLT40Count,
  kRefMapqLT30Count,
  kRefMapqLT20Count,
  kRefMapqLT60Ratio,
  kRefMapqLT40Ratio,
  kRefMapqLT30Ratio,
  kRefMapqLT20Ratio,
  kRefBaseqSum,
  kRefBaseqLT20Count,
  kRefDistanceMin,
  kRefDistanceMax,
  kRefDistanceSum,
  kRefDistanceSumLowbq,
  kRefDistanceSumSimplex,
  kRefDistanceMean,
  kRefDistanceMeanLowbq,
  kRefDistanceMeanSimplex,
  kRefDuplexLowbq,
  kRefDuplexConcordant,
  kRefDuplexSimplex,
  kRefDuplexDiscordant,
  kRefSimplex,
  kDuplexDP,
  kRefFamilysizeSum,
  kRefFamilysizeLT3Count,
  kRefFamilysizeLT5Count,
  kRefNonhomopolymerMapqMin,
  kRefNonhomopolymerMapqMax,
  kRefNonhomopolymerMapqSum,
  kRefNonhomopolymerMapqLT60Count,
  kRefNonhomopolymerMapqLT40Count,
  kRefNonhomopolymerMapqLT30Count,
  kRefNonhomopolymerMapqLT20Count,
  kRefNonhomopolymerMapqLT60Ratio,
  kRefNonhomopolymerMapqLT40Ratio,
  kRefNonhomopolymerMapqLT30Ratio,
  kRefNonhomopolymerMapqLT20Ratio,
  kRefNonhomopolymerBaseqMin,
  kRefNonhomopolymerBaseqMax,
  kRefNonhomopolymerBaseqSum,
  kBaseqLT20Ratio,
  kFamilysizeLT5Ratio,
  kFamilysizeLT3Ratio,
  kMQAF,
  kBQAF,
  kRefMQAF,
  kRefBQAF,
  kRefBaseqLT20Ratio,
  kRefFamilysizeMean,
  kRefFamilysizeLT3Ratio,
  kRefFamilysizeLT5Ratio,
  kRefNonhomopolymerMapqMean,
  kRefNonhomopolymerBaseqMean,
  kRefBaseqMean,
  kDel1bpAF,
  kDel2bpAF,
  kDel3bpAF,
  kDelGt3bpAF,
  kIns1bpAF,
  kIns2bpAF,
  kIns3bpAF,
  kInsGt3bpAF,
  kDel1bp2bpAFLogRatio,
  kDel2bp3bpAFLogRatio,
  kDel3bpGt3bpAFLogRatio,
  kIns1bp2bpAFLogRatio,
  kIns2bp3bpAFLogRatio,
  kIns3bpGt3bpAFLogRatio,
  kDel1bpIns1bpAFLogRatio,
  kSnvAF,
  kVcfNalod,
  kVcfNlod,
  kVcfTlod,
  kVcfMpos,
  kVcfMmqRef,
  kVcfMmqAlt,
  kVcfMbqRef,
  kVcfMbqAlt,
  kVariantType,
  kVcfPre2bContext,
  kVcfPost2bContext,
  kVcfPost30bContext,
  kVcfUnique3mers,
  kVcfUnique4mers,
  kVcfUnique5mers,
  kVcfUnique6mers,
  kVcfHomopolymer,
  kVcfDirepeat,
  kVcfTrirepeat,
  kVcfQuadrepeat,
  kVcfVariantDensity,
  kVcfRefAd,
  kVcfAltAd,
  kVcfAltAd2,
  kVcfRefAdAf,
  kVcfAltAdAf,
  kVcfAltAd2Af,
  kVcfTumorAltAd,
  kVcfNormalAltAd,
  kVcfTumorAf,
  kVcfNormalAf,
  kVcfTnAfRatio,
  kVcfTumorDp,
  kVcfNormalDp,
  kVcfPopAf,
  kVcfPonOcc,
  kVcfPonOccRatio,
  kVcfPonMeanDuplexAf,
  kVcfPonMinDuplexAf,
  kVcfPonMaxDuplexAf,
  kVcfGenotype,
  kVcfQual,
  kVcfGq,
  kVcfHapcomp,
  kVcfHapdom,
  kVcfMqRankSum,
  kVcfSor,
  kVcfRu,
  kVcfRpaRef,
  kVcfRpaAlt,
  kVcfStr,
  kVcfAtInterest,
  kVcfIsInRepeat,
  kVcfRepeatLen,
  kVcfRepeatUnitLen,
  kVcfIsUnitSizeMatch,
  kVcfInCatalog,
  kVcfGcContent,
  kNumAlt,
  kADT,
  kADTL,
  kAltLen,
  kIndelAf
};

/**
 * @brief Enum class indicating the feature name prefix.
 * kNone: No prefix, plain feature name indicating feature ignoring sample context
 * kTumor: "tumor_" prefix, indicates feature for the tumor sample only
 * kNormal: "normal_" prefix, indicates feature for the normal sample only
 */
enum class SampleContext : u8 {
  kNone,
  kTumor,
  kNormal
};

/**
 * @brief Structure encapsulating the information about a feature column, including its
 * enum representation and the sample context (tumor, normal, or none).
 * @see UnifiedFeatureCols
 * @see SampleContext
 */
struct FeatureColumn {
  UnifiedFeatureCols enum_val;
  SampleContext sample_context;

  auto operator<=>(const FeatureColumn&) const = default;
};

/**
 * Return the FeatureColumn for a given feature name, with an option to infer sample context from prefixes.
 * If infer_sample_context is true, the function will check for "tumor_" and "normal_" prefixes to determine sample
 * context. If infer_sample_context is false, the function will treat the feature name as is, without inferring sample
 * context.
 * @param feature_name Feature name string, possibly with "tumor_" or "normal_" prefix
 * @param infer_sample_context Whether to infer sample context from prefixes
 * @return a FeatureColumn that matches the specified string and inferred sample context
 * @throws error::Error if the feature name is not recognized
 * @see FeatureColumn
 */
FeatureColumn GetFeatureColumn(const std::string& feature_name, bool infer_sample_context);

/**
 * Return the feature name string that matches a given FeatureColumn
 * @param col a FeatureColumn structure
 * @return the feature name string, possibly with a prefix indicating sample context
 * @see FeatureColumn
 */
std::string GetFeatureName(const FeatureColumn& col);

// This is a list of feature columns that are either:
// 1. derived from read counts
// 2. proportional to read depth
using FeatureColSet = std::unordered_set<UnifiedFeatureCols>;
// clang-format off
static const FeatureColSet kNormalizableFeatureCols{
    UnifiedFeatureCols::kWeightedDepth,
    UnifiedFeatureCols::kSupport,
    UnifiedFeatureCols::kMapqSum,
    UnifiedFeatureCols::kMapqSumLowbq,
    UnifiedFeatureCols::kMapqSumSimplex,
    UnifiedFeatureCols::kMapqLT60Count,
    UnifiedFeatureCols::kMapqLT40Count,
    UnifiedFeatureCols::kMapqLT30Count,
    UnifiedFeatureCols::kMapqLT20Count,
    UnifiedFeatureCols::kBaseqSum,
    UnifiedFeatureCols::kBaseqLT20Count,
    UnifiedFeatureCols::kDistanceSum,
    UnifiedFeatureCols::kDistanceSumLowbq,
    UnifiedFeatureCols::kDistanceSumSimplex,
    UnifiedFeatureCols::kFamilysizeSum,
    UnifiedFeatureCols::kFamilysizeLT3Count,
    UnifiedFeatureCols::kFamilysizeLT5Count,
    UnifiedFeatureCols::kNonDuplex,
    UnifiedFeatureCols::kDuplex,
    UnifiedFeatureCols::kDuplexLowbq,
    UnifiedFeatureCols::kDuplexConcordant,
    UnifiedFeatureCols::kDuplexSimplex,
    UnifiedFeatureCols::kDuplexDiscordant,
    UnifiedFeatureCols::kDuplexOther,
    UnifiedFeatureCols::kDuplexGQ,
    UnifiedFeatureCols::kSimplex,
    UnifiedFeatureCols::kPlusOnly,
    UnifiedFeatureCols::kMinusOnly,
    UnifiedFeatureCols::kWeightedScore,
    UnifiedFeatureCols::kSupportReverse,
    UnifiedFeatureCols::kRefWeightedDepth,
    UnifiedFeatureCols::kRefNonhomopolymerWeightedDepth,
    UnifiedFeatureCols::kRefSupport,
    UnifiedFeatureCols::kRefNonhomopolymerSupport,
    UnifiedFeatureCols::kRefMapqSum,
    UnifiedFeatureCols::kRefMapqSumLowbq,
    UnifiedFeatureCols::kRefMapqSumSimplex,
    UnifiedFeatureCols::kRefMapqLT60Count,
    UnifiedFeatureCols::kRefMapqLT40Count,
    UnifiedFeatureCols::kRefMapqLT30Count,
    UnifiedFeatureCols::kRefMapqLT20Count,
    UnifiedFeatureCols::kRefBaseqSum,
    UnifiedFeatureCols::kRefBaseqLT20Count,
    UnifiedFeatureCols::kRefDistanceSum,
    UnifiedFeatureCols::kRefDistanceSumLowbq,
    UnifiedFeatureCols::kRefDistanceSumSimplex,
    UnifiedFeatureCols::kRefDuplexLowbq,
    UnifiedFeatureCols::kRefDuplexConcordant,
    UnifiedFeatureCols::kRefDuplexSimplex,
    UnifiedFeatureCols::kRefDuplexDiscordant,
    UnifiedFeatureCols::kRefSimplex,
    UnifiedFeatureCols::kDuplexDP,
    UnifiedFeatureCols::kRefFamilysizeSum,
    UnifiedFeatureCols::kRefFamilysizeLT3Count,
    UnifiedFeatureCols::kRefFamilysizeLT5Count,
    UnifiedFeatureCols::kRefNonhomopolymerMapqSum,
    UnifiedFeatureCols::kRefNonhomopolymerMapqLT60Count,
    UnifiedFeatureCols::kRefNonhomopolymerMapqLT40Count,
    UnifiedFeatureCols::kRefNonhomopolymerMapqLT30Count,
    UnifiedFeatureCols::kRefNonhomopolymerMapqLT20Count,
    UnifiedFeatureCols::kRefNonhomopolymerBaseqSum,
    UnifiedFeatureCols::kVcfRefAd,
    UnifiedFeatureCols::kVcfAltAd,
    UnifiedFeatureCols::kVcfAltAd2,
    UnifiedFeatureCols::kVcfTumorAltAd,
    UnifiedFeatureCols::kVcfNormalAltAd,
    UnifiedFeatureCols::kVcfTumorDp,
    UnifiedFeatureCols::kVcfNormalDp
};
// clang-format on

// Set of per-position aggregate BAM feature columns computed by UpdateAggregateFeatures (per-class allele fractions,
// AF log-ratios, and SNV AF). These are only computed when at least one of them is a selected feature column, so
// workflows that do not use them pay no cost on the hot path.
// clang-format off
static const FeatureColSet kAggregateFeatureCols{
    UnifiedFeatureCols::kDel1bpAF,
    UnifiedFeatureCols::kDel2bpAF,
    UnifiedFeatureCols::kDel3bpAF,
    UnifiedFeatureCols::kDelGt3bpAF,
    UnifiedFeatureCols::kIns1bpAF,
    UnifiedFeatureCols::kIns2bpAF,
    UnifiedFeatureCols::kIns3bpAF,
    UnifiedFeatureCols::kInsGt3bpAF,
    UnifiedFeatureCols::kDel1bp2bpAFLogRatio,
    UnifiedFeatureCols::kDel2bp3bpAFLogRatio,
    UnifiedFeatureCols::kDel3bpGt3bpAFLogRatio,
    UnifiedFeatureCols::kIns1bp2bpAFLogRatio,
    UnifiedFeatureCols::kIns2bp3bpAFLogRatio,
    UnifiedFeatureCols::kIns3bpGt3bpAFLogRatio,
    UnifiedFeatureCols::kDel1bpIns1bpAFLogRatio,
    UnifiedFeatureCols::kSnvAF
};
// clang-format on

// Set of VCF feature name enums
static const std::unordered_set<UnifiedFeatureCols> kVcfFeatureCols{UnifiedFeatureCols::kVcfNalod,
                                                                    UnifiedFeatureCols::kVcfNlod,
                                                                    UnifiedFeatureCols::kVcfTlod,
                                                                    UnifiedFeatureCols::kVcfMpos,
                                                                    UnifiedFeatureCols::kVcfMmqRef,
                                                                    UnifiedFeatureCols::kVcfMmqAlt,
                                                                    UnifiedFeatureCols::kVcfMbqRef,
                                                                    UnifiedFeatureCols::kVcfMbqAlt,
                                                                    UnifiedFeatureCols::kVariantType,
                                                                    UnifiedFeatureCols::kVcfPre2bContext,
                                                                    UnifiedFeatureCols::kVcfPost2bContext,
                                                                    UnifiedFeatureCols::kVcfPost30bContext,
                                                                    UnifiedFeatureCols::kVcfUnique3mers,
                                                                    UnifiedFeatureCols::kVcfUnique4mers,
                                                                    UnifiedFeatureCols::kVcfUnique5mers,
                                                                    UnifiedFeatureCols::kVcfUnique6mers,
                                                                    UnifiedFeatureCols::kVcfHomopolymer,
                                                                    UnifiedFeatureCols::kVcfDirepeat,
                                                                    UnifiedFeatureCols::kVcfTrirepeat,
                                                                    UnifiedFeatureCols::kVcfQuadrepeat,
                                                                    UnifiedFeatureCols::kVcfVariantDensity,
                                                                    UnifiedFeatureCols::kVcfRefAd,
                                                                    UnifiedFeatureCols::kVcfAltAd,
                                                                    UnifiedFeatureCols::kVcfAltAd2,
                                                                    UnifiedFeatureCols::kVcfRefAdAf,
                                                                    UnifiedFeatureCols::kVcfAltAdAf,
                                                                    UnifiedFeatureCols::kVcfAltAd2Af,
                                                                    UnifiedFeatureCols::kVcfTumorAltAd,
                                                                    UnifiedFeatureCols::kVcfNormalAltAd,
                                                                    UnifiedFeatureCols::kVcfTumorAf,
                                                                    UnifiedFeatureCols::kVcfNormalAf,
                                                                    UnifiedFeatureCols::kVcfTnAfRatio,
                                                                    UnifiedFeatureCols::kVcfTumorDp,
                                                                    UnifiedFeatureCols::kVcfNormalDp,
                                                                    UnifiedFeatureCols::kVcfPopAf,
                                                                    UnifiedFeatureCols::kVcfPonOcc,
                                                                    UnifiedFeatureCols::kVcfPonOccRatio,
                                                                    UnifiedFeatureCols::kVcfPonMeanDuplexAf,
                                                                    UnifiedFeatureCols::kVcfPonMinDuplexAf,
                                                                    UnifiedFeatureCols::kVcfPonMaxDuplexAf,
                                                                    UnifiedFeatureCols::kVcfGenotype,
                                                                    UnifiedFeatureCols::kVcfQual,
                                                                    UnifiedFeatureCols::kVcfGq,
                                                                    UnifiedFeatureCols::kVcfHapcomp,
                                                                    UnifiedFeatureCols::kVcfHapdom,
                                                                    UnifiedFeatureCols::kVcfMqRankSum,
                                                                    UnifiedFeatureCols::kVcfSor,
                                                                    UnifiedFeatureCols::kVcfRu,
                                                                    UnifiedFeatureCols::kVcfRpaRef,
                                                                    UnifiedFeatureCols::kVcfRpaAlt,
                                                                    UnifiedFeatureCols::kVcfStr,
                                                                    UnifiedFeatureCols::kVcfAtInterest,
                                                                    UnifiedFeatureCols::kVcfIsInRepeat,
                                                                    UnifiedFeatureCols::kVcfRepeatLen,
                                                                    UnifiedFeatureCols::kVcfRepeatUnitLen,
                                                                    UnifiedFeatureCols::kVcfIsUnitSizeMatch,
                                                                    UnifiedFeatureCols::kVcfInCatalog,
                                                                    UnifiedFeatureCols::kVcfGcContent,
                                                                    UnifiedFeatureCols::kADT,
                                                                    UnifiedFeatureCols::kADTL,
                                                                    UnifiedFeatureCols::kIndelAf};

// VCF features for a single variant
struct VcfFeature {
 public:
  std::string chrom;              // chromosome name
  u64 pos{0};                     // variant position
  std::string ref;                // REF allele
  std::string alt;                // ALT allele
  f32 nalod{0};                   // Negative log 10 odds of artifact in normal with same allele fraction as tumor
  f32 nlod{0};                    // Normal log 10 likelihood ratio of diploid het or hom alt genotypes
  f32 tlod{0};                    // Log 10 likelihood ratio score of variant existing versus not existing
  u64 mpos{0};                    // median distance from end of read
  u8 mmq_ref{0};                  // median mapping quality of REF allele
  u8 mmq_alt{0};                  // median mapping quality of ALT allele
  u8 mbq_ref{0};                  // median base quality of REF allele
  u8 mbq_alt{0};                  // median base quality of ALT allele
  std::string pre_2bp_context;    // extract sequence context 2-bp before variant
  std::string post_2bp_context;   // extract sequence context 2-bp after variant
  std::string post_30bp_context;  // extract sequence context 30-bp after variant
  u32 uniq_3mers{0};              // number of unique 3-mers in 30-bp post-context
  u32 uniq_4mers{0};              // number of unique 4-mers in 30-bp post-context
  u32 uniq_5mers{0};              // number of unique 5-mers in 30-bp post-context
  u32 uniq_6mers{0};              // number of unique 6-mers in 30-bp post-context
  u32 homopolymer{0};             // homopolymer length in 30-bp post-context
  u32 direpeat{0};                // occurrence of 2-mer repeat in 30-bp post-context
  u32 trirepeat{0};               // occurrence of 3-mer repeat in 30-bp post-context
  u32 quadrepeat{0};              // occurrence of 4-mer repeat in 30-bp post-context
  u32 variant_density{0};  // number of variants within a 201-bp window (100 bp up/downstream) at the variant site
  u32 ref_ad{0};           // REF allele depth
  u32 alt_ad{0};           // ALT allele depth
  u32 alt_ad2{0};          // other ALT allele depth
  f32 ref_ad_af{0};        // REF allele AF
  f32 alt_ad_af{0};        // ALT allele AF
  f32 alt_ad2_af{0};       // other ALT allele AF
  u32 tumor_alt_ad{0};     // ALT allele depth for tumor reads
  u32 normal_alt_ad{0};    // ALT allele depth for matched normal reads
  f32 tumor_af{0};         // variant allele frequency for tumor reads
  f32 normal_af{0};        // variant allele frequency for matched normal reads
  f32 tn_af_ratio{0};      // tumor-normal AF ratio
  u32 tumor_dp{0};         // total count for tumor reads
  u32 normal_dp{0};        // total count for matched normal reads
  f32 popaf{0};            // population allele frequency (e.g. gnomAD)
  // PON occurrence count (number of panel samples with this variant)
  u32 pon_occ{0};
  // PON occurrence ratio (fraction of panel samples with this variant)
  f32 pon_occ_ratio{0};
  // PON mean duplex allele frequency
  f32 pon_mean_duplex_af{0};
  // PON min duplex allele frequency
  f32 pon_min_duplex_af{0};
  // PON max duplex allele frequency
  f32 pon_max_duplex_af{0};
  Genotype genotype{};  // genotype string
  f32 qual{0};          // variant qual
  u32 gq{0};
  u32 hapcomp{0};      // Edit distances of each alt allele's most common supporting haplotype from closest germline
                       // haplotype, excluding differences at the site in question
  f32 hapdom{0};       // For each alt allele, fraction of read support that best fits the most-supported haplotype
                       // containing the allele
  f32 mq_rank_sum{0};  // GATK mapping-quality rank-sum statistic
  f32 sor{0};          // Symmetric Odds Ratio of 2x2 contingency table to detect strand bias
  std::string ru;
  s32 rpa_ref{0};
  s32 rpa_alt{0};
  bool str{false};
  bool at_interest{false};
  bool is_in_repeat{false};        // whether the variant is in a short tandem repeat (repeat_len >= threshold)
  u32 repeat_len{0};               // bp span of the dominant tandem repeat at the variant site
  u32 repeat_unit_len{0};          // size of the repeating motif (1-4 bp); 0 if not in a repeat
  bool is_unit_size_match{false};  // whether the indel size matches the repeat unit size
  // whether the variant is present in the supplied catalog VCF (--catalog-vcf)
  bool in_catalog{false};
  f32 gc_content{0};  // GC content in a 200-bp window centered on the variant's start position

  auto operator<=>(const VcfFeature&) const = default;
};

// Zero-initialized VcfFeature, this is used when no VcfFeature feature is found
static const VcfFeature kZeroVcfFeature{};

// Duplicate alignments can be produced by both bwa and GATK HaplotypeCaller/Mutect2,
// we do not want to count them twice in the variant support count. To avoid this
// we will use the read name to identify duplicates. Using the read name directly has
// a large impact on performance due to memory allocation and string hashing/comparison.
// To avoid this we assign each read a unique id and use this id to identify duplicates,
// this id is only unique within a region and not across regions, but this is sufficient.
using ReadId = u32;
using ReadIds = std::vector<ReadId>;

// unified features set for reference allele (REF) supporting reads
struct UnifiedReferenceFeature {
  ReadIds read_ids;
  f64 weighted_depth{0};                 // Sum of (baseq/138)*(mapq/60) across reads
  f64 nonhomopolymer_weighted_depth{0};  // Sum of (baseq/138)*(mapq/60) across reads not ending in homopolymers
  u32 support{0};                        // Count of reads
  u32 nonhomopolymer_support{0};         // Count of reads not ending in homopolymers
  u8 mapq_min{0};                        // Minimum mapq among reads
  u8 mapq_max{0};                        // Maximum mapq among reads
  u32 mapq_sum{0};                       // Sum of mapqs among reads
  u32 mapq_sum_lowbq{0};                 // Sum of mapq for low-baseq reads
  u32 mapq_sum_simplex{0};               // Sum of mapq for simplex reads
  f32 mapq_mean{0};                      // Mean of mapq
  f32 mapq_mean_lowbq{0};                // Mean of map for low-baseq reads
  f32 mapq_mean_simplex{0};              // Mean of map for simplex reads
  u32 mapq_lt60_count{0};                // Number of reads with mapq < 60
  u32 mapq_lt40_count{0};                // Number of reads with mapq < 40
  u32 mapq_lt30_count{0};                // Number of reads with mapq < 30
  u32 mapq_lt20_count{0};                // Number of reads with mapq < 20
  f32 mapq_lt60_ratio{0};                // Ratio of reads with mapq < 60
  f32 mapq_lt40_ratio{0};                // Ratio of reads with mapq < 40
  f32 mapq_lt30_ratio{0};                // Ratio of reads with mapq < 30
  f32 mapq_lt20_ratio{0};                // Ratio of reads with mapq < 20
  u32 baseq_sum{0};                      // Sum of baseqs among reads
  f32 baseq_mean{0};
  u32 baseq_lt20_count{0};  // Number of reads with baseq < 20
  f32 baseq_lt20_ratio{0};
  f32 mq_af{0};
  f32 bq_af{0};
  u64 distance_min{0};           // Minimum distance to alignment end among reads
  u64 distance_max{0};           // Maximum distance to alignment end among reads
  u64 distance_sum{0};           // Sum of distances to alignment end among reads
  u64 distance_sum_lowbq{0};     // Sum of distances to alignment end among low-baseq reads
  u64 distance_sum_simplex{0};   // Sum of distances to alignment end among simplex reads
  f64 distance_mean{0};          // Mean of distances to alignment end among reads
  f64 distance_mean_lowbq{0};    // Mean of distances to alignment end among low-baseq reads
  f64 distance_mean_simplex{0};  // Mean of distances to alignment end among simplex reads
  f64 duplex_lowbq{0};           // Number of low baseq duplex reads
  u32 duplex_concordant{0};      // Number of reads with concordant base type at reference position
  u32 duplex_simplex{0};         // Number of reads with simplex base type at reference position
  u32 duplex_discordant{0};      // Number of reads with discordant base type at reference position
  u32 simplex{0};                // Number of simplex reads
  f32 duplex_af{0};
  f64 duplex_dp{0};
  u32 familysize_sum{0};  // Sum of family sizes of reads
  f32 familysize_mean{0};
  u32 familysize_lt3_count{0};            // Number of reads with family size < 3
  u32 familysize_lt5_count{0};            // Number of reads with family size < 5
  f32 familysize_lt3_ratio{0};            // Fraction of reads with family size < 3
  f32 familysize_lt5_ratio{0};            // Fraction of reads with family size < 5
  u8 nonhomopolymer_mapq_min{0};          // Minimum mapq among reads not ending in homopolymers
  u8 nonhomopolymer_mapq_max{0};          // Maximum mapq among reads not ending in homopolymers
  u32 nonhomopolymer_mapq_sum{0};         // Sum of mapqs among reads not ending in homopolymers
  u32 nonhomopolymer_mapq_lt60_count{0};  // Number of reads not ending in homopolymers with mapq < 60
  u32 nonhomopolymer_mapq_lt40_count{0};  // Number of reads not ending in homopolymers with mapq < 40
  u32 nonhomopolymer_mapq_lt30_count{0};  // Number of reads not ending in homopolymers with mapq < 30
  u32 nonhomopolymer_mapq_lt20_count{0};  // Number of reads not ending in homopolymers with mapq < 20
  f32 nonhomopolymer_mapq_lt60_ratio{0};  // Ratio of reads not ending in homopolymers with mapq < 60
  f32 nonhomopolymer_mapq_lt40_ratio{0};  // Ratio of reads not ending in homopolymers with mapq < 40
  f32 nonhomopolymer_mapq_lt30_ratio{0};  // Ratio of reads not ending in homopolymers with mapq < 30
  f32 nonhomopolymer_mapq_lt20_ratio{0};  // Ratio of reads not ending in homopolymers with mapq < 20
  u8 nonhomopolymer_baseq_min{0};         // Minimum baseq among reads not ending in homopolymers
  u8 nonhomopolymer_baseq_max{0};         // Maximum baseq among reads not ending in homopolymers
  u32 nonhomopolymer_baseq_sum{0};        // Sum of baseqs among reads not ending in homopolymers
  f32 nonhomopolymer_baseq_mean{0};
  f32 nonhomopolymer_mapq_mean{0};
  u32 num_alt{0};
  // Per-position duplex depth aggregated across all alt alleles by indel class. Unlike duplex_dp, these aggregate
  // features use the duplex count only and never the duplex_lowbq reads, independent of the --include-duplex-lowbq
  // CLI option. The matching allele frequencies below are normalized by the total duplex depth at the position
  // (reference support + summed non-deletion alt support; deletion alt support is excluded to avoid double-counting
  // deletion reads already present in the reference term at the anchor base).
  // Support fields hold integer duplex read counts, so u32 is sufficient and halves their footprint vs f64.
  u32 del_1bp_support{0};    // Total duplex depth of 1bp deletion alt alleles at this position
  u32 del_2bp_support{0};    // Total duplex depth of 2bp deletion alt alleles at this position
  u32 del_3bp_support{0};    // Total duplex depth of 3bp deletion alt alleles at this position
  u32 del_gt3bp_support{0};  // Total duplex depth of >3bp deletion alt alleles at this position
  u32 ins_1bp_support{0};    // Total duplex depth of 1bp insertion alt alleles at this position
  u32 ins_2bp_support{0};    // Total duplex depth of 2bp insertion alt alleles at this position
  u32 ins_3bp_support{0};    // Total duplex depth of 3bp insertion alt alleles at this position
  u32 ins_gt3bp_support{0};  // Total duplex depth of >3bp insertion alt alleles at this position
  // Per-position allele frequencies of each indel class (class support / total duplex depth).
  // Stored as f32: these are bounded fractions/log-ratios and f32 precision is sufficient for the ML features.
  f32 del_1bp_af{0};    // AF of 1bp deletion alt alleles at this position
  f32 del_2bp_af{0};    // AF of 2bp deletion alt alleles at this position
  f32 del_3bp_af{0};    // AF of 3bp deletion alt alleles at this position
  f32 del_gt3bp_af{0};  // AF of >3bp deletion alt alleles at this position
  f32 ins_1bp_af{0};    // AF of 1bp insertion alt alleles at this position
  f32 ins_2bp_af{0};    // AF of 2bp insertion alt alleles at this position
  f32 ins_3bp_af{0};    // AF of 3bp insertion alt alleles at this position
  f32 ins_gt3bp_af{0};  // AF of >3bp insertion alt alleles at this position
  // Log ratios of adjacent indel-class AFs (log((a + eps) / (b + eps)))
  f32 del_1bp_2bp_af_logratio{0};      // log ratio of del_1bp_af to del_2bp_af
  f32 del_2bp_3bp_af_logratio{0};      // log ratio of del_2bp_af to del_3bp_af
  f32 del_3bp_gt3bp_af_logratio{0};    // log ratio of del_3bp_af to del_gt3bp_af
  f32 ins_1bp_2bp_af_logratio{0};      // log ratio of ins_1bp_af to ins_2bp_af
  f32 ins_2bp_3bp_af_logratio{0};      // log ratio of ins_2bp_af to ins_3bp_af
  f32 ins_3bp_gt3bp_af_logratio{0};    // log ratio of ins_3bp_af to ins_gt3bp_af
  f32 del_1bp_ins_1bp_af_logratio{0};  // log ratio of del_1bp_af to ins_1bp_af
  u32 snv_support{0};                  // Total duplex depth of SNV alt alleles at this position
  f32 snv_af{0};                       // AF of SNV alt alleles at this position (snv_support / total duplex depth)

  auto operator<=>(const UnifiedReferenceFeature&) const = default;
};

// zero-initialized UnifiedReferenceFeature, this is used when no reference BAM feature is found
static const UnifiedReferenceFeature kZeroUnifiedReferenceFeature{};

// duplex_lowbq is stored as a fractional count (0.5 per read); multiply by this factor to
// convert to an integer read count for VCF output.
static constexpr u32 kDuplexLowbqToCountFactor = 2;

// unified features set for variant (ALT) supporting reads
struct UnifiedVariantFeature {
  ReadIds read_ids;
  f64 weighted_depth{0};     // Sum of (baseq/138)*(mapq/60) across alt-supporting reads
  u32 support{0};            // Count of reads
  u8 mapq_min{0};            // Minimum mapq
  u8 mapq_max{0};            // Maximum mapq
  u32 mapq_sum{0};           // Sum of mapqs
  u32 mapq_sum_lowbq{0};     // Sum of mapq in low-baseq reads
  u32 mapq_sum_simplex{0};   // Sum of mapq in simplex reads
  f32 mapq_mean{0};          // Mean of mapqs
  f32 mapq_mean_lowbq{0};    // Mean of mapq in low-baseq reads
  f32 mapq_mean_simplex{0};  // Mean of mapq in simplex reads
  u32 mapq_lt60_count{0};    // Number of reads with mapq < 60
  f32 mapq_lt60_ratio{0};
  u32 mapq_lt40_count{0};  // Number of reads with mapq < 40
  f32 mapq_lt40_ratio{0};
  u32 mapq_lt30_count{0};  // Number of reads with mapq < 30
  f32 mapq_lt30_ratio{0};
  u32 mapq_lt20_count{0};  // Number of reads with mapq < 20
  f32 mapq_lt20_ratio{0};
  f64 baseq_min{0};         // Minimum (mean) baseq
  f64 baseq_max{0};         // Maximum (mean) baseq
  f64 baseq_sum{0};         // Sum of (mean) baseqs
  f32 baseq_mean{0};        // Mean of baseqs
  u32 baseq_lt20_count{0};  // Number of reads with baseq < 20
  f32 baseq_lt20_ratio{0};
  u64 distance_min{0};           // Minimum distance to alignment end among reads
  u64 distance_max{0};           // Maximum distance to alignment end among reads
  u64 distance_sum{0};           // Sum of distances to alignment end among reads
  u64 distance_sum_lowbq{0};     // Sum of distances to alignment end among low-baseq reads
  u64 distance_sum_simplex{0};   // Sum of distances to alignment end among simplex reads
  f64 distance_mean{0};          // Mean of distances to alignment end among reads
  f64 distance_mean_lowbq{0};    // Mean of distances to alignment end among low-baseq reads
  f64 distance_mean_simplex{0};  // Mean of distances to alignment end among simplex reads
  f32 mq_af{0};
  f32 bq_af{0};
  f32 duplex_af{0};
  f64 duplex_other{0};               // Duplex read count from other alleles at the same position
  f32 duplex_other_af{0};            // Allele frequency of other alleles' duplex reads
  f64 duplex_gq{0};                  // Genotype quality based on duplex AF distance to ideal GTs
  f64 afq_gt01{0};                   // 2.0 - manhattan distance to ideal GT=0/1 (0.5, 0.0, 0.5)
  f64 afq_gt11{0};                   // 2.0 - manhattan distance to ideal GT=1/1 (1.0, 0.0, 0.0)
  f64 afq_gt12{0};                   // 2.0 - manhattan distance to ideal GT=1/2 (0.5, 0.5, 0.0)
  Genotype afq_gt{Genotype::kGT00};  // Categorical genotype call from AF distances
  u32 familysize_sum{0};             // sum of family sizes of reads
  u32 familysize_lt3_count{0};       // Number of reads with family size < 3
  u32 familysize_lt5_count{0};       // Number of reads with family size < 5
  f32 familysize_lt3_ratio{0};       // Fraction of reads with family size < 3
  f32 familysize_lt5_ratio{0};       // Fraction of reads with family size < 5
  f32 familysize_mean{0};            // Mean of read family sizes
  u32 nonduplex{0};                  // Number of nonduplex reads
  u32 duplex{0};                     // Number of duplex reads
  f64 duplex_lowbq{0};               // Number of low-baseq duplex reads
  u32 duplex_concordant{0};          // Number of reads with concordant base type at variant position
  u32 duplex_simplex{0};             // Number of reads with simplex base type at variant position
  u32 duplex_discordant{0};          // Number of reads with discordant base type at variant position
  u32 simplex{0};                    // Number of simplex reads
  u32 plusonly{0};                   // Number of plus-only reads
  u32 minusonly{0};                  // Number of minus-only reads
  f64 weighted_score{0};             // Weighed score of duplex and non-duplex counts
  f64 strandbias{0};                 // Strand bias feature of variant
  std::string context;               // Sequence context of variant
  u32 context_index{0};              // Context Index of sequence context
  f32 tn_af_ratio{0};                // Ratio of tumor AF to sum of tumor and normal AF
  // `tn_af_ratio` should be computed using the `duplex_af` values from the tumor sample and normal sample.
  u32 support_reverse{0};                  // Number of supporting reads aligned in the reverse orientation
  f64 alignmentbias{0};                    // Alignment bias for supporting reads
  f64 ml_score{0};                         // ML model score for variant
  std::vector<std::string> filter_status;  // List of failure reasons for variant when filtered
  f64 adt{0};
  f64 adtl{0};
  f32 indel_af{0};

  auto operator<=>(const UnifiedVariantFeature&) const = default;

  static u32 ContextIndex(const std::string& context);
};

// Zero-initialized UnifiedVariantFeature, this is used when no variant BAM feature is found
static const UnifiedVariantFeature kZeroUnifiedVariantFeature{};

/**
 * @brief Mapping from UnifiedFeatureCols enum to feature name string
 */
extern const std::unordered_map<UnifiedFeatureCols, std::string> kUnifiedFeatureColsToString;

/**
 * @brief Mapping from feature name string to UnifiedFeatureCols enum
 */
extern const StrUnorderedMap<UnifiedFeatureCols> kStringToUnifiedFeatureCols;

/**
 * @brief Check if a UnifiedFeatureCols enum corresponds to a VCF feature column
 * @param col UnifiedFeatureCols enum value
 * @return true if the column is a VCF feature column, false otherwise
 */
bool IsVcfFeatureCol(UnifiedFeatureCols col);

/**
 * @brief Check if a FeatureColumn corresponds to a VCF feature column
 * @param col FeatureColumn structure
 * @return true if the column is a VCF feature column, false otherwise
 */
bool IsVcfFeatureColumn(const FeatureColumn& col);

/**
 * Return the UnifiedFeatureCols enum for a given feature name.
 * @param feature_name Feature name string
 * @return a UnifiedFeatureCols that matches the specified string
 */
UnifiedFeatureCols GetCol(std::string_view feature_name);

/**
 * Return the feature/column string that matches a given UnifiedFeatureCols
 * @param col a UnifiedFeatureCols enum value
 * @return a feature/column name string
 */
std::string GetFeatureName(UnifiedFeatureCols col);

using VarIdToVcfFeatures = std::unordered_map<VariantId, VcfFeature>;
using VarIdToVarBamFeatures = std::unordered_map<VariantId, UnifiedVariantFeature>;
using PosToRefBamFeatures = std::unordered_map<u64, UnifiedReferenceFeature>;
using ChromPosToRefBamFeatures = StrUnorderedMap<PosToRefBamFeatures>;

}  // namespace xoos::svc
