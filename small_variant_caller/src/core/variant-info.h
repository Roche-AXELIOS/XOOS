#pragma once

#include <string>
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
enum class UnifiedFeatureCols {
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
  kSimplex,
  kDuplexAF,
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
  kVcfGenotype,
  kVcfQual,
  kVcfGq,
  kVcfHapcomp,
  kVcfHapdom,
  kVcfRu,
  kVcfRpaRef,
  kVcfRpaAlt,
  kVcfStr,
  kVcfAtInterest,
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
enum class SampleContext {
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
 * Return the FeatureColumn for a given feature name.
 * @param feature_name Feature name string
 * @return a FeatureColumn that matches the specified string
 * @throws error::Error if the feature name is not recognized
 * @see FeatureColumn
 */
FeatureColumn GetFeatureColumn(const std::string& feature_name);

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
static const std::unordered_set<UnifiedFeatureCols> kNormalizableFeatureCols{
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
    // TODO : check whether values of kVcfQual and kVcfGq are proportional to read support
};

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
                                                                    UnifiedFeatureCols::kVcfGenotype,
                                                                    UnifiedFeatureCols::kVcfQual,
                                                                    UnifiedFeatureCols::kVcfGq,
                                                                    UnifiedFeatureCols::kVcfHapcomp,
                                                                    UnifiedFeatureCols::kVcfHapdom,
                                                                    UnifiedFeatureCols::kVcfRu,
                                                                    UnifiedFeatureCols::kVcfRpaRef,
                                                                    UnifiedFeatureCols::kVcfRpaAlt,
                                                                    UnifiedFeatureCols::kVcfStr,
                                                                    UnifiedFeatureCols::kVcfAtInterest,
                                                                    UnifiedFeatureCols::kADT,
                                                                    UnifiedFeatureCols::kADTL,
                                                                    UnifiedFeatureCols::kIndelAf};

// VCF features for a single variant
struct VcfFeature {
 public:
  std::string chrom{};              // chromosome name
  u64 pos{0};                       // variant position
  std::string ref{};                // REF allele
  std::string alt{};                // ALT allele
  f32 nalod{0};                     // Negative log 10 odds of artifact in normal with same allele fraction as tumor
  f32 nlod{0};                      // Normal log 10 likelihood ratio of diploid het or hom alt genotypes
  f32 tlod{0};                      // Log 10 likelihood ratio score of variant existing versus not existing
  u64 mpos{0};                      // median distance from end of read
  u8 mmq_ref{0};                    // median mapping quality of REF allele
  u8 mmq_alt{0};                    // median mapping quality of ALT allele
  u8 mbq_ref{0};                    // median base quality of REF allele
  u8 mbq_alt{0};                    // median base quality of ALT allele
  std::string pre_2bp_context{};    // extract sequence context 2-bp before variant
  std::string post_2bp_context{};   // extract sequence context 2-bp after variant
  std::string post_30bp_context{};  // extract sequence context 30-bp after variant
  u32 uniq_3mers{0};                // number of unique 3-mers in 30-bp post-context
  u32 uniq_4mers{0};                // number of unique 4-mers in 30-bp post-context
  u32 uniq_5mers{0};                // number of unique 5-mers in 30-bp post-context
  u32 uniq_6mers{0};                // number of unique 6-mers in 30-bp post-context
  u32 homopolymer{0};               // homopolymer length in 30-bp post-context
  u32 direpeat{0};                  // occurrence of 2-mer repeat in 30-bp post-context
  u32 trirepeat{0};                 // occurrence of 3-mer repeat in 30-bp post-context
  u32 quadrepeat{0};                // occurrence of 4-mer repeat in 30-bp post-context
  u32 variant_density{0};  // number of variants within a 201-bp window (100 bp up/downstream) at the variant site
  u32 ref_ad{0};           // REF allele depth
  u32 alt_ad{0};           // ALT allele depth
  u32 alt_ad2{0};          // other ALT allele depth
  f64 ref_ad_af{0};        // REF allele AF
  f64 alt_ad_af{0};        // ALT allele AF
  f64 alt_ad2_af{0};       // other ALT allele AF
  u32 tumor_alt_ad{0};     // ALT allele depth for tumor reads
  u32 normal_alt_ad{0};    // ALT allele depth for matched normal reads
  f32 tumor_af{0};         // variant allele frequency for tumor reads
  f32 normal_af{0};        // variant allele frequency for matched normal reads
  f32 tn_af_ratio{0};      // tumor-normal AF ratio
  u32 tumor_dp{0};         // total count for tumor reads
  u32 normal_dp{0};        // total count for matched normal reads
  f32 popaf{0};            // population allele frequency (e.g. gnomAD)
  Genotype genotype{};     // genotype string
  f32 qual{0};             // variant qual
  u32 gq{0};
  u32 hapcomp{0};  // Edit distances of each alt allele's most common supporting haplotype from closest germline
                   // haplotype, excluding differences at the site in question
  f32 hapdom{0};   // For each alt allele, fraction of read support that best fits the most-supported haplotype
                   // containing the allele
  std::string ru{};
  u32 rpa_ref{0};
  u32 rpa_alt{0};
  bool str{false};
  bool at_interest{false};

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
  ReadIds read_ids{};
  f64 weighted_depth{0};                 // Sum of (baseq/138)*(mapq/60) across reads
  f64 nonhomopolymer_weighted_depth{0};  // Sum of (baseq/138)*(mapq/60) across reads not ending in homopolymers
  u32 support{0};                        // Count of reads
  u32 nonhomopolymer_support{0};         // Count of reads not ending in homopolymers
  u8 mapq_min{0};                        // Minimum mapq among reads
  u8 mapq_max{0};                        // Maximum mapq among reads
  u32 mapq_sum{0};                       // Sum of mapqs among reads
  u32 mapq_sum_lowbq{0};                 // Sum of mapq for low-baseq reads
  u32 mapq_sum_simplex{0};               // Sum of mapq for simplex reads
  f64 mapq_mean{0};                      // Mean of mapq
  f64 mapq_mean_lowbq{0};                // Mean of map for low-baseq reads
  f64 mapq_mean_simplex{0};              // Mean of map for simplex reads
  u32 mapq_lt60_count{0};                // Number of reads with mapq < 60
  u32 mapq_lt40_count{0};                // Number of reads with mapq < 40
  u32 mapq_lt30_count{0};                // Number of reads with mapq < 30
  u32 mapq_lt20_count{0};                // Number of reads with mapq < 20
  f64 mapq_lt60_ratio{0};                // Ratio of reads with mapq < 60
  f64 mapq_lt40_ratio{0};                // Ratio of reads with mapq < 40
  f64 mapq_lt30_ratio{0};                // Ratio of reads with mapq < 30
  f64 mapq_lt20_ratio{0};                // Ratio of reads with mapq < 20
  u32 baseq_sum{0};                      // Sum of baseqs among reads
  f64 baseq_mean{0};
  u32 baseq_lt20_count{0};  // Number of reads with baseq < 20
  f64 baseq_lt20_ratio{0};
  f64 mq_af{0};
  f64 bq_af{0};
  u64 distance_min{0};           // Minimum distance to alignment end among reads
  u64 distance_max{0};           // Maximum distance to alignment end among reads
  u64 distance_sum{0};           // Sum of distances to alignment end among reads
  u64 distance_sum_lowbq{0};     // Sum of distances to alignment end among low-baseq reads
  u64 distance_sum_simplex{0};   // Sum of distances to alignment end among simplex reads
  f64 distance_mean{0};          // Mean of distances to alignment end among reads
  f64 distance_mean_lowbq{0};    // Mean of distances to alignment end among low-baseq reads
  f64 distance_mean_simplex{0};  // Mean of distances to alignment end among simplex reads
  f64 duplex_lowbq{0};           // Number of low baseq duplex reads
  u32 simplex{0};                // Number of simplex reads
  f64 duplex_af{0};
  f64 duplex_dp{0};
  u32 familysize_sum{0};  // Sum of family sizes of reads
  f64 familysize_mean{0};
  u32 familysize_lt3_count{0};            // Number of reads with family size < 3
  u32 familysize_lt5_count{0};            // Number of reads with family size < 5
  f64 familysize_lt3_ratio{0};            // Number of reads with family size < 3
  f64 familysize_lt5_ratio{0};            // Number of reads with family size < 5
  u8 nonhomopolymer_mapq_min{0};          // Minimum mapq among reads not ending in homopolymers
  u8 nonhomopolymer_mapq_max{0};          // Maximum mapq among reads not ending in homopolymers
  u32 nonhomopolymer_mapq_sum{0};         // Sum of mapqs among reads not ending in homopolymers
  u32 nonhomopolymer_mapq_lt60_count{0};  // Number of reads not ending in homopolymers with mapq < 60
  u32 nonhomopolymer_mapq_lt40_count{0};  // Number of reads not ending in homopolymers with mapq < 40
  u32 nonhomopolymer_mapq_lt30_count{0};  // Number of reads not ending in homopolymers with mapq < 30
  u32 nonhomopolymer_mapq_lt20_count{0};  // Number of reads not ending in homopolymers with mapq < 20
  f64 nonhomopolymer_mapq_lt60_ratio{0};  // Ratio of reads not ending in homopolymers with mapq < 60
  f64 nonhomopolymer_mapq_lt40_ratio{0};  // Ratio of reads not ending in homopolymers with mapq < 40
  f64 nonhomopolymer_mapq_lt30_ratio{0};  // Ratio of reads not ending in homopolymers with mapq < 30
  f64 nonhomopolymer_mapq_lt20_ratio{0};  // Ratio of reads not ending in homopolymers with mapq < 20
  u8 nonhomopolymer_baseq_min{0};         // Minimum baseq among reads not ending in homopolymers
  u8 nonhomopolymer_baseq_max{0};         // Maximum baseq among reads not ending in homopolymers
  u32 nonhomopolymer_baseq_sum{0};        // Sum of baseqs among reads not ending in homopolymers
  f64 nonhomopolymer_baseq_mean{0};
  f64 nonhomopolymer_mapq_mean{0};
  u32 num_alt{0};

  auto operator<=>(const UnifiedReferenceFeature&) const = default;
};

// zero-initialized UnifiedReferenceFeature, this is used when no reference BAM feature is found
static const UnifiedReferenceFeature kZeroUnifiedReferenceFeature{};

// unified features set for variant (ALT) supporting reads
struct UnifiedVariantFeature {
  ReadIds read_ids{};
  f64 weighted_depth{0};     // Sum of (baseq/138)*(mapq/60) across alt-supporting reads
  u32 support{0};            // Count of reads
  u8 mapq_min{0};            // Minimum mapq
  u8 mapq_max{0};            // Maximum mapq
  u32 mapq_sum{0};           // Sum of mapqs
  u32 mapq_sum_lowbq{0};     // Sum of mapq in low-baseq reads
  u32 mapq_sum_simplex{0};   // Sum of mapq in simplex reads
  f64 mapq_mean{0};          // Mean of mapqs
  f64 mapq_mean_lowbq{0};    // Mean of mapq in low-baseq reads
  f64 mapq_mean_simplex{0};  // Mean of mapq in simplex reads
  u32 mapq_lt60_count{0};    // Number of reads with mapq < 60
  f64 mapq_lt60_ratio{0};
  u32 mapq_lt40_count{0};  // Number of reads with mapq < 40
  f64 mapq_lt40_ratio{0};
  u32 mapq_lt30_count{0};  // Number of reads with mapq < 30
  f64 mapq_lt30_ratio{0};
  u32 mapq_lt20_count{0};  // Number of reads with mapq < 20
  f64 mapq_lt20_ratio{0};
  f64 baseq_min{0};         // Minimum (mean) baseq
  f64 baseq_max{0};         // Maximum (mean) baseq
  f64 baseq_sum{0};         // Sum of (mean) baseqs
  f64 baseq_mean{0};        // Mean of baseqs
  u32 baseq_lt20_count{0};  // Number of reads with baseq < 20
  f64 baseq_lt20_ratio{0};
  u64 distance_min{0};           // Minimum distance to alignment end among reads
  u64 distance_max{0};           // Maximum distance to alignment end among reads
  u64 distance_sum{0};           // Sum of distances to alignment end among reads
  u64 distance_sum_lowbq{0};     // Sum of distances to alignment end among low-baseq reads
  u64 distance_sum_simplex{0};   // Sum of distances to alignment end among simplex reads
  f64 distance_mean{0};          // Mean of distances to alignment end among reads
  f64 distance_mean_lowbq{0};    // Mean of distances to alignment end among low-baseq reads
  f64 distance_mean_simplex{0};  // Mean of distances to alignment end among simplex reads
  f64 mq_af{0};
  f64 bq_af{0};
  f64 duplex_af{0};
  u32 familysize_sum{0};        // sum of family sizes of reads
  u32 familysize_lt3_count{0};  // Number of reads with family size < 3
  u32 familysize_lt5_count{0};  // Number of reads with family size < 5
  f64 familysize_lt3_ratio{0};  // Number of reads with family size < 3
  f64 familysize_lt5_ratio{0};  // Number of reads with family size < 5
  f64 familysize_mean{0};       // Mean of read family sizes
  u32 nonduplex{0};             // Number of nonduplex reads
  u32 duplex{0};                // Number of duplex reads
  f64 duplex_lowbq{0};          // Number of low-baseq duplex reads
  u32 simplex{0};               // Number of simplex reads
  u32 plusonly{0};              // Number of plus-only reads
  u32 minusonly{0};             // Number of minus-only reads
  f64 weighted_score{0};        // Weighed score of duplex and non-duplex counts
  f64 strandbias{0};            // Strand bias feature of variant
  std::string context{};        // Sequence context of variant
  u32 context_index{0};         // Context Index of sequence context
  f64 tn_af_ratio{0};           // Ratio of tumor AF to sum of tumor and normal AF
  // `tn_af_ratio` should be computed using the `duplex_af` values from the tumor sample and normal sample.
  u32 support_reverse{0};                    // Number of supporting reads aligned in the reverse orientation
  f64 alignmentbias{0};                      // Alignment bias for supporting reads
  f64 ml_score{0};                           // ML model score for variant
  std::vector<std::string> filter_status{};  // List of failure reasons for variant when filtered
  f64 adt{0};
  f64 adtl{0};
  f64 indel_af{0};  // TODO : Generalize this feature name? Its calculation is not specific to indels.

  auto operator<=>(const UnifiedVariantFeature&) const = default;

  static s32 ContextIndex(const std::string& context);
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
UnifiedFeatureCols GetCol(const std::string& feature_name);

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
