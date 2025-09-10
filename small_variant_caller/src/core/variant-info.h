#pragma once

#include <map>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <xoos/types/int.h>
#include <xoos/types/str-container.h>

#include "genotype.h"
#include "variant-id.h"

namespace xoos::svc {

/**
 * Synopsis:
 * This file contains data structures and utility functions for handling variants and their feature column enums.
 */

int SubstIndex(const std::string& ref, const std::string& alt);

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
  kTumorSupport,
  kTumorDistanceSum,
  kTumorDistanceMean,
  kTumorBaseqSum,
  kTumorBaseqMean,
  kTumorMapqSum,
  kTumorMapqMean,
  kTumorRefSupport,
  kNormalSupport,
  kNormalDistanceSum,
  kNormalDistanceMean,
  kNormalBaseqSum,
  kNormalBaseqMean,
  kNormalMapqSum,
  kNormalMapqMean,
  kNormalRefSupport,
  kBamTumorAF,
  kBamNormalAF,
  kBamRAT,
  kSupportReverse,
  kTumorSupportReverse,
  kNormalSupportReverse,
  kAlignmentBias,
  kTumorAlignmentBias,
  kNormalAlignmentBias,
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
  kVcfTumorNormalAfRatio,
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
    UnifiedFeatureCols::kTumorSupport,
    UnifiedFeatureCols::kTumorDistanceSum,
    UnifiedFeatureCols::kTumorBaseqSum,
    UnifiedFeatureCols::kTumorMapqSum,
    UnifiedFeatureCols::kTumorRefSupport,
    UnifiedFeatureCols::kNormalSupport,
    UnifiedFeatureCols::kNormalDistanceSum,
    UnifiedFeatureCols::kNormalBaseqSum,
    UnifiedFeatureCols::kNormalMapqSum,
    UnifiedFeatureCols::kNormalRefSupport,
    UnifiedFeatureCols::kSupportReverse,
    UnifiedFeatureCols::kTumorSupportReverse,
    UnifiedFeatureCols::kNormalSupportReverse,
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
                                                                    UnifiedFeatureCols::kVcfTumorNormalAfRatio,
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
                                                                    UnifiedFeatureCols::kNumAlt,
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
  float nalod{0};                   // Negative log 10 odds of artifact in normal with same allele fraction as tumor
  float nlod{0};                    // Normal log 10 likelihood ratio of diploid het or hom alt genotypes
  float tlod{0};                    // Log 10 likelihood ratio score of variant existing versus not existing
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
  double ref_ad_af{0};     // REF allele AF
  double alt_ad_af{0};     // ALT allele AF
  double alt_ad2_af{0};    // other ALT allele AF
  u32 tumor_alt_ad{0};     // ALT allele depth for tumor reads
  u32 normal_alt_ad{0};    // ALT allele depth for matched normal reads
  float tumor_af{0};       // variant allele frequency for tumor reads
  float normal_af{0};      // variant allele frequency for matched normal reads
  float tumor_normal_af_ratio{0};  // tumor-normal AF ratio
  u32 tumor_dp{0};                 // total count for tumor reads
  u32 normal_dp{0};                // total count for matched normal reads
  float popaf{0};                  // population allele frequency (e.g. gnomAD)
  Genotype genotype{};             // genotype string
  float qual{0};                   // variant qual
  u32 gq{0};
  u32 hapcomp{0};   // Edit distances of each alt allele's most common supporting haplotype from closest germline
                    // haplotype, excluding differences at the site in question
  float hapdom{0};  // For each alt allele, fraction of read support that best fits the most-supported haplotype
                    // containing the allele
  std::string ru{};
  u32 rpa_ref{0};
  u32 rpa_alt{0};
  bool str{false};
  bool at_interest{false};
};

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
  double weighted_depth{0};                 // Sum of (baseq/138)*(mapq/60) across reads
  double nonhomopolymer_weighted_depth{0};  // Sum of (baseq/138)*(mapq/60) across reads not ending in homopolymers
  u32 support{0};                           // Count of reads
  u32 nonhomopolymer_support{0};            // Count of reads not ending in homopolymers
  u8 mapq_min{0};                           // Minimum mapq among reads
  u8 mapq_max{0};                           // Maximum mapq among reads
  u32 mapq_sum{0};                          // Sum of mapqs among reads
  u32 mapq_sum_lowbq{0};                    // Sum of mapq for low-baseq reads
  u32 mapq_sum_simplex{0};                  // Sum of mapq for simplex reads
  double mapq_mean{0};                      // Mean of mapq
  double mapq_mean_lowbq{0};                // Mean of map for low-baseq reads
  double mapq_mean_simplex{0};              // Mean of map for simplex reads
  u32 mapq_lt60_count{0};                   // Number of reads with mapq < 60
  u32 mapq_lt40_count{0};                   // Number of reads with mapq < 40
  u32 mapq_lt30_count{0};                   // Number of reads with mapq < 30
  u32 mapq_lt20_count{0};                   // Number of reads with mapq < 20
  double mapq_lt60_ratio{0};                // Ratio of reads with mapq < 60
  double mapq_lt40_ratio{0};                // Ratio of reads with mapq < 40
  double mapq_lt30_ratio{0};                // Ratio of reads with mapq < 30
  double mapq_lt20_ratio{0};                // Ratio of reads with mapq < 20
  u32 baseq_sum{0};                         // Sum of baseqs among reads
  double baseq_mean{0};
  u32 baseq_lt20_count{0};  // Number of reads with baseq < 20
  double baseq_lt20_ratio{0};
  double mq_af{0};
  double bq_af{0};
  u64 distance_min{0};              // Minimum distance to alignment end among reads
  u64 distance_max{0};              // Maximum distance to alignment end among reads
  u64 distance_sum{0};              // Sum of distances to alignment end among reads
  u64 distance_sum_lowbq{0};        // Sum of distances to alignment end among low-baseq reads
  u64 distance_sum_simplex{0};      // Sum of distances to alignment end among simplex reads
  double distance_mean{0};          // Mean of distances to alignment end among reads
  double distance_mean_lowbq{0};    // Mean of distances to alignment end among low-baseq reads
  double distance_mean_simplex{0};  // Mean of distances to alignment end among simplex reads
  double duplex_lowbq{0};           // Number of low baseq duplex reads
  u32 simplex{0};                   // Number of simplex reads
  double duplex_af{0};
  double duplex_dp{0};
  u32 familysize_sum{0};  // Sum of family sizes of reads
  double familysize_mean{0};
  u32 familysize_lt3_count{0};               // Number of reads with family size < 3
  u32 familysize_lt5_count{0};               // Number of reads with family size < 5
  double familysize_lt3_ratio{0};            // Number of reads with family size < 3
  double familysize_lt5_ratio{0};            // Number of reads with family size < 5
  u8 nonhomopolymer_mapq_min{0};             // Minimum mapq among reads not ending in homopolymers
  u8 nonhomopolymer_mapq_max{0};             // Maximum mapq among reads not ending in homopolymers
  u32 nonhomopolymer_mapq_sum{0};            // Sum of mapqs among reads not ending in homopolymers
  u32 nonhomopolymer_mapq_lt60_count{0};     // Number of reads not ending in homopolymers with mapq < 60
  u32 nonhomopolymer_mapq_lt40_count{0};     // Number of reads not ending in homopolymers with mapq < 40
  u32 nonhomopolymer_mapq_lt30_count{0};     // Number of reads not ending in homopolymers with mapq < 30
  u32 nonhomopolymer_mapq_lt20_count{0};     // Number of reads not ending in homopolymers with mapq < 20
  double nonhomopolymer_mapq_lt60_ratio{0};  // Ratio of reads not ending in homopolymers with mapq < 60
  double nonhomopolymer_mapq_lt40_ratio{0};  // Ratio of reads not ending in homopolymers with mapq < 40
  double nonhomopolymer_mapq_lt30_ratio{0};  // Ratio of reads not ending in homopolymers with mapq < 30
  double nonhomopolymer_mapq_lt20_ratio{0};  // Ratio of reads not ending in homopolymers with mapq < 20
  u8 nonhomopolymer_baseq_min{0};            // Minimum baseq among reads not ending in homopolymers
  u8 nonhomopolymer_baseq_max{0};            // Maximum baseq among reads not ending in homopolymers
  u32 nonhomopolymer_baseq_sum{0};           // Sum of baseqs among reads not ending in homopolymers
  double nonhomopolymer_baseq_mean{0};
  double nonhomopolymer_mapq_mean{0};
  u32 tumor_support{0};   // Count of reads in tumor sample
  u32 normal_support{0};  // Count of reads in normal sample
  u32 num_alt{0};
};

// Zero-initialized UnifiedReferenceFeature, this is used when no reference feature is found
static const UnifiedReferenceFeature kZeroUnifiedReferenceFeature{};

using UnifiedReferenceFeatures = std::unordered_map<u64, UnifiedReferenceFeature>;

// unified features set for variant (ALT) supporting reads
struct UnifiedVariantFeature {
  ReadIds read_ids{};
  double weighted_depth{0};     // Sum of (baseq/138)*(mapq/60) across alt-supporting reads
  u32 support{0};               // Count of reads
  u8 mapq_min{0};               // Minimum mapq
  u8 mapq_max{0};               // Maximum mapq
  u32 mapq_sum{0};              // Sum of mapqs
  u32 mapq_sum_lowbq{0};        // Sum of mapq in low-baseq reads
  u32 mapq_sum_simplex{0};      // Sum of mapq in simplex reads
  double mapq_mean{0};          // Mean of mapqs
  double mapq_mean_lowbq{0};    // Mean of mapq in low-baseq reads
  double mapq_mean_simplex{0};  // Mean of mapq in simplex reads
  u32 mapq_lt60_count{0};       // Number of reads with mapq < 60
  double mapq_lt60_ratio{0};
  u32 mapq_lt40_count{0};  // Number of reads with mapq < 40
  double mapq_lt40_ratio{0};
  u32 mapq_lt30_count{0};  // Number of reads with mapq < 30
  double mapq_lt30_ratio{0};
  u32 mapq_lt20_count{0};  // Number of reads with mapq < 20
  double mapq_lt20_ratio{0};
  double baseq_min{0};      // Minimum (mean) baseq
  double baseq_max{0};      // Maximum (mean) baseq
  double baseq_sum{0};      // Sum of (mean) baseqs
  double baseq_mean{0};     // Mean of baseqs
  u32 baseq_lt20_count{0};  // Number of reads with baseq < 20
  double baseq_lt20_ratio{0};
  u64 distance_min{0};              // Minimum distance to alignment end among reads
  u64 distance_max{0};              // Maximum distance to alignment end among reads
  u64 distance_sum{0};              // Sum of distances to alignment end among reads
  u64 distance_sum_lowbq{0};        // Sum of distances to alignment end among low-baseq reads
  u64 distance_sum_simplex{0};      // Sum of distances to alignment end among simplex reads
  double distance_mean{0};          // Mean of distances to alignment end among reads
  double distance_mean_lowbq{0};    // Mean of distances to alignment end among low-baseq reads
  double distance_mean_simplex{0};  // Mean of distances to alignment end among simplex reads
  double mq_af{0};
  double bq_af{0};
  double duplex_af{0};
  u32 familysize_sum{0};                     // sum of family sizes of reads
  u32 familysize_lt3_count{0};               // Number of reads with family size < 3
  u32 familysize_lt5_count{0};               // Number of reads with family size < 5
  double familysize_lt3_ratio{0};            // Number of reads with family size < 3
  double familysize_lt5_ratio{0};            // Number of reads with family size < 5
  double familysize_mean{0};                 // Mean of read family sizes
  u32 nonduplex{0};                          // Number of nonduplex reads
  u32 duplex{0};                             // Number of duplex reads
  double duplex_lowbq{0};                    // Number of low-baseq duplex reads
  u32 simplex{0};                            // Number of simplex reads
  u32 plusonly{0};                           // Number of plus-only reads
  u32 minusonly{0};                          // Number of minus-only reads
  double weighted_score{0};                  // Weighed score of duplex and non-duplex counts
  double strandbias{0};                      // Strand bias feature of variant
  std::string context{};                     // Sequence context of variant
  u32 context_index{0};                      // Context Index of sequence context
  u32 tumor_support{0};                      // total variant molecule count for tumor reads only
  u64 tumor_distance_sum{0};                 // Sum of distances from read end of variant site among tumor reads only
  double tumor_distance_mean{0};             // Mean of distances from read end of variant site among tumor reads only
  double tumor_baseq_sum{0};                 // Sum of baseq among tumor reads only
  double tumor_baseq_mean{0};                // Mean of baseq among tumor reads only
  u32 tumor_mapq_sum{0};                     // Sum of mapq among tumor reads only
  double tumor_mapq_mean{0};                 // Mean of mapq among tumor reads only
  u32 normal_support{0};                     // total variant molecule count for normal reads only
  u64 normal_distance_sum{0};                // Sum of distances from read end of variant site among normal reads only
  double normal_distance_mean{0};            // Mean of distances from read end of variant site among normal reads only
  double normal_baseq_sum{0};                // Sum of baseq among normal reads only
  double normal_baseq_mean{0};               // Mean of baseq among normal reads only
  u32 normal_mapq_sum{0};                    // Sum of mapq among normal reads only
  double normal_mapq_mean{0};                // Mean of mapq among normal reads only
  double tumor_af{0};                        // Ratio of tumor support to sum of tumor and normal support
  double normal_af{0};                       // Ratio of normal support to sum of tumor and normal support
  double rat{0};                             // Ratio of tumor AF to sum of tumor and normal AF
  u32 support_reverse{0};                    // Number of supporting reads aligned in the reverse orientation
  u32 tumor_support_reverse{0};              // Number of tumor supporting reads aligned in the reverse orientation
  u32 normal_support_reverse{0};             // Number of normal supporting reads aligned in the reverse orientation
  double alignmentbias{0};                   // Alignment bias for supporting reads
  double tumor_alignmentbias{0};             // Alignment bias for tumor supporting reads
  double normal_alignmentbias{0};            // Alignment bias for normal supporting reads
  double ml_score{0};                        // ML model score for variant
  std::vector<std::string> filter_status{};  // List of failure reasons for variant when filtered
  double adt{0};
  double adtl{0};
  double indel_af{0};  // TODO : Generalize this feature name? Its calculation is not specific to indels.

  static int ContextIndex(const std::string& context);
};

using UnifiedVariantFeatures = std::map<VariantId, UnifiedVariantFeature>;

/**
 * Builds a map between UnifiedFeatureCols enums to strings representing column names
 * @return a map between UnifiedFeatureCols and string
 */
std::unordered_map<UnifiedFeatureCols, std::string> BuildColToString();

/**
 * Builds a map between strings ( column names ) and UnifiedFeatureCols enums
 * @return a map between sting and UnifiedFeatureCols
 */
StrUnorderedMap<UnifiedFeatureCols> BuildStringToCol();

extern const std::unordered_map<UnifiedFeatureCols, std::string> kUnifiedFeatureColsToString;
extern const StrUnorderedMap<UnifiedFeatureCols> kStringToUnifiedFeatureCols;
bool IsVcfFeatureCol(UnifiedFeatureCols col);

UnifiedFeatureCols GetCol(const std::string& feature_name);
std::string GetFeatureName(UnifiedFeatureCols col);

// TODO: review the usage of these maps and consider renaming them for clarity
using PositionToVariantInfosMap = std::unordered_map<u64, UnifiedVariantFeatures>;
using ChromToVariantInfoMap = StrUnorderedMap<PositionToVariantInfosMap>;
using RefInfoMap = StrUnorderedMap<std::unordered_map<u64, UnifiedReferenceFeature>>;

using PositionToVariantInfosMapWithLabel =
    std::unordered_map<u64, std::map<VariantId, std::vector<std::tuple<UnifiedVariantFeature, std::string, u32>>>>;
using ChromToVariantInfoMapWithLabel = StrUnorderedMap<PositionToVariantInfosMapWithLabel>;
using RefInfoMapMultiSample = StrUnorderedMap<std::unordered_map<u64, UnifiedReferenceFeatures>>;

using PositionToVcfFeaturesMap = std::unordered_map<u64, std::map<VariantId, VcfFeature>>;
using ChromToVcfFeaturesMap = StrUnorderedMap<PositionToVcfFeaturesMap>;

using PositionToVcfFeaturesMapMultiSample =
    std::unordered_map<u64, std::map<VariantId, std::vector<std::pair<VcfFeature, u32>>>>;
using ChromToVcfFeaturesMapMultiSample = StrUnorderedMap<PositionToVcfFeaturesMapMultiSample>;
}  // namespace xoos::svc
