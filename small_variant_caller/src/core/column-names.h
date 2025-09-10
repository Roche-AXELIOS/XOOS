#pragma once

#include <string>

#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>

namespace xoos::svc {

/**
 * Synopsis:
 * Consolidate all feature column names here to simplify future changes and to avoid spelling mistakes.
 */

// TODO: review the usage of column names and consider renaming them for clarity
// TODO: remove `ref_int` and `alt_int` and simply use `ref` and `alt` instead?

// variant ID
static const std::string kNameChrom = "chrom";
static const std::string kNamePos = "pos";
static const std::string kNameRef = "ref";
static const std::string kNameAlt = "alt";
static const std::string kNameAltLen = "alt_len";
static const std::string kNameSubtypeIndex = "subtype_index";
static const std::string kNameVariantType = "variant_type";

// unified BAM features
static const std::string kNameBaseqAF = "baseq_af";
static const std::string kNameBaseqMin = "baseq_min";
static const std::string kNameBaseqMax = "baseq_max";
static const std::string kNameBaseqSum = "baseq_sum";
static const std::string kNameBaseqMean = "baseq_mean";
static const std::string kNameBaseqLT20Count = "baseq_lt20_count";
static const std::string kNameBaseqLT20Ratio = "baseq_lt20_ratio";
static const std::string kNameContext = "context";
static const std::string kNameContextIndex = "context_index";
static const std::string kNameDistanceMin = "distance_min";
static const std::string kNameDistanceMax = "distance_max";
static const std::string kNameDistanceSum = "distance_sum";
static const std::string kNameDistanceSumLowbq = "distance_sum_lowbq";
static const std::string kNameDistanceSumSimplex = "distance_sum_simplex";
static const std::string kNameDistanceMean = "distance_mean";
static const std::string kNameDistanceMeanLowbq = "distance_mean_lowbq";
static const std::string kNameDistanceMeanSimplex = "distance_mean_simplex";
static const std::string kNameDuplex = "duplex";
static const std::string kNameDuplexLowbq = "duplex_lowbq";
static const std::string kNameDuplexAF = "duplex_af";
static const std::string kNameDuplexDP = "duplex_dp";
static const std::string kNameSimplex = "simplex";
static const std::string kNameFamilysizeSum = "familysize_sum";
static const std::string kNameFamilysizeMean = "familysize_mean";
static const std::string kNameFamilysizeLT3Count = "familysize_lt3_count";
static const std::string kNameFamilysizeLT5Count = "familysize_lt5_count";
static const std::string kNameFamilysizeLT3Ratio = "familysize_lt3_ratio";
static const std::string kNameFamilysizeLT5Ratio = "familysize_lt5_ratio";
static const std::string kNameMapqAF = "mapq_af";
static const std::string kNameMapqMin = "mapq_min";
static const std::string kNameMapqMax = "mapq_max";
static const std::string kNameMapqSum = "mapq_sum";
static const std::string kNameMapqSumLowbq = "mapq_sum_lowbq";
static const std::string kNameMapqSumSimplex = "mapq_sum_simplex";
static const std::string kNameMapqMean = "mapq_mean";
static const std::string kNameMapqMeanLowbq = "mapq_mean_lowbq";
static const std::string kNameMapqMeanSimplex = "mapq_mean_simplex";
static const std::string kNameMapqLT60Count = "mapq_lt60_count";
static const std::string kNameMapqLT40Count = "mapq_lt40_count";
static const std::string kNameMapqLT30Count = "mapq_lt30_count";
static const std::string kNameMapqLT20Count = "mapq_lt20_count";
static const std::string kNameMapqLT60Ratio = "mapq_lt60_ratio";
static const std::string kNameMapqLT40Ratio = "mapq_lt40_ratio";
static const std::string kNameMapqLT30Ratio = "mapq_lt30_ratio";
static const std::string kNameMapqLT20Ratio = "mapq_lt20_ratio";
static const std::string kNameMinusOnly = "minusonly";
static const std::string kNameNonDuplex = "nonduplex";
static const std::string kNameNormalSupport = "normal_support";
static const std::string kNameNormalMapqMean = "normal_mapq_mean";
static const std::string kNameNormalBaseqMean = "normal_baseq_mean";
static const std::string kNameNormalDistanceMean = "normal_distance_mean";
static const std::string kNameNormalRefSupport = "normal_ref_support";
static const std::string kNamePlusOnly = "plusonly";
static const std::string kNameRefBaseqSum = "ref_baseq_sum";
static const std::string kNameRefBaseqMean = "ref_baseq_mean";
static const std::string kNameRefBaseqLT20Count = "ref_baseq_lt20_count";
static const std::string kNameRefBaseqLT20Ratio = "ref_baseq_lt20_ratio";
static const std::string kNameRefBaseqAF = "ref_baseq_af";
static const std::string kNameRefDistanceMin = "ref_distance_min";
static const std::string kNameRefDistanceMax = "ref_distance_max";
static const std::string kNameRefDistanceSum = "ref_distance_sum";
static const std::string kNameRefDistanceSumLowbq = "ref_distance_sum_lowbq";
static const std::string kNameRefDistanceSumSimplex = "ref_distance_sum_simplex";
static const std::string kNameRefDistanceMean = "ref_distance_mean";
static const std::string kNameRefDistanceMeanLowbq = "ref_distance_mean_lowbq";
static const std::string kNameRefDistanceMeanSimplex = "ref_distance_mean_simplex";
static const std::string kNameRefDuplexLowbq = "ref_duplex_lowbq";
static const std::string kNameRefSimplex = "ref_simplex";
static const std::string kNameRefDuplexAF = "ref_duplex_af";
static const std::string kNameRefFamilysizeSum = "ref_familysize_sum";
static const std::string kNameRefFamilysizeMean = "ref_familysize_mean";
static const std::string kNameRefFamilysizeLT3Count = "ref_familysize_lt3_count";
static const std::string kNameRefFamilysizeLT5Count = "ref_familysize_lt5_count";
static const std::string kNameRefFamilysizeLT3Ratio = "ref_familysize_lt3_ratio";
static const std::string kNameRefFamilysizeLT5Ratio = "ref_familysize_lt5_ratio";
static const std::string kNameRefMapqMin = "ref_mapq_min";
static const std::string kNameRefMapqMax = "ref_mapq_max";
static const std::string kNameRefMapqSum = "ref_mapq_sum";
static const std::string kNameRefMapqSumLowbq = "ref_mapq_sum_lowbq";
static const std::string kNameRefMapqSumSimplex = "ref_mapq_sum_simplex";
static const std::string kNameRefMapqMean = "ref_mapq_mean";
static const std::string kNameRefMapqMeanLowbq = "ref_mapq_mean_lowbq";
static const std::string kNameRefMapqMeanSimplex = "ref_mapq_mean_simplex";
static const std::string kNameRefMapqAF = "ref_mapq_af";
static const std::string kNameRefMapqLT60Count = "ref_mapq_lt60_count";
static const std::string kNameRefMapqLT40Count = "ref_mapq_lt40_count";
static const std::string kNameRefMapqLT30Count = "ref_mapq_lt30_count";
static const std::string kNameRefMapqLT20Count = "ref_mapq_lt20_count";
static const std::string kNameRefMapqLT60Ratio = "ref_mapq_lt60_ratio";
static const std::string kNameRefMapqLT40Ratio = "ref_mapq_lt40_ratio";
static const std::string kNameRefMapqLT30Ratio = "ref_mapq_lt30_ratio";
static const std::string kNameRefMapqLT20Ratio = "ref_mapq_lt20_ratio";
static const std::string kNameRefNonHpBaseqMin = "ref_nonhp_baseq_min";
static const std::string kNameRefNonHpBaseqMax = "ref_nonhp_baseq_max";
static const std::string kNameRefNonHpBaseqSum = "ref_nonhp_baseq_sum";
static const std::string kNameRefNonHpBaseqMean = "ref_nonhp_baseq_mean";
static const std::string kNameRefNonHpMapqMin = "ref_nonhp_mapq_min";
static const std::string kNameRefNonHpMapqMax = "ref_nonhp_mapq_max";
static const std::string kNameRefNonHpMapqSum = "ref_nonhp_mapq_sum";
static const std::string kNameRefNonHpMapqMean = "ref_nonhp_mapq_mean";
static const std::string kNameRefNonHpMapqLT60Count = "ref_nonhp_mapq_lt60_count";
static const std::string kNameRefNonHpMapqLT40Count = "ref_nonhp_mapq_lt40_count";
static const std::string kNameRefNonHpMapqLT30Count = "ref_nonhp_mapq_lt30_count";
static const std::string kNameRefNonHpMapqLT20Count = "ref_nonhp_mapq_lt20_count";
static const std::string kNameRefNonHpMapqLT60Ratio = "ref_nonhp_mapq_lt60_ratio";
static const std::string kNameRefNonHpMapqLT40Ratio = "ref_nonhp_mapq_lt40_ratio";
static const std::string kNameRefNonHpMapqLT30Ratio = "ref_nonhp_mapq_lt30_ratio";
static const std::string kNameRefNonHpMapqLT20Ratio = "ref_nonhp_mapq_lt20_ratio";
static const std::string kNameRefNonHpSupport = "ref_nonhp_support";
static const std::string kNameRefNonHpWeightedDepth = "ref_nonhp_weighted_depth";
static const std::string kNameRefSupport = "ref_support";
static const std::string kNameRefWeightedDepth = "ref_weighted_depth";
static const std::string kNameStrandBias = "strandbias";
static const std::string kNameSupport = "support";
static const std::string kNameTumorSupport = "tumor_support";
static const std::string kNameTumorMapqMean = "tumor_mapq_mean";
static const std::string kNameTumorBaseqMean = "tumor_baseq_mean";
static const std::string kNameTumorDistanceMean = "tumor_distance_mean";
static const std::string kNameTumorRefSupport = "tumor_ref_support";
static const std::string kNameBamTumorAF = "bam_tumor_af";
static const std::string kNameBamNormalAF = "bam_normal_af";
static const std::string kNameBamRAT = "bam_rat";
static const std::string kNameSupportReverse = "support_reverse";
static const std::string kNameTumorSupportReverse = "tumor_support_reverse";
static const std::string kNameNormalSupportReverse = "normal_support_reverse";
static const std::string kNameAlignmentBias = "alignmentbias";
static const std::string kNameTumorAlignmentBias = "tumor_alignmentbias";
static const std::string kNameNormalAlignmentBias = "normal_alignmentbias";
static const std::string kNameWeightedDepth = "weighted_depth";
static const std::string kNameWeightedScore = "weightedscore";

// VCF and reference sequence features
static const std::string kNameNalod = "nalod";
static const std::string kNameNlod = "nlod";
static const std::string kNameTlod = "tlod";
static const std::string kNameMpos = "mpos";
static const std::string kNameMmqRef = "mmq_ref";
static const std::string kNameMmqAlt = "mmq_alt";
static const std::string kNameMbqRef = "mbq_ref";
static const std::string kNameMbqAlt = "mbq_alt";
static const std::string kNamePre2BpContext = "pre_2bp_context";
static const std::string kNamePost2BpContext = "post_2bp_context";
static const std::string kNamePost30BpContext = "post_30bp_context";
static const std::string kNameUniq3mers = "uniq_3mers";
static const std::string kNameUniq4mers = "uniq_4mers";
static const std::string kNameUniq5mers = "uniq_5mers";
static const std::string kNameUniq6mers = "uniq_6mers";
static const std::string kNameHomopolymer = "homopolymer";
static const std::string kNameDirepeat = "direpeat";
static const std::string kNameTrirepeat = "trirepeat";
static const std::string kNameQuadrepeat = "quadrepeat";
static const std::string kNameVariantDensity = "variant_density";
static const std::string kNameRefAD = "ref_ad";
static const std::string kNameAltAD = "alt_ad";
static const std::string kNameAltAD2 = "alt_ad2";
static const std::string kNameRefAdAF = "ref_ad_af";
static const std::string kNameAltAdAF = "alt_ad_af";
static const std::string kNameAltAd2AF = "alt_ad2_af";
static const std::string kNameTumorAltAD = "tumor_alt_ad";
static const std::string kNameNormalAltAD = "normal_alt_ad";
static const std::string kNameTumorAF = "tumor_af";
static const std::string kNameNormalAF = "normal_af";
static const std::string kNameTumorNormalAfRatio = "tumor_normal_af_ratio";
static const std::string kNameTumorDP = "tumor_dp";
static const std::string kNameNormalDP = "normal_dp";
static const std::string kNamePopAF = "popaf";
static const std::string kNameGenotype = "genotype";
static const std::string kNameQual = "vcf_variant_qual";
static const std::string kNameGq = "vcf_variant_gq";
static const std::string kNameHapcomp = "hapcomp";
static const std::string kNameHapdom = "hapdom";
static const std::string kNameRu = "ru";
static const std::string kNameRpaRef = "rpa_ref";
static const std::string kNameRpaAlt = "rpa_alt";
static const std::string kNameStr = "str";
static const std::string kNameAtInterest = "at_interest";

// Computed Germline InDel Fields
static const std::string kNameADT = "adt";
static const std::string kNameADTL = "adtl";
static const std::string kNameNumAlt = "num_alt";
static const std::string kNameIndelAF = "indel_af";

// output fields
static const std::string kNameMLScore = "ml_score";
static const std::string kNameFilterStatus = "filter_status";

static const vec<std::string> kVariantIDFeatureNames{
    kNameChrom, kNamePos, kNameRef, kNameAlt, kNameAltLen, kNameSubtypeIndex, kNameVariantType};
static const vec<std::string> kBamFeatureNames{kNameADT,
                                               kNameADTL,
                                               kNameBaseqAF,
                                               kNameBaseqMin,
                                               kNameBaseqMax,
                                               kNameBaseqSum,
                                               kNameBaseqMean,
                                               kNameBaseqLT20Count,
                                               kNameBaseqLT20Ratio,
                                               kNameContext,
                                               kNameContextIndex,
                                               kNameDistanceMin,
                                               kNameDistanceMax,
                                               kNameDistanceSum,
                                               kNameDistanceSumLowbq,
                                               kNameDistanceSumSimplex,
                                               kNameDistanceMean,
                                               kNameDistanceMeanLowbq,
                                               kNameDistanceMeanSimplex,
                                               kNameDuplex,
                                               kNameDuplexLowbq,
                                               kNameSimplex,
                                               kNameDuplexDP,
                                               kNameDuplexAF,
                                               kNameFamilysizeSum,
                                               kNameFamilysizeMean,
                                               kNameFamilysizeLT3Count,
                                               kNameFamilysizeLT5Count,
                                               kNameFamilysizeLT3Ratio,
                                               kNameFamilysizeLT5Ratio,
                                               kNameIndelAF,
                                               kNameMapqAF,
                                               kNameMapqMin,
                                               kNameMapqMax,
                                               kNameMapqSum,
                                               kNameMapqSumLowbq,
                                               kNameMapqSumSimplex,
                                               kNameMapqMean,
                                               kNameMapqMeanLowbq,
                                               kNameMapqMeanSimplex,
                                               kNameMapqLT60Count,
                                               kNameMapqLT40Count,
                                               kNameMapqLT30Count,
                                               kNameMapqLT20Count,
                                               kNameMapqLT60Ratio,
                                               kNameMapqLT40Ratio,
                                               kNameMapqLT30Ratio,
                                               kNameMapqLT20Ratio,
                                               kNameMinusOnly,
                                               kNameNonDuplex,
                                               kNameNormalSupport,
                                               kNameNormalMapqMean,
                                               kNameNormalBaseqMean,
                                               kNameNormalDistanceMean,
                                               kNameNormalRefSupport,
                                               kNameNumAlt,
                                               kNamePlusOnly,
                                               kNameRefBaseqAF,
                                               kNameRefBaseqSum,
                                               kNameRefBaseqMean,
                                               kNameRefBaseqLT20Count,
                                               kNameRefBaseqLT20Ratio,
                                               kNameRefDistanceMin,
                                               kNameRefDistanceMax,
                                               kNameRefDistanceSum,
                                               kNameRefDistanceSumLowbq,
                                               kNameRefDistanceSumSimplex,
                                               kNameRefDistanceMean,
                                               kNameRefDistanceMeanLowbq,
                                               kNameRefDistanceMeanSimplex,
                                               kNameRefDuplexLowbq,
                                               kNameRefSimplex,
                                               kNameRefDuplexAF,
                                               kNameRefFamilysizeSum,
                                               kNameRefFamilysizeMean,
                                               kNameRefFamilysizeLT3Count,
                                               kNameRefFamilysizeLT5Count,
                                               kNameRefFamilysizeLT3Ratio,
                                               kNameRefFamilysizeLT5Ratio,
                                               kNameRefMapqAF,
                                               kNameRefMapqMin,
                                               kNameRefMapqMax,
                                               kNameRefMapqSum,
                                               kNameRefMapqSumLowbq,
                                               kNameRefMapqSumSimplex,
                                               kNameRefMapqMean,
                                               kNameRefMapqMeanLowbq,
                                               kNameRefMapqMeanSimplex,
                                               kNameRefMapqLT60Count,
                                               kNameRefMapqLT40Count,
                                               kNameRefMapqLT30Count,
                                               kNameRefMapqLT20Count,
                                               kNameRefMapqLT60Ratio,
                                               kNameRefMapqLT40Ratio,
                                               kNameRefMapqLT30Ratio,
                                               kNameRefMapqLT20Ratio,
                                               kNameRefNonHpBaseqMin,
                                               kNameRefNonHpBaseqMax,
                                               kNameRefNonHpBaseqSum,
                                               kNameRefNonHpBaseqMean,
                                               kNameRefNonHpMapqMin,
                                               kNameRefNonHpMapqMax,
                                               kNameRefNonHpMapqSum,
                                               kNameRefNonHpMapqMean,
                                               kNameRefNonHpMapqLT60Count,
                                               kNameRefNonHpMapqLT40Count,
                                               kNameRefNonHpMapqLT30Count,
                                               kNameRefNonHpMapqLT20Count,
                                               kNameRefNonHpMapqLT60Ratio,
                                               kNameRefNonHpMapqLT40Ratio,
                                               kNameRefNonHpMapqLT30Ratio,
                                               kNameRefNonHpMapqLT20Ratio,
                                               kNameRefNonHpSupport,
                                               kNameRefNonHpWeightedDepth,
                                               kNameRefSupport,
                                               kNameRefWeightedDepth,
                                               kNameStrandBias,
                                               kNameSupport,
                                               kNameTumorSupport,
                                               kNameTumorMapqMean,
                                               kNameTumorBaseqMean,
                                               kNameTumorDistanceMean,
                                               kNameTumorRefSupport,
                                               kNameSupportReverse,
                                               kNameTumorSupportReverse,
                                               kNameNormalSupportReverse,
                                               kNameBamTumorAF,
                                               kNameBamNormalAF,
                                               kNameBamRAT,
                                               kNameAlignmentBias,
                                               kNameTumorAlignmentBias,
                                               kNameNormalAlignmentBias,
                                               kNameWeightedDepth,
                                               kNameWeightedScore};
static const vec<std::string> kVcfFeatureNames{kNameNalod,
                                               kNameNlod,
                                               kNameTlod,
                                               kNameMpos,
                                               kNameMmqRef,
                                               kNameMmqAlt,
                                               kNameMbqRef,
                                               kNameMbqAlt,
                                               kNamePre2BpContext,
                                               kNamePost2BpContext,
                                               kNamePost30BpContext,
                                               kNameUniq3mers,
                                               kNameUniq4mers,
                                               kNameUniq5mers,
                                               kNameUniq6mers,
                                               kNameHomopolymer,
                                               kNameDirepeat,
                                               kNameTrirepeat,
                                               kNameQuadrepeat,
                                               kNameVariantDensity,
                                               kNameRefAD,
                                               kNameAltAD,
                                               kNameAltAD2,
                                               kNameRefAdAF,
                                               kNameAltAdAF,
                                               kNameAltAd2AF,
                                               kNameTumorAltAD,
                                               kNameNormalAltAD,
                                               kNameTumorAF,
                                               kNameNormalAF,
                                               kNameTumorNormalAfRatio,
                                               kNameTumorDP,
                                               kNameNormalDP,
                                               kNamePopAF,
                                               kNameGenotype,
                                               kNameQual,
                                               kNameGq,
                                               kNameHapcomp,
                                               kNameHapdom,
                                               kNameRu,
                                               kNameRpaRef,
                                               kNameRpaAlt,
                                               kNameStr,
                                               kNameAtInterest,
                                               kNameADT,
                                               kNameADTL,
                                               kNameNumAlt,
                                               kNameIndelAF};

vec<std::string> FindUnsupportedFeatureNames(const vec<std::string>& names,
                                             const StrSet& supported_names,
                                             bool allow_num_suffix = false);
vec<std::string> FindUnsupportedBamFeatureNames(const vec<std::string>& names, bool allow_num_suffix = false);
vec<std::string> FindUnsupportedVcfFeatureNames(const vec<std::string>& names, bool allow_num_suffix = false);
vec<std::string> FindUnsupportedScoringFeatureNames(const vec<std::string>& names, bool allow_num_suffix = false);

}  // namespace xoos::svc
