#pragma once

#include <array>
#include <string>
#include <string_view>

#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>

namespace xoos::svc {

/**
 * Synopsis:
 * Consolidate all feature column names here to simplify future changes and to avoid spelling mistakes.
 */

// feature name prefixes for tumor and normal samples
inline constexpr std::string_view kTumorPrefix = "tumor_";
inline constexpr std::string_view kNormalPrefix = "normal_";

// variant ID
inline constexpr std::string_view kNameChrom = "chrom";
inline constexpr std::string_view kNamePos = "pos";
inline constexpr std::string_view kNameRef = "ref";
inline constexpr std::string_view kNameAlt = "alt";
inline constexpr std::string_view kNameAltLen = "alt_len";
inline constexpr std::string_view kNameSubtypeIndex = "subtype_index";
inline constexpr std::string_view kNameVariantType = "variant_type";

// unified BAM features
inline constexpr std::string_view kNameBaseqAF = "baseq_af";
inline constexpr std::string_view kNameBaseqMin = "baseq_min";
inline constexpr std::string_view kNameBaseqMax = "baseq_max";
inline constexpr std::string_view kNameBaseqSum = "baseq_sum";
inline constexpr std::string_view kNameBaseqMean = "baseq_mean";
inline constexpr std::string_view kNameBaseqLT20Count = "baseq_lt20_count";
inline constexpr std::string_view kNameBaseqLT20Ratio = "baseq_lt20_ratio";
inline constexpr std::string_view kNameContext = "context";
inline constexpr std::string_view kNameContextIndex = "context_index";
inline constexpr std::string_view kNameDistanceMin = "distance_min";
inline constexpr std::string_view kNameDistanceMax = "distance_max";
inline constexpr std::string_view kNameDistanceSum = "distance_sum";
inline constexpr std::string_view kNameDistanceSumLowbq = "distance_sum_lowbq";
inline constexpr std::string_view kNameDistanceSumSimplex = "distance_sum_simplex";
inline constexpr std::string_view kNameDistanceMean = "distance_mean";
inline constexpr std::string_view kNameDistanceMeanLowbq = "distance_mean_lowbq";
inline constexpr std::string_view kNameDistanceMeanSimplex = "distance_mean_simplex";
inline constexpr std::string_view kNameDuplex = "duplex";
inline constexpr std::string_view kNameDuplexLowbq = "duplex_lowbq";
inline constexpr std::string_view kNameDuplexConcordant = "duplex_concordant";
inline constexpr std::string_view kNameDuplexSimplex = "duplex_simplex";
inline constexpr std::string_view kNameDuplexDiscordant = "duplex_discordant";
inline constexpr std::string_view kNameDuplexAF = "duplex_af";
inline constexpr std::string_view kNameDuplexDP = "duplex_dp";
inline constexpr std::string_view kNameDuplexOther = "duplex_other";
inline constexpr std::string_view kNameDuplexOtherAF = "duplex_other_af";
inline constexpr std::string_view kNameDuplexGQ = "duplex_gq";
inline constexpr std::string_view kNameAfqGT01 = "afq_gt01";
inline constexpr std::string_view kNameAfqGT11 = "afq_gt11";
inline constexpr std::string_view kNameAfqGT12 = "afq_gt12";
inline constexpr std::string_view kNameAfqGT = "afq_gt";
inline constexpr std::string_view kNameSimplex = "simplex";
inline constexpr std::string_view kNameFamilysizeSum = "familysize_sum";
inline constexpr std::string_view kNameFamilysizeMean = "familysize_mean";
inline constexpr std::string_view kNameFamilysizeLT3Count = "familysize_lt3_count";
inline constexpr std::string_view kNameFamilysizeLT5Count = "familysize_lt5_count";
inline constexpr std::string_view kNameFamilysizeLT3Ratio = "familysize_lt3_ratio";
inline constexpr std::string_view kNameFamilysizeLT5Ratio = "familysize_lt5_ratio";
inline constexpr std::string_view kNameMapqAF = "mapq_af";
inline constexpr std::string_view kNameMapqMin = "mapq_min";
inline constexpr std::string_view kNameMapqMax = "mapq_max";
inline constexpr std::string_view kNameMapqSum = "mapq_sum";
inline constexpr std::string_view kNameMapqSumLowbq = "mapq_sum_lowbq";
inline constexpr std::string_view kNameMapqSumSimplex = "mapq_sum_simplex";
inline constexpr std::string_view kNameMapqMean = "mapq_mean";
inline constexpr std::string_view kNameMapqMeanLowbq = "mapq_mean_lowbq";
inline constexpr std::string_view kNameMapqMeanSimplex = "mapq_mean_simplex";
inline constexpr std::string_view kNameMapqLT60Count = "mapq_lt60_count";
inline constexpr std::string_view kNameMapqLT40Count = "mapq_lt40_count";
inline constexpr std::string_view kNameMapqLT30Count = "mapq_lt30_count";
inline constexpr std::string_view kNameMapqLT20Count = "mapq_lt20_count";
inline constexpr std::string_view kNameMapqLT60Ratio = "mapq_lt60_ratio";
inline constexpr std::string_view kNameMapqLT40Ratio = "mapq_lt40_ratio";
inline constexpr std::string_view kNameMapqLT30Ratio = "mapq_lt30_ratio";
inline constexpr std::string_view kNameMapqLT20Ratio = "mapq_lt20_ratio";
inline constexpr std::string_view kNameMinusOnly = "minusonly";
inline constexpr std::string_view kNameNonDuplex = "nonduplex";
inline constexpr std::string_view kNamePlusOnly = "plusonly";
inline constexpr std::string_view kNameRefBaseqSum = "ref_baseq_sum";
inline constexpr std::string_view kNameRefBaseqMean = "ref_baseq_mean";
inline constexpr std::string_view kNameRefBaseqLT20Count = "ref_baseq_lt20_count";
inline constexpr std::string_view kNameRefBaseqLT20Ratio = "ref_baseq_lt20_ratio";
inline constexpr std::string_view kNameRefBaseqAF = "ref_baseq_af";
inline constexpr std::string_view kNameRefDistanceMin = "ref_distance_min";
inline constexpr std::string_view kNameRefDistanceMax = "ref_distance_max";
inline constexpr std::string_view kNameRefDistanceSum = "ref_distance_sum";
inline constexpr std::string_view kNameRefDistanceSumLowbq = "ref_distance_sum_lowbq";
inline constexpr std::string_view kNameRefDistanceSumSimplex = "ref_distance_sum_simplex";
inline constexpr std::string_view kNameRefDistanceMean = "ref_distance_mean";
inline constexpr std::string_view kNameRefDistanceMeanLowbq = "ref_distance_mean_lowbq";
inline constexpr std::string_view kNameRefDistanceMeanSimplex = "ref_distance_mean_simplex";
inline constexpr std::string_view kNameRefDuplexLowbq = "ref_duplex_lowbq";
inline constexpr std::string_view kNameRefDuplexConcordant = "ref_duplex_concordant";
inline constexpr std::string_view kNameRefDuplexSimplex = "ref_duplex_simplex";
inline constexpr std::string_view kNameRefDuplexDiscordant = "ref_duplex_discordant";
inline constexpr std::string_view kNameRefSimplex = "ref_simplex";
inline constexpr std::string_view kNameRefDuplexAF = "ref_duplex_af";
inline constexpr std::string_view kNameRefFamilysizeSum = "ref_familysize_sum";
inline constexpr std::string_view kNameRefFamilysizeMean = "ref_familysize_mean";
inline constexpr std::string_view kNameRefFamilysizeLT3Count = "ref_familysize_lt3_count";
inline constexpr std::string_view kNameRefFamilysizeLT5Count = "ref_familysize_lt5_count";
inline constexpr std::string_view kNameRefFamilysizeLT3Ratio = "ref_familysize_lt3_ratio";
inline constexpr std::string_view kNameRefFamilysizeLT5Ratio = "ref_familysize_lt5_ratio";
inline constexpr std::string_view kNameRefMapqMin = "ref_mapq_min";
inline constexpr std::string_view kNameRefMapqMax = "ref_mapq_max";
inline constexpr std::string_view kNameRefMapqSum = "ref_mapq_sum";
inline constexpr std::string_view kNameRefMapqSumLowbq = "ref_mapq_sum_lowbq";
inline constexpr std::string_view kNameRefMapqSumSimplex = "ref_mapq_sum_simplex";
inline constexpr std::string_view kNameRefMapqMean = "ref_mapq_mean";
inline constexpr std::string_view kNameRefMapqMeanLowbq = "ref_mapq_mean_lowbq";
inline constexpr std::string_view kNameRefMapqMeanSimplex = "ref_mapq_mean_simplex";
inline constexpr std::string_view kNameRefMapqAF = "ref_mapq_af";
inline constexpr std::string_view kNameRefMapqLT60Count = "ref_mapq_lt60_count";
inline constexpr std::string_view kNameRefMapqLT40Count = "ref_mapq_lt40_count";
inline constexpr std::string_view kNameRefMapqLT30Count = "ref_mapq_lt30_count";
inline constexpr std::string_view kNameRefMapqLT20Count = "ref_mapq_lt20_count";
inline constexpr std::string_view kNameRefMapqLT60Ratio = "ref_mapq_lt60_ratio";
inline constexpr std::string_view kNameRefMapqLT40Ratio = "ref_mapq_lt40_ratio";
inline constexpr std::string_view kNameRefMapqLT30Ratio = "ref_mapq_lt30_ratio";
inline constexpr std::string_view kNameRefMapqLT20Ratio = "ref_mapq_lt20_ratio";
inline constexpr std::string_view kNameRefNonHpBaseqMin = "ref_nonhp_baseq_min";
inline constexpr std::string_view kNameRefNonHpBaseqMax = "ref_nonhp_baseq_max";
inline constexpr std::string_view kNameRefNonHpBaseqSum = "ref_nonhp_baseq_sum";
inline constexpr std::string_view kNameRefNonHpBaseqMean = "ref_nonhp_baseq_mean";
inline constexpr std::string_view kNameRefNonHpMapqMin = "ref_nonhp_mapq_min";
inline constexpr std::string_view kNameRefNonHpMapqMax = "ref_nonhp_mapq_max";
inline constexpr std::string_view kNameRefNonHpMapqSum = "ref_nonhp_mapq_sum";
inline constexpr std::string_view kNameRefNonHpMapqMean = "ref_nonhp_mapq_mean";
inline constexpr std::string_view kNameRefNonHpMapqLT60Count = "ref_nonhp_mapq_lt60_count";
inline constexpr std::string_view kNameRefNonHpMapqLT40Count = "ref_nonhp_mapq_lt40_count";
inline constexpr std::string_view kNameRefNonHpMapqLT30Count = "ref_nonhp_mapq_lt30_count";
inline constexpr std::string_view kNameRefNonHpMapqLT20Count = "ref_nonhp_mapq_lt20_count";
inline constexpr std::string_view kNameRefNonHpMapqLT60Ratio = "ref_nonhp_mapq_lt60_ratio";
inline constexpr std::string_view kNameRefNonHpMapqLT40Ratio = "ref_nonhp_mapq_lt40_ratio";
inline constexpr std::string_view kNameRefNonHpMapqLT30Ratio = "ref_nonhp_mapq_lt30_ratio";
inline constexpr std::string_view kNameRefNonHpMapqLT20Ratio = "ref_nonhp_mapq_lt20_ratio";
inline constexpr std::string_view kNameRefNonHpSupport = "ref_nonhp_support";
inline constexpr std::string_view kNameRefNonHpWeightedDepth = "ref_nonhp_weighted_depth";
inline constexpr std::string_view kNameRefSupport = "ref_support";
inline constexpr std::string_view kNameRefWeightedDepth = "ref_weighted_depth";
inline constexpr std::string_view kNameStrandBias = "strandbias";
inline constexpr std::string_view kNameSupport = "support";
inline constexpr std::string_view kNameSupportReverse = "support_reverse";
inline constexpr std::string_view kNameAlignmentBias = "alignmentbias";
inline constexpr std::string_view kNameWeightedDepth = "weighted_depth";
inline constexpr std::string_view kNameWeightedScore = "weightedscore";
// Features to reduce multi-allelic false positives in tumor-normal-wgs
inline constexpr std::string_view kNameDel1bpAF = "del_1bp_af";
inline constexpr std::string_view kNameDel2bpAF = "del_2bp_af";
inline constexpr std::string_view kNameDel3bpAF = "del_3bp_af";
inline constexpr std::string_view kNameDelGt3bpAF = "del_gt3bp_af";
inline constexpr std::string_view kNameIns1bpAF = "ins_1bp_af";
inline constexpr std::string_view kNameIns2bpAF = "ins_2bp_af";
inline constexpr std::string_view kNameIns3bpAF = "ins_3bp_af";
inline constexpr std::string_view kNameInsGt3bpAF = "ins_gt3bp_af";
inline constexpr std::string_view kNameDel1bp2bpAFLogRatio = "del_1bp_2bp_af_logratio";
inline constexpr std::string_view kNameDel2bp3bpAFLogRatio = "del_2bp_3bp_af_logratio";
inline constexpr std::string_view kNameDel3bpGt3bpAFLogRatio = "del_3bp_gt3bp_af_logratio";
inline constexpr std::string_view kNameIns1bp2bpAFLogRatio = "ins_1bp_2bp_af_logratio";
inline constexpr std::string_view kNameIns2bp3bpAFLogRatio = "ins_2bp_3bp_af_logratio";
inline constexpr std::string_view kNameIns3bpGt3bpAFLogRatio = "ins_3bp_gt3bp_af_logratio";
inline constexpr std::string_view kNameDel1bpIns1bpAFLogRatio = "del_1bp_ins_1bp_af_logratio";
inline constexpr std::string_view kNameSnvAF = "snv_af";
inline constexpr std::string_view kNameIsInRepeat = "is_in_repeat";
inline constexpr std::string_view kNameRepeatLen = "repeat_len";
inline constexpr std::string_view kNameRepeatUnitLen = "repeat_unit_len";
inline constexpr std::string_view kNameIsUnitSizeMatch = "is_unit_size_match";

// tumor-normal specific BAM features
inline constexpr std::string_view kNameBamTnAfRatio = "bam_tn_af_ratio";

// VCF and reference sequence features
inline constexpr std::string_view kNameNalod = "nalod";
inline constexpr std::string_view kNameNlod = "nlod";
inline constexpr std::string_view kNameTlod = "tlod";
inline constexpr std::string_view kNameMpos = "mpos";
inline constexpr std::string_view kNameMmqRef = "mmq_ref";
inline constexpr std::string_view kNameMmqAlt = "mmq_alt";
inline constexpr std::string_view kNameMbqRef = "mbq_ref";
inline constexpr std::string_view kNameMbqAlt = "mbq_alt";
inline constexpr std::string_view kNamePre2BpContext = "pre_2bp_context";
inline constexpr std::string_view kNamePost2BpContext = "post_2bp_context";
inline constexpr std::string_view kNamePost30BpContext = "post_30bp_context";
inline constexpr std::string_view kNameUniq3mers = "uniq_3mers";
inline constexpr std::string_view kNameUniq4mers = "uniq_4mers";
inline constexpr std::string_view kNameUniq5mers = "uniq_5mers";
inline constexpr std::string_view kNameUniq6mers = "uniq_6mers";
inline constexpr std::string_view kNameHomopolymer = "homopolymer";
inline constexpr std::string_view kNameDirepeat = "direpeat";
inline constexpr std::string_view kNameTrirepeat = "trirepeat";
inline constexpr std::string_view kNameQuadrepeat = "quadrepeat";
inline constexpr std::string_view kNameVariantDensity = "variant_density";
inline constexpr std::string_view kNameRefAD = "ref_ad";
inline constexpr std::string_view kNameAltAD = "alt_ad";
inline constexpr std::string_view kNameAltAD2 = "alt_ad2";
inline constexpr std::string_view kNameRefAdAF = "ref_ad_af";
inline constexpr std::string_view kNameAltAdAF = "alt_ad_af";
inline constexpr std::string_view kNameAltAd2AF = "alt_ad2_af";
inline constexpr std::string_view kNamePopAF = "popaf";
inline constexpr std::string_view kNameGenotype = "genotype";
inline constexpr std::string_view kNameQual = "vcf_variant_qual";
inline constexpr std::string_view kNameGq = "vcf_variant_gq";
inline constexpr std::string_view kNameHapcomp = "hapcomp";
inline constexpr std::string_view kNameHapdom = "hapdom";
inline constexpr std::string_view kNameMqRankSum = "mq_rank_sum";
inline constexpr std::string_view kNameSor = "sor";
inline constexpr std::string_view kNameRu = "ru";
inline constexpr std::string_view kNameRpaRef = "rpa_ref";
inline constexpr std::string_view kNameRpaAlt = "rpa_alt";
inline constexpr std::string_view kNameStr = "str";
inline constexpr std::string_view kNameAtInterest = "at_interest";
inline constexpr std::string_view kNameInCatalog = "in_catalog";
inline constexpr std::string_view kNameGcContent = "gc_content";
// tumor-normal specific VCF features
inline constexpr std::string_view kNameTumorAltAD = "tumor_alt_ad";
inline constexpr std::string_view kNameNormalAltAD = "normal_alt_ad";
inline constexpr std::string_view kNameTumorAF = "tumor_af";
inline constexpr std::string_view kNameNormalAF = "normal_af";
inline constexpr std::string_view kNameVcfTnAfRatio = "vcf_tn_af_ratio";
inline constexpr std::string_view kNameTumorDP = "tumor_dp";
inline constexpr std::string_view kNameNormalDP = "normal_dp";

// Panel of Normals features
inline constexpr std::string_view kNamePonOcc = "pon_occ";
inline constexpr std::string_view kNamePonOccRatio = "pon_occ_ratio";
inline constexpr std::string_view kNamePonMeanDuplexAf = "pon_mean_duplex_af";
inline constexpr std::string_view kNamePonMinDuplexAf = "pon_min_duplex_af";
inline constexpr std::string_view kNamePonMaxDuplexAf = "pon_max_duplex_af";

// Computed Germline InDel Fields
inline constexpr std::string_view kNameADT = "adt";
inline constexpr std::string_view kNameADTL = "adtl";
inline constexpr std::string_view kNameNumAlt = "num_alt";
inline constexpr std::string_view kNameIndelAF = "indel_af";

// output fields
inline constexpr std::string_view kNameMLScore = "ml_score";
inline constexpr std::string_view kNameFilterStatus = "filter_status";

inline constexpr auto kVariantIDFeatureNames = std::to_array<std::string_view>(
    {kNameChrom, kNamePos, kNameRef, kNameAlt, kNameAltLen, kNameSubtypeIndex, kNameVariantType});
inline constexpr auto kBamFeatureNames = std::to_array<std::string_view>({kNameADT,
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
                                                                          kNameDuplexConcordant,
                                                                          kNameDuplexSimplex,
                                                                          kNameDuplexDiscordant,
                                                                          kNameSimplex,
                                                                          kNameDuplexDP,
                                                                          kNameDuplexAF,
                                                                          kNameDuplexOther,
                                                                          kNameDuplexOtherAF,
                                                                          kNameDuplexGQ,
                                                                          kNameAfqGT01,
                                                                          kNameAfqGT11,
                                                                          kNameAfqGT12,
                                                                          kNameAfqGT,
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
                                                                          kNameSupportReverse,
                                                                          kNameBamTnAfRatio,
                                                                          kNameAlignmentBias,
                                                                          kNameWeightedDepth,
                                                                          kNameWeightedScore,
                                                                          kNameDel1bpAF,
                                                                          kNameDel2bpAF,
                                                                          kNameDel3bpAF,
                                                                          kNameDelGt3bpAF,
                                                                          kNameIns1bpAF,
                                                                          kNameIns2bpAF,
                                                                          kNameIns3bpAF,
                                                                          kNameInsGt3bpAF,
                                                                          kNameDel1bp2bpAFLogRatio,
                                                                          kNameDel2bp3bpAFLogRatio,
                                                                          kNameDel3bpGt3bpAFLogRatio,
                                                                          kNameIns1bp2bpAFLogRatio,
                                                                          kNameIns2bp3bpAFLogRatio,
                                                                          kNameIns3bpGt3bpAFLogRatio,
                                                                          kNameDel1bpIns1bpAFLogRatio,
                                                                          kNameSnvAF});

inline constexpr auto kVcfFeatureNames = std::to_array<std::string_view>({kNameNalod,
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
                                                                          kNameVcfTnAfRatio,
                                                                          kNameTumorDP,
                                                                          kNameNormalDP,
                                                                          kNamePopAF,
                                                                          kNamePonOcc,
                                                                          kNamePonOccRatio,
                                                                          kNamePonMeanDuplexAf,
                                                                          kNamePonMinDuplexAf,
                                                                          kNamePonMaxDuplexAf,
                                                                          kNameGenotype,
                                                                          kNameQual,
                                                                          kNameGq,
                                                                          kNameHapcomp,
                                                                          kNameHapdom,
                                                                          kNameMqRankSum,
                                                                          kNameSor,
                                                                          kNameRu,
                                                                          kNameRpaRef,
                                                                          kNameRpaAlt,
                                                                          kNameStr,
                                                                          kNameAtInterest,
                                                                          kNameInCatalog,
                                                                          kNameGcContent,
                                                                          kNameADT,
                                                                          kNameADTL,
                                                                          kNameNumAlt,
                                                                          kNameIndelAF,
                                                                          kNameIsInRepeat,
                                                                          kNameRepeatLen,
                                                                          kNameRepeatUnitLen,
                                                                          kNameIsUnitSizeMatch});

/**
 * @brief Find feature names not supported.
 * @param names Vector of feature names to be evaluated
 * @param supported_names Set of feature names supported
 * @param allow_tumor_normal_prefix Allow "tumor_" or "normal_" prefix in feature names
 * @return Vector of feature names not supported
 */
vec<std::string> FindUnsupportedFeatureNames(const vec<std::string>& names,
                                             const StrUnorderedSet& supported_names,
                                             bool allow_tumor_normal_prefix);

/**
 * @brief Find BAM feature names not supported.
 * @param names Vector of feature names to be evaluated
 * @param allow_tumor_normal_prefix Allow "tumor_" or "normal_" prefix in feature names
 * @return Vector of feature names not supported
 */
vec<std::string> FindUnsupportedBamFeatureNames(const vec<std::string>& names, bool allow_tumor_normal_prefix);

/**
 * @brief Find VCF feature names not supported.
 * @param names Vector of feature names to be evaluated
 * @return Vector of feature names not supported
 * @note VCF feature names already have hard-coded sample context prefixes. So, we check for unsupported feature names
 * without inferring sample context.
 */
vec<std::string> FindUnsupportedVcfFeatureNames(const vec<std::string>& names);

/**
 * @brief Find scoring feature names not supported.
 * @param names Vector of feature names to be evaluated
 * @param allow_tumor_normal_prefix Allow "tumor_" or "normal_" prefix in feature names
 * @return Vector of feature names not supported
 */
vec<std::string> FindUnsupportedScoringFeatureNames(const vec<std::string>& names, bool allow_tumor_normal_prefix);

}  // namespace xoos::svc
