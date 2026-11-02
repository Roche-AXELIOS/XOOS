#pragma once

#include <string>

namespace xoos::svc {

/**
 * Synopsis:
 * Consolidate all VCF field names here to simplify future changes and to avoid spelling mistakes.
 */

// Standard VCF field names
static const std::string kFieldAc{"AC"};
static const std::string kFieldAd{"AD"};
static const std::string kFieldAf{"AF"};
static const std::string kFieldAn{"AN"};
static const std::string kFieldDp{"DP"};
static const std::string kFieldGq{"GQ"};
static const std::string kFieldGt{"GT"};

// GATK-specific VCF field names
static const std::string kFieldPopaf{"POPAF"};
static const std::string kFieldNalod{"NALOD"};
static const std::string kFieldNlod{"NLOD"};
static const std::string kFieldTlod{"TLOD"};
static const std::string kFieldMpos{"MPOS"};
static const std::string kFieldMmq{"MMQ"};
static const std::string kFieldMbq{"MBQ"};
static const std::string kFieldHapcomp{"HAPCOMP"};
static const std::string kFieldHapdom{"HAPDOM"};
static const std::string kFieldMqRankSum{"MQRankSum"};
static const std::string kFieldSor{"SOR"};
static const std::string kFieldMleac{"MLEAC"};
static const std::string kFieldMleaf{"MLEAF"};
static const std::string kFieldStr{"STR"};
static const std::string kFieldRu{"RU"};
static const std::string kFieldRpa{"RPA"};

// Panel of Normals VCF field names
static const std::string kFieldPonOcc{"PON_OCC"};
static const std::string kFieldPonOccRatio{"PON_OCC_RATIO"};
static const std::string kFieldPonMeanDuplexAf{"PON_MEAN_DUPLEX_AF"};
static const std::string kFieldPonMinDuplexAf{"PON_MIN_DUPLEX_AF"};
static const std::string kFieldPonMaxDuplexAf{"PON_MAX_DUPLEX_AF"};

// Output VCF field names
static const std::string kWeightedCountsId = "WC";
static const std::string kMachineLearningId = "PRED_ML";
static const std::string kNonDuplexCountsId = "ND";
static const std::string kDuplexCountsId = "DC";
static const std::string kStrandBiasMetricId = "SBM";
static const std::string kPlusOnlyCountsId = "PO";
static const std::string kMinusOnlyCountsId = "MO";
static const std::string kSequenceContextId = "SC";
static const std::string kAlleleDepthId = "AD";
static const std::string kPhysicalPhasedId = "PID";
static const std::string kHotspotId = "HS";
static const std::string kGermlineMLId = "ML_PROCESSED";
static const std::string kDuplexConcordantCountsId = "ADC";
static const std::string kDuplexSimplexCountsId = "ADS";
static const std::string kDuplexDiscordantCountsId = "ADD";
static const std::string kDuplexLowbqCountsId = "ADL";
static const std::string kGermlineGnomadAFId = "GNOMAD_AF";
static const std::string kGermlineRefAvgMapqId = "REF_AVG_MAPQ";
static const std::string kGermlineAltAvgMapqId = "ALT_AVG_MAPQ";
static const std::string kGermlineRefAvgDistId = "REF_AVG_DIST";
static const std::string kGermlineAltAvgDistId = "ALT_AVG_DIST";
static const std::string kGermlineDensity100BPId = "DENSITY_100BP";
static const std::string kBaseqQualId = "BQ";
static const std::string kMapQualId = "MQ";
static const std::string kDistanceId = "DT";
static const std::string kRefBQId = "REFBQ";
static const std::string kRefMQId = "REFMQ";
static const std::string kAltBQId = "ALTBQ";
static const std::string kAltMQId = "ALTMQ";
static const std::string kSubtypeId = "SUBTYPE";
static const std::string kContextId = "CONTEXT";
static const std::string kGermlineTaggingInfoGermlineId = "GERMLINE";
static const std::string kGermlineTaggingInfoSomaticId = "SOMATIC";

// Descriptions for output VCF fields
static const std::string kWeightedCountsDesc = "Weighted counts of molecules supporting alternate allele";
static const std::string kMachineLearningDesc = "Machine learning probability score for the predicted genotype";
static const std::string kNonDuplexCountsDesc = "Counts of non-duplex molecules supporting alternate allele";
static const std::string kDuplexCountsDesc = "Counts of duplex molecule supporting alternate allele";
static const std::string kStrandBiasMetricDesc = "Strand Bias metric";
static const std::string kPlusOnlyCountsDesc = "Counts of alternate allele with only support from + strand";
static const std::string kMinusOnlyCountsDesc = "Counts of alternate alleles with only support from - strand";
static const std::string kSequenceContextDesc = "2 bp sequence context";
static const std::string kPhysicalPhasedDesc = "Physical Phasing ID";
static const std::string kHotspotDesc = "Hotspot Variant";
static const std::string kGermlineMLDesc = "1 if variant was processed by ML and 0 otherwise";
static const std::string kDuplexConcordantCountsDesc =
    "Allelic depths of concordant molecules for the ref and alt alleles in the order listed";
static const std::string kDuplexSimplexCountsDesc =
    "Allelic depths of simplex molecules for the ref and alt alleles in the order listed";
static const std::string kDuplexDiscordantCountsDesc =
    "Allelic depths of discordant molecules for the ref and alt alleles in the order listed";
static const std::string kDuplexLowbqCountsDesc =
    "Allelic depths of lowbq molecules for the ref and alt alleles in the order listed";
static const std::string kGermlineGnomadAFDesc = "Gnomad allele frequency or -1 if not used by ML";
static const std::string kGermlineRefAvgMapqDesc =
    "Average mapping quality for reads supporting REF or -1 if not used in ML";
static const std::string kGermlineAltAvgMapqDesc =
    "Average mapping quality for reads supporting ALT or -1 if not used in ML";
static const std::string kGermlineRefAvgDistDesc =
    "Average distance from variant site to read end for reads supporting REF or -1 if not used in ML";
static const std::string kGermlineAltAvgDistDesc =
    "Average distance from variant site to read end for reads supporting ALT or -1 if not used in ML";
static const std::string kGermlineDensity100BPDesc =
    "Number of other variants called by GATK within 100bp of the call site or -1 if not used in ML";
static const std::string kBaseQualDesc = "Mean base quality of alternate allele";
static const std::string kMapQualDesc = "Mean mapping quality of alternate allele";
static const std::string kDistanceDesc = "Mean distance of alternate allele";
static const std::string kRefBQDesc = "Mean base quality of reference allele";
static const std::string kRefMQDesc = "Mean mapping quality of reference allele";
static const std::string kAltBQDesc = "Mean base quality of alternate allele";
static const std::string kAltMQDesc = "Mean mapping quality of alternate allele";
static const std::string kSubtypeDesc = "Type of substitution";
static const std::string kContextDesc = "Sequence context around the variant";
static const std::string kPopafDesc = "Population allele frequency";
static const std::string kGermlineTaggingInfoGermlineDesc = "Alternate allele is a germline variant";
static const std::string kGermlineTaggingInfoSomaticDesc = "Alternate allele is a somatic variant";

}  // namespace xoos::svc
