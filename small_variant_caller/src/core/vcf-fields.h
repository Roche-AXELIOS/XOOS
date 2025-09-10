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
static const std::string kFieldMleac{"MLEAC"};
static const std::string kFieldMleaf{"MLEAF"};
static const std::string kFieldStr{"STR"};
static const std::string kFieldRu{"RU"};
static const std::string kFieldRpa{"RPA"};

// Output VCF field names
const std::string kWeightedCountsId = "WC";
const std::string kMachineLearningId = "ML";
const std::string kNonDuplexCountsId = "ND";
const std::string kDuplexCountsId = "DC";
const std::string kStrandBiasMetricId = "SBM";
const std::string kPlusOnlyCountsId = "PO";
const std::string kMinusOnlyCountsId = "MO";
const std::string kSequenceContextId = "SC";
const std::string kAlleleDepthId = "AD";
const std::string kPhysicalPhasedId = "PID";
const std::string kHotspotId = "HS";
const std::string kGermlineMLId = "ML_PROCESSED";
const std::string kGermlineGATKDPId = "GATK_DP";
const std::string kGermlineGATKGTId = "GATK_GT";
const std::string kGermlineGATKADId = "GATK_AD";
const std::string kGermlineGATKAltId = "GATK_ALT";
const std::string kGermlineGnomadAFId = "GNOMAD_AF";
const std::string kGermlineRefAvgMapqId = "REF_AVG_MAPQ";
const std::string kGermlineAltAvgMapqId = "ALT_AVG_MAPQ";
const std::string kGermlineRefAvgDistId = "REF_AVG_DIST";
const std::string kGermlineAltAvgDistId = "ALT_AVG_DIST";
const std::string kGermlineDensity100BPId = "DENSITY_100BP";
const std::string kBaseqQualId = "BQ";
const std::string kMapQualId = "MQ";
const std::string kDistanceId = "DT";
const std::string kPredId = "PRED";
const std::string kRefBQId = "REFBQ";
const std::string kRefMQId = "REFMQ";
const std::string kAltBQId = "ALTBQ";
const std::string kAltMQId = "ALTMQ";
const std::string kSubtypeId = "SUBTYPE";
const std::string kContextId = "CONTEXT";

// Descriptions for output VCF fields
const std::string kWeightedCountsDesc = "Weighted counts of molecules supporting alternate allele";
const std::string kMachineLearningDesc = "Machine learning probability score";
const std::string kNonDuplexCountsDesc = "Counts of non-duplex molecules supporting alternate allele";
const std::string kDuplexCountsDesc = "Counts of duplex molecule supporting alternate allele";
const std::string kStrandBiasMetricDesc = "Strand Bias metric";
const std::string kPlusOnlyCountsDesc = "Counts of alternate allele with only support from + strand";
const std::string kMinusOnlyCountsDesc = "Counts of alternate alleles with only support from - strand";
const std::string kSequenceContextDesc = "2 bp sequence context";
const std::string kPhysicalPhasedDesc = "Physical Phasing ID";
const std::string kHotspotDesc = "Hotspot Variant";
const std::string kGermlineMLDesc = "1 if variant was processed by ML and 0 otherwise";
const std::string kGermlineGATKDPDesc = "DP reported by GATK";
const std::string kGermlineGATKGTDesc = "GT reported by GATK";
const std::string kGermlineGATKADDesc = "AD reported by GATK";
const std::string kGermlineGATKAltDesc = "ALT reported by GATK";
const std::string kGermlineGnomadAFDesc = "Gnomad allele frequency or -1 if not used by ML";
const std::string kGermlineRefAvgMapqDesc = "Average mapping quality for reads supporting REF or -1 if not used in ML";
const std::string kGermlineAltAvgMapqDesc = "Average mapping quality for reads supporting ALT or -1 if not used in ML";
const std::string kGermlineRefAvgDistDesc =
    "Average distance from variant site to read end for reads supporting REF or -1 if not used in ML";
const std::string kGermlineAltAvgDistDesc =
    "Average distance from variant site to read end for reads supporting ALT or -1 if not used in ML";
const std::string kGermlineDensity100BPDesc =
    "Number of other variants called by GATK within 100bp of the call site or -1 if not used in ML";
const std::string kBaseQualDesc = "Mean base quality of alternate allele";
const std::string kMapQualDesc = "Mean mapping quality of alternate allele";
const std::string kDistanceDesc = "Mean distance of alternate allele";
const std::string kPredDesc = "Predicted probability score from ML model";
const std::string kRefBQDesc = "Mean base quality of reference allele";
const std::string kRefMQDesc = "Mean mapping quality of reference allele";
const std::string kAltBQDesc = "Mean base quality of alternate allele";
const std::string kAltMQDesc = "Mean mapping quality of alternate allele";
const std::string kSubtypeDesc = "Type of substitution";
const std::string kContextDesc = "Sequence context around the variant";
const std::string kPopafDesc = "Population allele frequency";

}  // namespace xoos::svc
