#include "variant-info.h"

#include <algorithm>
#include <array>

#include "column-names.h"
#include "util/seq-util.h"
#include "xoos/error/error.h"

namespace xoos::svc {

// base context at the variant's start position
static constexpr auto kContexts = {
    "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT", ""};

/**
 * Converts a two basepair context into a numeric index
 * @param context a two basepair string of A,C,G or T
 * @return an int representing the context
 */
int UnifiedVariantFeature::ContextIndex(const std::string& context) {
  const auto* iter = std::find(std::cbegin(kContexts), std::cend(kContexts), context);
  if (iter == std::cend(kContexts)) {
    throw error::Error("Invalid context: " + context);
  }
  return static_cast<int>(std::distance(std::cbegin(kContexts), iter) + 1);
}

static constexpr auto kSubTypes = {"AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "ID", "TA", "TC", "TG"};
static constexpr auto kDefaultSubType = "ID";

static std::string DetermineSubType(const std::string& ref, const std::string& alt) {
  // TODO: add better handling of non-ACTG IUPAC codes, for now just map them to ID
  // TODO: alt == "I" or "D" is for legacy support; to be removed if feature file is no longer required
  if (ref.size() != 1 || IsAnyNotACTG(ref)) {
    return kDefaultSubType;
  }
  if (alt.size() != 1 || IsAnyNotACTG(alt)) {
    return kDefaultSubType;
  }
  return ref + alt;
}

/**
 * @brief Returns the 1-based substitution index for a given reference and alt base
 * @param ref Reference base
 * @param alt Alt base
 * @return Substitution index
 */
int SubstIndex(const std::string& ref, const std::string& alt) {
  auto sub_type = DetermineSubType(ref, alt);
  const auto* iter = std::find(std::cbegin(kSubTypes), std::cend(kSubTypes), sub_type);
  if (iter == std::cend(kSubTypes)) {
    throw error::Error("Invalid substitution type: " + sub_type);
  }
  return static_cast<int>(std::distance(std::cbegin(kSubTypes), iter) + 1);
}

/**
 * Return the UnifiedFeatureCols enum value for a given column name
 * @param feature_name a column/feature name as a string
 * @return a UnifiedFeatureCols enum value that matches the requested string
 */
UnifiedFeatureCols GetCol(const std::string& feature_name) {
  if (!kStringToUnifiedFeatureCols.contains(feature_name)) {
    throw error::Error("Requested Feature: " + feature_name + " doesn't match any known feature");
  }
  return kStringToUnifiedFeatureCols.at(feature_name);
}

/**
 * Return the feature/column string that matches a given UnifiedFeatureCols
 * @param col a UnifiedFeatureCols enum value
 * @return a feature/column name string
 */
std::string GetFeatureName(const UnifiedFeatureCols col) {
  return kUnifiedFeatureColsToString.at(col);
}

/**
 * Builds a map between UnifiedFeatureCols enums to BAM/VCF feature column name strings
 * @return resulting map between UnifiedFeatureCols and string
 */
std::unordered_map<UnifiedFeatureCols, std::string> BuildColToString() {
  using enum UnifiedFeatureCols;

  std::unordered_map<UnifiedFeatureCols, std::string> col_to_string;
  col_to_string[kChrom] = kNameChrom;
  col_to_string[kPos] = kNamePos;
  col_to_string[kRef] = kNameRef;
  col_to_string[kAlt] = kNameAlt;
  col_to_string[kContext] = kNameContext;
  col_to_string[kContextIndex] = kNameContextIndex;
  col_to_string[kSubTypeIndex] = kNameSubtypeIndex;
  col_to_string[kWeightedDepth] = kNameWeightedDepth;
  col_to_string[kSupport] = kNameSupport;
  col_to_string[kMapqMin] = kNameMapqMin;
  col_to_string[kMapqMax] = kNameMapqMax;
  col_to_string[kMapqSum] = kNameMapqSum;
  col_to_string[kMapqSumLowbq] = kNameMapqSumLowbq;
  col_to_string[kMapqSumSimplex] = kNameMapqSumSimplex;
  col_to_string[kMapqMean] = kNameMapqMean;
  col_to_string[kMapqMeanLowbq] = kNameMapqMeanLowbq;
  col_to_string[kMapqMeanSimplex] = kNameMapqMeanSimplex;
  col_to_string[kMapqLT60Count] = kNameMapqLT60Count;
  col_to_string[kMapqLT40Count] = kNameMapqLT40Count;
  col_to_string[kMapqLT30Count] = kNameMapqLT30Count;
  col_to_string[kMapqLT20Count] = kNameMapqLT20Count;
  col_to_string[kBaseqMin] = kNameBaseqMin;
  col_to_string[kBaseqMax] = kNameBaseqMax;
  col_to_string[kBaseqSum] = kNameBaseqSum;
  col_to_string[kBaseqMean] = kNameBaseqMean;
  col_to_string[kBaseqLT20Count] = kNameBaseqLT20Count;
  col_to_string[kBaseqLT20Ratio] = kNameBaseqLT20Ratio;
  col_to_string[kDistanceMin] = kNameDistanceMin;
  col_to_string[kDistanceMax] = kNameDistanceMax;
  col_to_string[kDistanceSum] = kNameDistanceSum;
  col_to_string[kDistanceSumLowbq] = kNameDistanceSumLowbq;
  col_to_string[kDistanceSumSimplex] = kNameDistanceSumSimplex;
  col_to_string[kDistanceMean] = kNameDistanceMean;
  col_to_string[kDistanceMeanLowbq] = kNameDistanceMeanLowbq;
  col_to_string[kDistanceMeanSimplex] = kNameDistanceMeanSimplex;
  col_to_string[kFamilysizeSum] = kNameFamilysizeSum;
  col_to_string[kFamilysizeMean] = kNameFamilysizeMean;
  col_to_string[kFamilysizeLT3Count] = kNameFamilysizeLT3Count;
  col_to_string[kFamilysizeLT5Count] = kNameFamilysizeLT5Count;
  col_to_string[kFamilysizeLT3Ratio] = kNameFamilysizeLT3Ratio;
  col_to_string[kFamilysizeLT5Ratio] = kNameFamilysizeLT5Ratio;
  col_to_string[kNonDuplex] = kNameNonDuplex;
  col_to_string[kDuplex] = kNameDuplex;
  col_to_string[kDuplexLowbq] = kNameDuplexLowbq;
  col_to_string[kSimplex] = kNameSimplex;
  col_to_string[kDuplexDP] = kNameDuplexDP;
  col_to_string[kDuplexAF] = kNameDuplexAF;
  col_to_string[kPlusOnly] = kNamePlusOnly;
  col_to_string[kMinusOnly] = kNameMinusOnly;
  col_to_string[kWeightedScore] = kNameWeightedScore;
  col_to_string[kStrandBias] = kNameStrandBias;
  col_to_string[kMLScore] = kNameMLScore;
  col_to_string[kFilterStatus] = kNameFilterStatus;
  col_to_string[kRefWeightedDepth] = kNameRefWeightedDepth;
  col_to_string[kRefNonhomopolymerWeightedDepth] = kNameRefNonHpWeightedDepth;
  col_to_string[kRefSupport] = kNameRefSupport;
  col_to_string[kTumorRefSupport] = kNameTumorRefSupport;
  col_to_string[kNormalRefSupport] = kNameNormalRefSupport;
  col_to_string[kRefNonhomopolymerSupport] = kNameRefNonHpSupport;
  col_to_string[kRefMapqMin] = kNameRefMapqMin;
  col_to_string[kRefMapqMax] = kNameRefMapqMax;
  col_to_string[kRefMapqSum] = kNameRefMapqSum;
  col_to_string[kRefMapqSumLowbq] = kNameRefMapqSumLowbq;
  col_to_string[kRefMapqSumSimplex] = kNameRefMapqSumSimplex;
  col_to_string[kRefMapqMean] = kNameRefMapqMean;
  col_to_string[kRefMapqMeanLowbq] = kNameRefMapqMeanLowbq;
  col_to_string[kRefMapqMeanSimplex] = kNameRefMapqMeanSimplex;
  col_to_string[kRefMapqLT60Count] = kNameRefMapqLT60Count;
  col_to_string[kRefMapqLT40Count] = kNameRefMapqLT40Count;
  col_to_string[kRefMapqLT30Count] = kNameRefMapqLT30Count;
  col_to_string[kRefMapqLT20Count] = kNameRefMapqLT20Count;
  col_to_string[kRefMapqLT60Ratio] = kNameRefMapqLT60Ratio;
  col_to_string[kRefMapqLT40Ratio] = kNameRefMapqLT40Ratio;
  col_to_string[kRefMapqLT30Ratio] = kNameRefMapqLT30Ratio;
  col_to_string[kRefMapqLT20Ratio] = kNameRefMapqLT20Ratio;
  col_to_string[kRefBaseqSum] = kNameRefBaseqSum;
  col_to_string[kRefBaseqMean] = kNameRefBaseqMean;
  col_to_string[kRefBaseqLT20Count] = kNameRefBaseqLT20Count;
  col_to_string[kRefBaseqLT20Ratio] = kNameRefBaseqLT20Ratio;
  col_to_string[kRefDistanceMin] = kNameRefDistanceMin;
  col_to_string[kRefDistanceMax] = kNameRefDistanceMax;
  col_to_string[kRefDistanceSum] = kNameRefDistanceSum;
  col_to_string[kRefDistanceSumLowbq] = kNameRefDistanceSumLowbq;
  col_to_string[kRefDistanceSumSimplex] = kNameRefDistanceSumSimplex;
  col_to_string[kRefDistanceMean] = kNameRefDistanceMean;
  col_to_string[kRefDistanceMeanLowbq] = kNameRefDistanceMeanLowbq;
  col_to_string[kRefDistanceMeanSimplex] = kNameRefDistanceMeanSimplex;
  col_to_string[kRefDuplexLowbq] = kNameRefDuplexLowbq;
  col_to_string[kRefSimplex] = kNameRefSimplex;
  col_to_string[kRefDuplexAF] = kNameRefDuplexAF;
  col_to_string[kRefFamilysizeSum] = kNameRefFamilysizeSum;
  col_to_string[kRefFamilysizeMean] = kNameRefFamilysizeMean;
  col_to_string[kRefFamilysizeLT3Count] = kNameRefFamilysizeLT3Count;
  col_to_string[kRefFamilysizeLT5Count] = kNameRefFamilysizeLT5Count;
  col_to_string[kRefFamilysizeLT3Ratio] = kNameRefFamilysizeLT3Ratio;
  col_to_string[kRefFamilysizeLT5Ratio] = kNameRefFamilysizeLT5Ratio;
  col_to_string[kRefNonhomopolymerMapqMin] = kNameRefNonHpMapqMin;
  col_to_string[kRefNonhomopolymerMapqMax] = kNameRefNonHpMapqMax;
  col_to_string[kRefNonhomopolymerMapqSum] = kNameRefNonHpMapqSum;
  col_to_string[kRefNonhomopolymerMapqMean] = kNameRefNonHpMapqMean;
  col_to_string[kRefNonhomopolymerMapqLT60Count] = kNameRefNonHpMapqLT60Count;
  col_to_string[kRefNonhomopolymerMapqLT40Count] = kNameRefNonHpMapqLT40Count;
  col_to_string[kRefNonhomopolymerMapqLT30Count] = kNameRefNonHpMapqLT30Count;
  col_to_string[kRefNonhomopolymerMapqLT20Count] = kNameRefNonHpMapqLT20Count;
  col_to_string[kRefNonhomopolymerMapqLT60Ratio] = kNameRefNonHpMapqLT60Ratio;
  col_to_string[kRefNonhomopolymerMapqLT40Ratio] = kNameRefNonHpMapqLT40Ratio;
  col_to_string[kRefNonhomopolymerMapqLT30Ratio] = kNameRefNonHpMapqLT30Ratio;
  col_to_string[kRefNonhomopolymerMapqLT20Ratio] = kNameRefNonHpMapqLT20Ratio;
  col_to_string[kRefNonhomopolymerBaseqMin] = kNameRefNonHpBaseqMin;
  col_to_string[kRefNonhomopolymerBaseqMax] = kNameRefNonHpBaseqMax;
  col_to_string[kRefNonhomopolymerBaseqSum] = kNameRefNonHpBaseqSum;
  col_to_string[kRefNonhomopolymerBaseqMean] = kNameRefNonHpBaseqMean;
  col_to_string[kMapqLT60Ratio] = kNameMapqLT60Ratio;
  col_to_string[kMapqLT40Ratio] = kNameMapqLT40Ratio;
  col_to_string[kMapqLT30Ratio] = kNameMapqLT30Ratio;
  col_to_string[kMapqLT20Ratio] = kNameMapqLT20Ratio;
  col_to_string[kMQAF] = kNameMapqAF;
  col_to_string[kBQAF] = kNameBaseqAF;
  col_to_string[kRefMQAF] = kNameRefMapqAF;
  col_to_string[kRefBQAF] = kNameRefBaseqAF;
  col_to_string[kTumorSupport] = kNameTumorSupport;
  col_to_string[kTumorMapqMean] = kNameTumorMapqMean;
  col_to_string[kTumorBaseqMean] = kNameTumorBaseqMean;
  col_to_string[kTumorDistanceMean] = kNameTumorDistanceMean;
  col_to_string[kNormalSupport] = kNameNormalSupport;
  col_to_string[kNormalMapqMean] = kNameNormalMapqMean;
  col_to_string[kNormalBaseqMean] = kNameNormalBaseqMean;
  col_to_string[kNormalDistanceMean] = kNameNormalDistanceMean;
  col_to_string[kBamTumorAF] = kNameBamTumorAF;
  col_to_string[kBamNormalAF] = kNameBamNormalAF;
  col_to_string[kBamRAT] = kNameBamRAT;
  col_to_string[kSupportReverse] = kNameSupportReverse;
  col_to_string[kTumorSupportReverse] = kNameTumorSupportReverse;
  col_to_string[kNormalSupportReverse] = kNameNormalSupportReverse;
  col_to_string[kAlignmentBias] = kNameAlignmentBias;
  col_to_string[kTumorAlignmentBias] = kNameTumorAlignmentBias;
  col_to_string[kNormalAlignmentBias] = kNameNormalAlignmentBias;

  // VCF features
  col_to_string[kVcfNalod] = kNameNalod;
  col_to_string[kVcfNlod] = kNameNlod;
  col_to_string[kVcfTlod] = kNameTlod;
  col_to_string[kVcfMpos] = kNameMpos;
  col_to_string[kVcfMmqRef] = kNameMmqRef;
  col_to_string[kVcfMmqAlt] = kNameMmqAlt;
  col_to_string[kVcfMbqRef] = kNameMbqRef;
  col_to_string[kVcfMbqAlt] = kNameMbqAlt;
  col_to_string[kVariantType] = kNameVariantType;
  col_to_string[kVcfPre2bContext] = kNamePre2BpContext;
  col_to_string[kVcfPost2bContext] = kNamePost2BpContext;
  col_to_string[kVcfPost30bContext] = kNamePost30BpContext;
  col_to_string[kVcfUnique3mers] = kNameUniq3mers;
  col_to_string[kVcfUnique4mers] = kNameUniq4mers;
  col_to_string[kVcfUnique5mers] = kNameUniq5mers;
  col_to_string[kVcfUnique6mers] = kNameUniq6mers;
  col_to_string[kVcfHomopolymer] = kNameHomopolymer;
  col_to_string[kVcfDirepeat] = kNameDirepeat;
  col_to_string[kVcfTrirepeat] = kNameTrirepeat;
  col_to_string[kVcfQuadrepeat] = kNameQuadrepeat;
  col_to_string[kVcfVariantDensity] = kNameVariantDensity;
  col_to_string[kVcfRefAd] = kNameRefAD;
  col_to_string[kVcfAltAd] = kNameAltAD;
  col_to_string[kVcfAltAd2] = kNameAltAD2;
  col_to_string[kVcfRefAdAf] = kNameRefAdAF;
  col_to_string[kVcfAltAdAf] = kNameAltAdAF;
  col_to_string[kVcfAltAd2Af] = kNameAltAd2AF;
  col_to_string[kVcfTumorAltAd] = kNameTumorAltAD;
  col_to_string[kVcfNormalAltAd] = kNameNormalAltAD;
  col_to_string[kVcfTumorAf] = kNameTumorAF;
  col_to_string[kVcfNormalAf] = kNameNormalAF;
  col_to_string[kVcfTumorNormalAfRatio] = kNameTumorNormalAfRatio;
  col_to_string[kVcfTumorDp] = kNameTumorDP;
  col_to_string[kVcfNormalDp] = kNameNormalDP;
  col_to_string[kVcfPopAf] = kNamePopAF;
  col_to_string[kVcfGenotype] = kNameGenotype;
  col_to_string[kVcfQual] = kNameQual;
  col_to_string[kVcfGq] = kNameGq;
  col_to_string[kVcfHapcomp] = kNameHapcomp;
  col_to_string[kVcfHapdom] = kNameHapdom;
  col_to_string[kVcfRu] = kNameRu;
  col_to_string[kVcfRpaRef] = kNameRpaRef;
  col_to_string[kVcfRpaAlt] = kNameRpaAlt;
  col_to_string[kVcfStr] = kNameStr;
  col_to_string[kVcfAtInterest] = kNameAtInterest;
  col_to_string[kADT] = kNameADT;
  col_to_string[kADTL] = kNameADTL;
  col_to_string[kAltLen] = kNameAltLen;
  col_to_string[kNumAlt] = kNameNumAlt;
  col_to_string[kIndelAf] = kNameIndelAF;

  return col_to_string;
}

/**
 * Builds a map of BAM/VCF feature column name strings to UnifiedFeatureCols enums.
 * @return resulting map between string and UnifiedFeatureCols
 */
StrUnorderedMap<UnifiedFeatureCols> BuildStringToCol() {
  using enum UnifiedFeatureCols;

  StrUnorderedMap<UnifiedFeatureCols> string_to_col;
  string_to_col[kNameChrom] = kChrom;
  string_to_col[kNamePos] = kPos;
  string_to_col[kNameRef] = kRef;
  string_to_col[kNameAlt] = kAlt;
  string_to_col[kNameContext] = kContext;
  string_to_col[kNameContextIndex] = kContextIndex;
  string_to_col[kNameSubtypeIndex] = kSubTypeIndex;
  string_to_col[kNameWeightedDepth] = kWeightedDepth;
  string_to_col[kNameSupport] = kSupport;
  string_to_col[kNameMapqMin] = kMapqMin;
  string_to_col[kNameMapqMax] = kMapqMax;
  string_to_col[kNameMapqSum] = kMapqSum;
  string_to_col[kNameMapqSumLowbq] = kMapqSumLowbq;
  string_to_col[kNameMapqSumSimplex] = kMapqSumSimplex;
  string_to_col[kNameMapqMean] = kMapqMean;
  string_to_col[kNameMapqMeanLowbq] = kMapqMeanLowbq;
  string_to_col[kNameMapqMeanSimplex] = kMapqMeanSimplex;
  string_to_col[kNameMapqLT60Count] = kMapqLT60Count;
  string_to_col[kNameMapqLT40Count] = kMapqLT40Count;
  string_to_col[kNameMapqLT30Count] = kMapqLT30Count;
  string_to_col[kNameMapqLT20Count] = kMapqLT20Count;
  string_to_col[kNameMapqLT60Ratio] = kMapqLT60Ratio;
  string_to_col[kNameMapqLT40Ratio] = kMapqLT40Ratio;
  string_to_col[kNameMapqLT30Ratio] = kMapqLT30Ratio;
  string_to_col[kNameMapqLT20Ratio] = kMapqLT20Ratio;
  string_to_col[kNameBaseqMin] = kBaseqMin;
  string_to_col[kNameBaseqMax] = kBaseqMax;
  string_to_col[kNameBaseqSum] = kBaseqSum;
  string_to_col[kNameBaseqMean] = kBaseqMean;
  string_to_col[kNameBaseqLT20Count] = kBaseqLT20Count;
  string_to_col[kNameBaseqLT20Ratio] = kBaseqLT20Ratio;
  string_to_col[kNameDistanceMin] = kDistanceMin;
  string_to_col[kNameDistanceMax] = kDistanceMax;
  string_to_col[kNameDistanceSum] = kDistanceSum;
  string_to_col[kNameDistanceSumLowbq] = kDistanceSumLowbq;
  string_to_col[kNameDistanceSumSimplex] = kDistanceSumSimplex;
  string_to_col[kNameDistanceMean] = kDistanceMean;
  string_to_col[kNameDistanceMeanLowbq] = kDistanceMeanLowbq;
  string_to_col[kNameDistanceMeanSimplex] = kDistanceMeanSimplex;
  string_to_col[kNameFamilysizeSum] = kFamilysizeSum;
  string_to_col[kNameFamilysizeMean] = kFamilysizeMean;
  string_to_col[kNameFamilysizeLT3Count] = kFamilysizeLT3Count;
  string_to_col[kNameFamilysizeLT5Count] = kFamilysizeLT5Count;
  string_to_col[kNameFamilysizeLT3Ratio] = kFamilysizeLT3Ratio;
  string_to_col[kNameFamilysizeLT5Ratio] = kFamilysizeLT5Ratio;
  string_to_col[kNameNonDuplex] = kNonDuplex;
  string_to_col[kNameDuplex] = kDuplex;
  string_to_col[kNameDuplexLowbq] = kDuplexLowbq;
  string_to_col[kNameSimplex] = kSimplex;
  string_to_col[kNameDuplexDP] = kDuplexDP;
  string_to_col[kNameDuplexAF] = kDuplexAF;
  string_to_col[kNamePlusOnly] = kPlusOnly;
  string_to_col[kNameMinusOnly] = kMinusOnly;
  string_to_col[kNameWeightedScore] = kWeightedScore;
  string_to_col[kNameStrandBias] = kStrandBias;
  string_to_col[kNameMLScore] = kMLScore;
  string_to_col[kNameMapqAF] = kMQAF;
  string_to_col[kNameBaseqAF] = kBQAF;
  string_to_col[kNameRefBaseqAF] = kRefBQAF;
  string_to_col[kNameRefMapqAF] = kRefMQAF;
  string_to_col[kNameTumorSupport] = kTumorSupport;
  string_to_col[kNameTumorMapqMean] = kTumorMapqMean;
  string_to_col[kNameTumorBaseqMean] = kTumorBaseqMean;
  string_to_col[kNameTumorDistanceMean] = kTumorDistanceMean;
  string_to_col[kNameNormalSupport] = kNormalSupport;
  string_to_col[kNameNormalMapqMean] = kNormalMapqMean;
  string_to_col[kNameNormalBaseqMean] = kNormalBaseqMean;
  string_to_col[kNameNormalDistanceMean] = kNormalDistanceMean;
  string_to_col[kNameBamTumorAF] = kBamTumorAF;
  string_to_col[kNameBamNormalAF] = kBamNormalAF;
  string_to_col[kNameBamRAT] = kBamRAT;
  string_to_col[kNameSupportReverse] = kSupportReverse;
  string_to_col[kNameTumorSupportReverse] = kTumorSupportReverse;
  string_to_col[kNameNormalSupportReverse] = kNormalSupportReverse;
  string_to_col[kNameAlignmentBias] = kAlignmentBias;
  string_to_col[kNameTumorAlignmentBias] = kTumorAlignmentBias;
  string_to_col[kNameNormalAlignmentBias] = kNormalAlignmentBias;
  string_to_col[kNameFilterStatus] = kFilterStatus;
  string_to_col[kNameRefWeightedDepth] = kRefWeightedDepth;
  string_to_col[kNameRefNonHpWeightedDepth] = kRefNonhomopolymerWeightedDepth;
  string_to_col[kNameRefSupport] = kRefSupport;
  string_to_col[kNameTumorRefSupport] = kTumorRefSupport;
  string_to_col[kNameNormalRefSupport] = kNormalRefSupport;
  string_to_col[kNameRefNonHpSupport] = kRefNonhomopolymerSupport;
  string_to_col[kNameRefMapqMin] = kRefMapqMin;
  string_to_col[kNameRefMapqMax] = kRefMapqMax;
  string_to_col[kNameRefMapqSum] = kRefMapqSum;
  string_to_col[kNameRefMapqSumLowbq] = kRefMapqSumLowbq;
  string_to_col[kNameRefMapqSumSimplex] = kRefMapqSumSimplex;
  string_to_col[kNameRefMapqMean] = kRefMapqMean;
  string_to_col[kNameRefMapqMeanLowbq] = kRefMapqMeanLowbq;
  string_to_col[kNameRefMapqMeanSimplex] = kRefMapqMeanSimplex;
  string_to_col[kNameRefMapqLT60Count] = kRefMapqLT60Count;
  string_to_col[kNameRefMapqLT40Count] = kRefMapqLT40Count;
  string_to_col[kNameRefMapqLT30Count] = kRefMapqLT30Count;
  string_to_col[kNameRefMapqLT20Count] = kRefMapqLT20Count;
  string_to_col[kNameRefMapqLT60Ratio] = kRefMapqLT60Ratio;
  string_to_col[kNameRefMapqLT40Ratio] = kRefMapqLT40Ratio;
  string_to_col[kNameRefMapqLT30Ratio] = kRefMapqLT30Ratio;
  string_to_col[kNameRefMapqLT20Ratio] = kRefMapqLT20Ratio;
  string_to_col[kNameRefBaseqSum] = kRefBaseqSum;
  string_to_col[kNameRefBaseqMean] = kRefBaseqMean;
  string_to_col[kNameRefBaseqLT20Count] = kRefBaseqLT20Count;
  string_to_col[kNameRefBaseqLT20Ratio] = kRefBaseqLT20Ratio;
  string_to_col[kNameRefDistanceMin] = kRefDistanceMin;
  string_to_col[kNameRefDistanceMax] = kRefDistanceMax;
  string_to_col[kNameRefDistanceSum] = kRefDistanceSum;
  string_to_col[kNameRefDistanceSumLowbq] = kRefDistanceSumLowbq;
  string_to_col[kNameRefDistanceSumSimplex] = kRefDistanceSumSimplex;
  string_to_col[kNameRefDistanceMean] = kRefDistanceMean;
  string_to_col[kNameRefDistanceMeanLowbq] = kRefDistanceMeanLowbq;
  string_to_col[kNameRefDistanceMeanSimplex] = kRefDistanceMeanSimplex;
  string_to_col[kNameRefDuplexLowbq] = kRefDuplexLowbq;
  string_to_col[kNameRefSimplex] = kRefSimplex;
  string_to_col[kNameRefDuplexAF] = kRefDuplexAF;
  string_to_col[kNameRefFamilysizeSum] = kRefFamilysizeSum;
  string_to_col[kNameRefFamilysizeMean] = kRefFamilysizeMean;
  string_to_col[kNameRefFamilysizeLT3Count] = kRefFamilysizeLT3Count;
  string_to_col[kNameRefFamilysizeLT5Count] = kRefFamilysizeLT5Count;
  string_to_col[kNameRefFamilysizeLT3Ratio] = kRefFamilysizeLT3Ratio;
  string_to_col[kNameRefFamilysizeLT5Ratio] = kRefFamilysizeLT5Ratio;
  string_to_col[kNameRefNonHpMapqMin] = kRefNonhomopolymerMapqMin;
  string_to_col[kNameRefNonHpMapqMax] = kRefNonhomopolymerMapqMax;
  string_to_col[kNameRefNonHpMapqSum] = kRefNonhomopolymerMapqSum;
  string_to_col[kNameRefNonHpMapqMean] = kRefNonhomopolymerMapqMean;
  string_to_col[kNameRefNonHpMapqLT60Count] = kRefNonhomopolymerMapqLT60Count;
  string_to_col[kNameRefNonHpMapqLT40Count] = kRefNonhomopolymerMapqLT40Count;
  string_to_col[kNameRefNonHpMapqLT30Count] = kRefNonhomopolymerMapqLT30Count;
  string_to_col[kNameRefNonHpMapqLT20Count] = kRefNonhomopolymerMapqLT20Count;
  string_to_col[kNameRefNonHpMapqLT60Ratio] = kRefNonhomopolymerMapqLT60Ratio;
  string_to_col[kNameRefNonHpMapqLT40Ratio] = kRefNonhomopolymerMapqLT40Ratio;
  string_to_col[kNameRefNonHpMapqLT30Ratio] = kRefNonhomopolymerMapqLT30Ratio;
  string_to_col[kNameRefNonHpMapqLT20Ratio] = kRefNonhomopolymerMapqLT20Ratio;
  string_to_col[kNameRefNonHpBaseqMin] = kRefNonhomopolymerBaseqMin;
  string_to_col[kNameRefNonHpBaseqMax] = kRefNonhomopolymerBaseqMax;
  string_to_col[kNameRefNonHpBaseqSum] = kRefNonhomopolymerBaseqSum;
  string_to_col[kNameRefNonHpBaseqMean] = kRefNonhomopolymerBaseqMean;
  string_to_col[kNameNalod] = kVcfNalod;
  string_to_col[kNameNlod] = kVcfNlod;
  string_to_col[kNameTlod] = kVcfTlod;
  string_to_col[kNameMpos] = kVcfMpos;
  string_to_col[kNameMmqRef] = kVcfMmqRef;
  string_to_col[kNameMmqAlt] = kVcfMmqAlt;
  string_to_col[kNameMbqRef] = kVcfMbqRef;
  string_to_col[kNameMbqAlt] = kVcfMbqAlt;
  string_to_col[kNameVariantType] = kVariantType;
  string_to_col[kNamePre2BpContext] = kVcfPre2bContext;
  string_to_col[kNamePost2BpContext] = kVcfPost2bContext;
  string_to_col[kNamePost30BpContext] = kVcfPost30bContext;
  string_to_col[kNameUniq3mers] = kVcfUnique3mers;
  string_to_col[kNameUniq4mers] = kVcfUnique4mers;
  string_to_col[kNameUniq5mers] = kVcfUnique5mers;
  string_to_col[kNameUniq6mers] = kVcfUnique6mers;
  string_to_col[kNameHomopolymer] = kVcfHomopolymer;
  string_to_col[kNameDirepeat] = kVcfDirepeat;
  string_to_col[kNameTrirepeat] = kVcfTrirepeat;
  string_to_col[kNameQuadrepeat] = kVcfQuadrepeat;
  string_to_col[kNameVariantDensity] = kVcfVariantDensity;
  string_to_col[kNameRefAD] = kVcfRefAd;
  string_to_col[kNameAltAD] = kVcfAltAd;
  string_to_col[kNameAltAD2] = kVcfAltAd2;
  string_to_col[kNameRefAdAF] = kVcfRefAdAf;
  string_to_col[kNameAltAdAF] = kVcfAltAdAf;
  string_to_col[kNameAltAd2AF] = kVcfAltAd2Af;
  string_to_col[kNameTumorAltAD] = kVcfTumorAltAd;
  string_to_col[kNameNormalAltAD] = kVcfNormalAltAd;
  string_to_col[kNameTumorAF] = kVcfTumorAf;
  string_to_col[kNameNormalAF] = kVcfNormalAf;
  string_to_col[kNameTumorNormalAfRatio] = kVcfTumorNormalAfRatio;
  string_to_col[kNameTumorDP] = kVcfTumorDp;
  string_to_col[kNameNormalDP] = kVcfNormalDp;
  string_to_col[kNamePopAF] = kVcfPopAf;
  string_to_col[kNameGenotype] = kVcfGenotype;
  string_to_col[kNameQual] = kVcfQual;
  string_to_col[kNameGq] = kVcfGq;
  string_to_col[kNameHapcomp] = kVcfHapcomp;
  string_to_col[kNameHapdom] = kVcfHapdom;
  string_to_col[kNameRu] = kVcfRu;
  string_to_col[kNameRpaRef] = kVcfRpaRef;
  string_to_col[kNameRpaAlt] = kVcfRpaAlt;
  string_to_col[kNameStr] = kVcfStr;
  string_to_col[kNameAtInterest] = kVcfAtInterest;
  string_to_col[kNameADT] = kADT;
  string_to_col[kNameADTL] = kADTL;
  string_to_col[kNameAltLen] = kAltLen;
  string_to_col[kNameNumAlt] = kNumAlt;
  string_to_col[kNameIndelAF] = kIndelAf;
  return string_to_col;
}

const std::unordered_map<UnifiedFeatureCols, std::string> kUnifiedFeatureColsToString = BuildColToString();
const StrUnorderedMap<UnifiedFeatureCols> kStringToUnifiedFeatureCols = BuildStringToCol();

bool IsVcfFeatureCol(const UnifiedFeatureCols col) {
  return kVcfFeatureCols.contains(col);
}

}  // namespace xoos::svc
