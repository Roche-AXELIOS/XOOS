#include "variant-info.h"

#include <algorithm>
#include <array>
#include <string_view>

#include <fmt/format.h>

#include "column-names.h"
#include "util/seq-util.h"
#include "xoos/error/error.h"

namespace xoos::svc {

// base context at the variant's start position
static constexpr std::array<std::string_view, 17> kContexts = {
    "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT", ""};

/**
 * Converts a two basepair context into a numeric index
 * @param context a two basepair string of A,C,G or T
 * @return an int representing the context
 */
u32 UnifiedVariantFeature::ContextIndex(const std::string& context) {
  const auto* const iter = std::ranges::find(kContexts, context);
  if (iter == kContexts.end()) {
    throw error::Error("Invalid context: " + context);
  }
  return static_cast<u32>(std::distance(kContexts.begin(), iter) + 1);
}

static constexpr std::array<std::string_view, 13> kSubTypes = {
    "AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "ID", "TA", "TC", "TG"};
static constexpr std::string_view kDefaultSubType = "ID";

namespace {
std::string DetermineSubType(const std::string& ref, const std::string& alt) {
  if (ref.size() != 1 || IsAnyNotACTG(ref)) {
    return std::string{kDefaultSubType};
  }
  if (alt.size() != 1 || IsAnyNotACTG(alt)) {
    return std::string{kDefaultSubType};
  }
  return ref + alt;
}
}  // namespace

/**
 * @brief Returns the 1-based substitution index for a given reference and alt base
 * @param ref Reference base
 * @param alt Alt base
 * @return Substitution index
 */
s32 SubstIndex(const std::string& ref, const std::string& alt) {
  const auto sub_type = DetermineSubType(ref, alt);
  const auto* const iter = std::ranges::find(kSubTypes, sub_type);
  if (iter == kSubTypes.end()) {
    throw error::Error("Invalid substitution type: " + sub_type);
  }
  return static_cast<s32>(std::distance(kSubTypes.begin(), iter) + 1);
}

/**
 * Return the UnifiedFeatureCols enum value for a given column name
 * @param feature_name a column/feature name as a string
 * @return a UnifiedFeatureCols enum value that matches the requested string
 */
UnifiedFeatureCols GetCol(const std::string_view feature_name) {
  const auto itr = kStringToUnifiedFeatureCols.find(feature_name);
  if (itr == kStringToUnifiedFeatureCols.end()) {
    throw error::Error(fmt::format("Requested Feature: {} doesn't match any known feature", feature_name));
  }
  return itr->second;
}

FeatureColumn GetFeatureColumn(const std::string& feature_name, const bool infer_sample_context) {
  using enum SampleContext;
  FeatureColumn col{};
  // Check for sample context prefixes
  const bool has_tumor_prefix = infer_sample_context && feature_name.starts_with(kTumorPrefix);
  const bool has_normal_prefix = infer_sample_context && !has_tumor_prefix && feature_name.starts_with(kNormalPrefix);
  if (has_tumor_prefix) {
    col.sample_context = kTumor;
  } else if (has_normal_prefix) {
    col.sample_context = kNormal;
  } else {
    col.sample_context = kNone;
  }
  // Without a sample-context prefix, resolution is a plain name lookup (throwing on miss), which is exactly GetCol.
  if (!has_tumor_prefix && !has_normal_prefix) {
    col.enum_val = GetCol(feature_name);
    return col;
  }
  // Try to find the feature name in the map, first with the original name, then with prefixes stripped if applicable
  auto itr = kStringToUnifiedFeatureCols.find(feature_name);
  if (itr != kStringToUnifiedFeatureCols.end()) {
    col.enum_val = itr->second;
    return col;
  }
  if (has_tumor_prefix) {
    itr = kStringToUnifiedFeatureCols.find(feature_name.substr(kTumorPrefix.size()));
    if (itr != kStringToUnifiedFeatureCols.end()) {
      col.enum_val = itr->second;
      return col;
    }
  }
  if (has_normal_prefix) {
    itr = kStringToUnifiedFeatureCols.find(feature_name.substr(kNormalPrefix.size()));
    if (itr != kStringToUnifiedFeatureCols.end()) {
      col.enum_val = itr->second;
      return col;
    }
  }
  throw error::Error("Requested Feature: " + feature_name + " doesn't match any known feature");
}

std::string GetFeatureName(const UnifiedFeatureCols col) {
  return kUnifiedFeatureColsToString.at(col);
}

std::string GetFeatureName(const FeatureColumn& col) {
  using enum SampleContext;
  std::string feature_name = kUnifiedFeatureColsToString.at(col.enum_val);
  if (IsVcfFeatureCol(col.enum_val)) {
    // Do not add sample prefixes to VCF feature name, even if the column has sample context.
    return feature_name;
  }
  // For non-VCF features, add sample context prefixes to the feature name if applicable.
  switch (col.sample_context) {
    case kTumor:
      return fmt::format("{}{}", kTumorPrefix, feature_name);
    case kNormal:
      return fmt::format("{}{}", kNormalPrefix, feature_name);
    case kNone:
    default:
      return feature_name;
  }
}

namespace {
/**
 * Builds a map between UnifiedFeatureCols enums to BAM/VCF feature column name strings
 * @return resulting map between UnifiedFeatureCols and string
 */
std::unordered_map<UnifiedFeatureCols, std::string> BuildColToString() noexcept {
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
  col_to_string[kDuplexConcordant] = kNameDuplexConcordant;
  col_to_string[kDuplexSimplex] = kNameDuplexSimplex;
  col_to_string[kDuplexDiscordant] = kNameDuplexDiscordant;
  col_to_string[kSimplex] = kNameSimplex;
  col_to_string[kDuplexDP] = kNameDuplexDP;
  col_to_string[kDuplexAF] = kNameDuplexAF;
  col_to_string[kDuplexOther] = kNameDuplexOther;
  col_to_string[kDuplexOtherAF] = kNameDuplexOtherAF;
  col_to_string[kDuplexGQ] = kNameDuplexGQ;
  col_to_string[kAfqGT01] = kNameAfqGT01;
  col_to_string[kAfqGT11] = kNameAfqGT11;
  col_to_string[kAfqGT12] = kNameAfqGT12;
  col_to_string[kAfqGT] = kNameAfqGT;
  col_to_string[kPlusOnly] = kNamePlusOnly;
  col_to_string[kMinusOnly] = kNameMinusOnly;
  col_to_string[kWeightedScore] = kNameWeightedScore;
  col_to_string[kStrandBias] = kNameStrandBias;
  col_to_string[kMLScore] = kNameMLScore;
  col_to_string[kFilterStatus] = kNameFilterStatus;
  col_to_string[kRefWeightedDepth] = kNameRefWeightedDepth;
  col_to_string[kRefNonhomopolymerWeightedDepth] = kNameRefNonHpWeightedDepth;
  col_to_string[kRefSupport] = kNameRefSupport;
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
  col_to_string[kDel1bpAF] = kNameDel1bpAF;
  col_to_string[kDel2bpAF] = kNameDel2bpAF;
  col_to_string[kDel3bpAF] = kNameDel3bpAF;
  col_to_string[kDelGt3bpAF] = kNameDelGt3bpAF;
  col_to_string[kIns1bpAF] = kNameIns1bpAF;
  col_to_string[kIns2bpAF] = kNameIns2bpAF;
  col_to_string[kIns3bpAF] = kNameIns3bpAF;
  col_to_string[kInsGt3bpAF] = kNameInsGt3bpAF;
  col_to_string[kDel1bp2bpAFLogRatio] = kNameDel1bp2bpAFLogRatio;
  col_to_string[kDel2bp3bpAFLogRatio] = kNameDel2bp3bpAFLogRatio;
  col_to_string[kDel3bpGt3bpAFLogRatio] = kNameDel3bpGt3bpAFLogRatio;
  col_to_string[kIns1bp2bpAFLogRatio] = kNameIns1bp2bpAFLogRatio;
  col_to_string[kIns2bp3bpAFLogRatio] = kNameIns2bp3bpAFLogRatio;
  col_to_string[kIns3bpGt3bpAFLogRatio] = kNameIns3bpGt3bpAFLogRatio;
  col_to_string[kDel1bpIns1bpAFLogRatio] = kNameDel1bpIns1bpAFLogRatio;
  col_to_string[kSnvAF] = kNameSnvAF;
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
  col_to_string[kRefDuplexConcordant] = kNameRefDuplexConcordant;
  col_to_string[kRefDuplexSimplex] = kNameRefDuplexSimplex;
  col_to_string[kRefDuplexDiscordant] = kNameRefDuplexDiscordant;
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
  col_to_string[kBamTnAfRatio] = kNameBamTnAfRatio;
  col_to_string[kSupportReverse] = kNameSupportReverse;
  col_to_string[kAlignmentBias] = kNameAlignmentBias;

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
  col_to_string[kVcfTnAfRatio] = kNameVcfTnAfRatio;
  col_to_string[kVcfTumorDp] = kNameTumorDP;
  col_to_string[kVcfNormalDp] = kNameNormalDP;
  col_to_string[kVcfPopAf] = kNamePopAF;
  col_to_string[kVcfPonOcc] = kNamePonOcc;
  col_to_string[kVcfPonOccRatio] = kNamePonOccRatio;
  col_to_string[kVcfPonMeanDuplexAf] = kNamePonMeanDuplexAf;
  col_to_string[kVcfPonMinDuplexAf] = kNamePonMinDuplexAf;
  col_to_string[kVcfPonMaxDuplexAf] = kNamePonMaxDuplexAf;
  col_to_string[kVcfGenotype] = kNameGenotype;
  col_to_string[kVcfQual] = kNameQual;
  col_to_string[kVcfGq] = kNameGq;
  col_to_string[kVcfHapcomp] = kNameHapcomp;
  col_to_string[kVcfHapdom] = kNameHapdom;
  col_to_string[kVcfMqRankSum] = kNameMqRankSum;
  col_to_string[kVcfSor] = kNameSor;
  col_to_string[kVcfRu] = kNameRu;
  col_to_string[kVcfRpaRef] = kNameRpaRef;
  col_to_string[kVcfRpaAlt] = kNameRpaAlt;
  col_to_string[kVcfStr] = kNameStr;
  col_to_string[kVcfAtInterest] = kNameAtInterest;
  col_to_string[kVcfIsInRepeat] = kNameIsInRepeat;
  col_to_string[kVcfRepeatLen] = kNameRepeatLen;
  col_to_string[kVcfRepeatUnitLen] = kNameRepeatUnitLen;
  col_to_string[kVcfIsUnitSizeMatch] = kNameIsUnitSizeMatch;
  col_to_string[kVcfInCatalog] = kNameInCatalog;
  col_to_string[kVcfGcContent] = kNameGcContent;
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
StrUnorderedMap<UnifiedFeatureCols> BuildStringToCol() noexcept {
  // Inverse of kUnifiedFeatureColsToString. The enum->name map is a bijection, so inverting it keeps both
  // directions in sync from a single source of truth.
  StrUnorderedMap<UnifiedFeatureCols> string_to_col;
  string_to_col.reserve(kUnifiedFeatureColsToString.size());
  for (const auto& [col, name] : kUnifiedFeatureColsToString) {
    string_to_col[name] = col;
  }
  return string_to_col;
}
}  // namespace

const std::unordered_map<UnifiedFeatureCols, std::string> kUnifiedFeatureColsToString = BuildColToString();
const StrUnorderedMap<UnifiedFeatureCols> kStringToUnifiedFeatureCols = BuildStringToCol();

bool IsVcfFeatureCol(const UnifiedFeatureCols col) {
  return kVcfFeatureCols.contains(col);
}

bool IsVcfFeatureColumn(const FeatureColumn& col) {
  return kVcfFeatureCols.contains(col.enum_val);
}

}  // namespace xoos::svc
