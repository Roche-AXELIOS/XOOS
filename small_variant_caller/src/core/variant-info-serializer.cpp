#include "variant-info-serializer.h"

#include <cmath>
#include <fstream>
#include <vector>

#include <xoos/error/error.h>
#include <xoos/util/hash.h>
#include <xoos/util/string-functions.h>

#include "core/column-names.h"
#include "util/seq-util.h"

namespace xoos::svc {

/**
 * @brief Convert string to a dot if it is empty, otherwise return the string as is.
 * @param value Input string value
 * @return Converted string, empty strings become ".", otherwise the original string
 */
static std::string EmptyStringToDot(const std::string& value) {
  return value.empty() ? "." : value;
}

/**
 * @brief Convert a dot to an empty string if the input is ".", otherwise return the string as is.
 * @param value Input string value
 * @return Converted string, "." becomes empty string, otherwise the original string
 */
static std::string DotToEmptyString(const std::string& value) {
  return value == "." ? "" : value;
}

/**
 * @brief Serialize a VariantId struct into a map of feature enum to string value
 * @param vid a VariantId struct
 * @return a map of UnifiedFeatureCols to std::strings
 */
std::map<UnifiedFeatureCols, std::string> VariantInfoSerializer::SerializeVariantId(const VariantId& vid) {
  using enum UnifiedFeatureCols;
  std::map<UnifiedFeatureCols, std::string> result;
  result[kChrom] = vid.chrom;
  result[kPos] = std::to_string(vid.pos + 1);
  result[kRef] = vid.ref;
  result[kAlt] = vid.alt;
  result[kSubTypeIndex] = std::to_string(SubstIndex(vid.ref, vid.alt));
  result[kVariantType] = std::to_string(static_cast<u8>(vid.type));
  result[kAltLen] = std::to_string(static_cast<s64>(vid.alt.size()) - static_cast<s64>(vid.ref.size()));
  return result;
}

/**
 * @brief Serialize a UnifiedVariantFeature struct into a map of feature enum to string value
 * @param variant_info a UnifiedVariantFeature struct
 * @return a map of UnifiedFeatureCols to std::strings
 */
std::map<UnifiedFeatureCols, std::string> VariantInfoSerializer::SerializeVariantFeature(
    const UnifiedVariantFeature& variant_info) {
  using enum UnifiedFeatureCols;
  std::map<UnifiedFeatureCols, std::string> result;
  result[kWeightedDepth] = std::to_string(variant_info.weighted_depth);
  result[kSupport] = std::to_string(variant_info.support);
  result[kMapqMin] = std::to_string(variant_info.mapq_min);
  result[kMapqMax] = std::to_string(variant_info.mapq_max);
  result[kMapqSum] = std::to_string(variant_info.mapq_sum);
  result[kMapqMean] = std::to_string(variant_info.mapq_mean);
  result[kMapqSumLowbq] = std::to_string(variant_info.mapq_sum_lowbq);
  result[kMapqMeanLowbq] = std::to_string(variant_info.mapq_mean_lowbq);
  result[kMapqSumSimplex] = std::to_string(variant_info.mapq_sum_simplex);
  result[kMapqMeanSimplex] = std::to_string(variant_info.mapq_mean_simplex);
  result[kMapqLT60Count] = std::to_string(variant_info.mapq_lt60_count);
  result[kMapqLT40Count] = std::to_string(variant_info.mapq_lt40_count);
  result[kMapqLT30Count] = std::to_string(variant_info.mapq_lt30_count);
  result[kMapqLT20Count] = std::to_string(variant_info.mapq_lt20_count);
  result[kMapqLT60Ratio] = std::to_string(variant_info.mapq_lt60_ratio);
  result[kMapqLT40Ratio] = std::to_string(variant_info.mapq_lt40_ratio);
  result[kMapqLT30Ratio] = std::to_string(variant_info.mapq_lt30_ratio);
  result[kMapqLT20Ratio] = std::to_string(variant_info.mapq_lt20_ratio);
  result[kBaseqMin] = std::to_string(variant_info.baseq_min);
  result[kBaseqMax] = std::to_string(variant_info.baseq_max);
  result[kBaseqSum] = std::to_string(variant_info.baseq_sum);
  result[kBaseqMean] = std::to_string(variant_info.baseq_mean);
  result[kBaseqLT20Count] = std::to_string(variant_info.baseq_lt20_count);
  result[kBaseqLT20Ratio] = std::to_string(variant_info.baseq_lt20_ratio);
  result[kDistanceMin] = std::to_string(variant_info.distance_min);
  result[kDistanceMax] = std::to_string(variant_info.distance_max);
  result[kDistanceSum] = std::to_string(variant_info.distance_sum);
  result[kDistanceMean] = std::to_string(variant_info.distance_mean);
  result[kDistanceSumLowbq] = std::to_string(variant_info.distance_sum_lowbq);
  result[kDistanceMeanLowbq] = std::to_string(variant_info.distance_mean_lowbq);
  result[kDistanceSumSimplex] = std::to_string(variant_info.distance_sum_simplex);
  result[kDistanceMeanSimplex] = std::to_string(variant_info.distance_mean_simplex);
  result[kFamilysizeSum] = std::to_string(variant_info.familysize_sum);
  result[kFamilysizeMean] = std::to_string(variant_info.familysize_mean);
  result[kFamilysizeLT3Count] = std::to_string(variant_info.familysize_lt3_count);
  result[kFamilysizeLT5Count] = std::to_string(variant_info.familysize_lt5_count);
  result[kFamilysizeLT3Ratio] = std::to_string(variant_info.familysize_lt3_ratio);
  result[kFamilysizeLT5Ratio] = std::to_string(variant_info.familysize_lt5_ratio);
  result[kNonDuplex] = std::to_string(variant_info.nonduplex);
  result[kDuplex] = std::to_string(variant_info.duplex);
  result[kDuplexLowbq] = std::to_string(variant_info.duplex_lowbq);
  result[kSimplex] = std::to_string(variant_info.simplex);
  result[kDuplexAF] = std::to_string(variant_info.duplex_af);
  result[kPlusOnly] = std::to_string(variant_info.plusonly);
  result[kMinusOnly] = std::to_string(variant_info.minusonly);
  result[kMQAF] = std::to_string(variant_info.mq_af);
  result[kBQAF] = std::to_string(variant_info.bq_af);
  result[kTumorSupport] = std::to_string(variant_info.tumor_support);
  result[kTumorMapqSum] = std::to_string(variant_info.tumor_mapq_sum);
  result[kTumorMapqMean] = std::to_string(variant_info.tumor_mapq_mean);
  result[kTumorBaseqSum] = std::to_string(variant_info.tumor_baseq_sum);
  result[kTumorBaseqMean] = std::to_string(variant_info.tumor_baseq_mean);
  result[kTumorDistanceSum] = std::to_string(variant_info.tumor_distance_sum);
  result[kTumorDistanceMean] = std::to_string(variant_info.tumor_distance_mean);
  result[kNormalSupport] = std::to_string(variant_info.normal_support);
  result[kNormalMapqSum] = std::to_string(variant_info.normal_mapq_sum);
  result[kNormalMapqMean] = std::to_string(variant_info.normal_mapq_mean);
  result[kNormalBaseqSum] = std::to_string(variant_info.normal_baseq_sum);
  result[kNormalBaseqMean] = std::to_string(variant_info.normal_baseq_mean);
  result[kNormalDistanceSum] = std::to_string(variant_info.normal_distance_sum);
  result[kNormalDistanceMean] = std::to_string(variant_info.normal_distance_mean);
  result[kBamTumorAF] = std::to_string(variant_info.tumor_af);
  result[kBamNormalAF] = std::to_string(variant_info.normal_af);
  result[kBamRAT] = std::to_string(variant_info.rat);
  result[kSupportReverse] = std::to_string(variant_info.support_reverse);
  result[kTumorSupportReverse] = std::to_string(variant_info.tumor_support_reverse);
  result[kNormalSupportReverse] = std::to_string(variant_info.normal_support_reverse);
  result[kAlignmentBias] = std::to_string(variant_info.alignmentbias);
  result[kTumorAlignmentBias] = std::to_string(variant_info.tumor_alignmentbias);
  result[kNormalAlignmentBias] = std::to_string(variant_info.normal_alignmentbias);
  result[kWeightedScore] = std::to_string(variant_info.weighted_score);
  result[kStrandBias] = std::to_string(variant_info.strandbias);
  result[kContext] = EmptyStringToDot(variant_info.context);
  result[kContextIndex] = std::to_string(variant_info.context_index);
  result[kMLScore] = std::to_string(variant_info.ml_score);
  result[kFilterStatus] = string::Join(variant_info.filter_status, ",");
  result[kADT] = std::to_string(variant_info.adt);
  result[kADTL] = std::to_string(variant_info.adtl);
  result[kIndelAf] = std::to_string(variant_info.indel_af);
  return result;
}

/**
 * @brief Serialize a UnifiedReferenceFeature struct into a map of feature enum to string value
 * @param ref_info a UnifiedReferenceFeature struct
 * @return a map of UnifiedFeatureCols to std::strings
 */
std::map<UnifiedFeatureCols, std::string> VariantInfoSerializer::SerializeReferenceFeature(
    const UnifiedReferenceFeature& ref_info) {
  using enum UnifiedFeatureCols;
  std::map<UnifiedFeatureCols, std::string> result;
  result[kRefSupport] = std::to_string(ref_info.support);
  result[kRefWeightedDepth] = std::to_string(ref_info.weighted_depth);
  result[kRefNonhomopolymerWeightedDepth] = std::to_string(ref_info.nonhomopolymer_weighted_depth);
  result[kRefNonhomopolymerSupport] = std::to_string(ref_info.nonhomopolymer_support);
  result[kRefMapqMin] = std::to_string(ref_info.mapq_min);
  result[kRefMapqMax] = std::to_string(ref_info.mapq_max);
  result[kRefMapqSum] = std::to_string(ref_info.mapq_sum);
  result[kRefMapqMean] = std::to_string(ref_info.mapq_mean);
  result[kRefMapqSumLowbq] = std::to_string(ref_info.mapq_sum_lowbq);
  result[kRefMapqMeanLowbq] = std::to_string(ref_info.mapq_mean_lowbq);
  result[kRefMapqSumSimplex] = std::to_string(ref_info.mapq_sum_simplex);
  result[kRefMapqMeanSimplex] = std::to_string(ref_info.mapq_mean_simplex);
  result[kRefMapqLT60Count] = std::to_string(ref_info.mapq_lt60_count);
  result[kRefMapqLT40Count] = std::to_string(ref_info.mapq_lt40_count);
  result[kRefMapqLT30Count] = std::to_string(ref_info.mapq_lt30_count);
  result[kRefMapqLT20Count] = std::to_string(ref_info.mapq_lt20_count);
  result[kRefBaseqSum] = std::to_string(ref_info.baseq_sum);
  result[kRefBaseqMean] = std::to_string(ref_info.baseq_mean);
  result[kRefBaseqLT20Count] = std::to_string(ref_info.baseq_lt20_count);
  result[kRefBaseqLT20Ratio] = std::to_string(ref_info.baseq_lt20_ratio);
  result[kRefDistanceMin] = std::to_string(ref_info.distance_min);
  result[kRefDistanceMax] = std::to_string(ref_info.distance_max);
  result[kRefDistanceSum] = std::to_string(ref_info.distance_sum);
  result[kRefDistanceMean] = std::to_string(ref_info.distance_mean);
  result[kRefDistanceSumLowbq] = std::to_string(ref_info.distance_sum_lowbq);
  result[kRefDistanceMeanLowbq] = std::to_string(ref_info.distance_mean_lowbq);
  result[kRefDistanceSumSimplex] = std::to_string(ref_info.distance_sum_simplex);
  result[kRefDistanceMeanSimplex] = std::to_string(ref_info.distance_mean_simplex);
  result[kRefDuplexLowbq] = std::to_string(ref_info.duplex_lowbq);
  result[kRefSimplex] = std::to_string(ref_info.simplex);
  result[kRefDuplexAF] = std::to_string(ref_info.duplex_af);
  result[kDuplexDP] = std::to_string(ref_info.duplex_dp);
  result[kRefFamilysizeSum] = std::to_string(ref_info.familysize_sum);
  result[kRefFamilysizeMean] = std::to_string(ref_info.familysize_mean);
  result[kRefFamilysizeLT3Count] = std::to_string(ref_info.familysize_lt3_count);
  result[kRefFamilysizeLT5Count] = std::to_string(ref_info.familysize_lt5_count);
  result[kRefFamilysizeLT3Ratio] = std::to_string(ref_info.familysize_lt3_ratio);
  result[kRefFamilysizeLT5Ratio] = std::to_string(ref_info.familysize_lt5_ratio);
  result[kRefNonhomopolymerMapqMin] = std::to_string(ref_info.nonhomopolymer_mapq_min);
  result[kRefNonhomopolymerMapqMax] = std::to_string(ref_info.nonhomopolymer_mapq_max);
  result[kRefNonhomopolymerMapqSum] = std::to_string(ref_info.nonhomopolymer_mapq_sum);
  result[kRefNonhomopolymerMapqMean] = std::to_string(ref_info.nonhomopolymer_mapq_mean);
  result[kRefNonhomopolymerMapqLT60Count] = std::to_string(ref_info.nonhomopolymer_mapq_lt60_count);
  result[kRefNonhomopolymerMapqLT40Count] = std::to_string(ref_info.nonhomopolymer_mapq_lt40_count);
  result[kRefNonhomopolymerMapqLT30Count] = std::to_string(ref_info.nonhomopolymer_mapq_lt30_count);
  result[kRefNonhomopolymerMapqLT20Count] = std::to_string(ref_info.nonhomopolymer_mapq_lt20_count);
  result[kRefNonhomopolymerMapqLT60Ratio] = std::to_string(ref_info.nonhomopolymer_mapq_lt60_ratio);
  result[kRefNonhomopolymerMapqLT40Ratio] = std::to_string(ref_info.nonhomopolymer_mapq_lt40_ratio);
  result[kRefNonhomopolymerMapqLT30Ratio] = std::to_string(ref_info.nonhomopolymer_mapq_lt30_ratio);
  result[kRefNonhomopolymerMapqLT20Ratio] = std::to_string(ref_info.nonhomopolymer_mapq_lt20_ratio);
  result[kRefNonhomopolymerBaseqMin] = std::to_string(ref_info.nonhomopolymer_baseq_min);
  result[kRefNonhomopolymerBaseqMax] = std::to_string(ref_info.nonhomopolymer_baseq_max);
  result[kRefNonhomopolymerBaseqSum] = std::to_string(ref_info.nonhomopolymer_baseq_sum);
  result[kRefNonhomopolymerBaseqMean] = std::to_string(ref_info.nonhomopolymer_baseq_mean);
  result[kRefMapqLT60Ratio] = std::to_string(ref_info.mapq_lt60_ratio);
  result[kRefMapqLT40Ratio] = std::to_string(ref_info.mapq_lt40_ratio);
  result[kRefMapqLT30Ratio] = std::to_string(ref_info.mapq_lt30_ratio);
  result[kRefMapqLT20Ratio] = std::to_string(ref_info.mapq_lt20_ratio);
  result[kTumorRefSupport] = std::to_string(ref_info.tumor_support);
  result[kNormalRefSupport] = std::to_string(ref_info.normal_support);
  result[kRefMQAF] = std::to_string(ref_info.mq_af);
  result[kRefBQAF] = std::to_string(ref_info.bq_af);
  result[kNumAlt] = std::to_string(ref_info.num_alt);
  return result;
}

/**
 * @brief Deserialize a map of UnifiedFeatureCols to string into a VariantId struct
 * @param serialized a map of UnifiedFeatureCols enums to strings
 * @return a VariantId Struct
 */
VariantId VariantInfoSerializer::DeserializeVariantId(const std::map<UnifiedFeatureCols, std::string>& serialized) {
  using enum UnifiedFeatureCols;
  auto get_value = [&](const UnifiedFeatureCols col) {
    const auto itr = serialized.find(col);
    return itr != serialized.end() ? itr->second : "";
  };
  // convert position from 1-based (string) to 0-based (unsigned 64-bit integer)
  const u64 pos = get_value(kPos).empty() ? 0 : std::stol(get_value(kPos)) - 1;
  const VariantId vid(get_value(kChrom), pos, get_value(kRef), get_value(kAlt));
  return vid;
}

/**
 * @brief Deserialize a map of UnifiedFeatureCols to string into a UnifiedVariantFeature struct
 * @param serialized a map of UnifiedFeatureCols enums to strings
 * @return a UnifiedVariantFeature Struct
 */
UnifiedVariantFeature VariantInfoSerializer::DeserializeVariantFeature(
    const std::map<UnifiedFeatureCols, std::string>& serialized) {
  using enum UnifiedFeatureCols;
  UnifiedVariantFeature feat{};
  for (const auto& [col, value] : serialized) {
    switch (col) {
      case kWeightedDepth:
        feat.weighted_depth = std::stod(value);
        break;
      case kSupport:
        feat.support = std::stoul(value);
        break;
      case kMapqMax:
        feat.mapq_max = static_cast<u8>(std::stoi(value));
        break;
      case kMapqMin:
        feat.mapq_min = static_cast<u8>(std::stoi(value));
        break;
      case kMapqSum:
        feat.mapq_sum = std::stoul(value);
        break;
      case kMapqSumLowbq:
        feat.mapq_sum_lowbq = std::stoul(value);
        break;
      case kMapqSumSimplex:
        feat.mapq_sum_simplex = std::stoul(value);
        break;
      case kMapqMean:
        feat.mapq_mean = std::stod(value);
        break;
      case kMapqMeanLowbq:
        feat.mapq_mean_lowbq = std::stod(value);
        break;
      case kMapqMeanSimplex:
        feat.mapq_mean_simplex = std::stod(value);
        break;
      case kMapqLT60Count:
        feat.mapq_lt60_count = std::stoul(value);
        break;
      case kMapqLT40Count:
        feat.mapq_lt40_count = std::stoul(value);
        break;
      case kMapqLT30Count:
        feat.mapq_lt30_count = std::stoul(value);
        break;
      case kMapqLT20Count:
        feat.mapq_lt20_count = std::stoul(value);
        break;
      case kBaseqMin:
        feat.baseq_min = std::stod(value);
        break;
      case kBaseqMax:
        feat.baseq_max = std::stod(value);
        break;
      case kBaseqSum:
        feat.baseq_sum = std::stod(value);
        break;
      case kBaseqMean:
        feat.baseq_mean = std::stod(value);
        break;
      case kBaseqLT20Count:
        feat.baseq_lt20_count = std::stoul(value);
        break;
      case kDistanceMin:
        feat.distance_min = std::stoul(value);
        break;
      case kDistanceMax:
        feat.distance_max = std::stoul(value);
        break;
      case kDistanceSum:
        feat.distance_sum = std::stoul(value);
        break;
      case kDistanceSumLowbq:
        feat.distance_sum_lowbq = std::stoul(value);
        break;
      case kDistanceSumSimplex:
        feat.distance_sum_simplex = std::stoul(value);
        break;
      case kDistanceMean:
        feat.distance_mean = std::stod(value);
        break;
      case kDistanceMeanLowbq:
        feat.distance_mean_lowbq = std::stod(value);
        break;
      case kDistanceMeanSimplex:
        feat.distance_mean_simplex = std::stod(value);
        break;
      case kFamilysizeSum:
        feat.familysize_sum = std::stoul(value);
        break;
      case kFamilysizeMean:
        feat.familysize_mean = std::stod(value);
        break;
      case kFamilysizeLT3Count:
        feat.familysize_lt3_count = std::stoul(value);
        break;
      case kFamilysizeLT5Count:
        feat.familysize_lt5_count = std::stoul(value);
        break;
      case kNonDuplex:
        feat.nonduplex = std::stoul(value);
        break;
      case kDuplex:
        feat.duplex = std::stoul(value);
        break;
      case kDuplexLowbq:
        feat.duplex_lowbq = std::stod(value);
        break;
      case kSimplex:
        feat.simplex = std::stoul(value);
        break;
      case kDuplexAF:
        feat.duplex_af = std::stod(value);
        break;
      case kPlusOnly:
        feat.plusonly = std::stoul(value);
        break;
      case kMinusOnly:
        feat.minusonly = std::stoul(value);
        break;
      case kWeightedScore:
        feat.weighted_score = std::stod(value);
        break;
      case kStrandBias:
        feat.strandbias = std::stod(value);
        break;
      case kContext:
        feat.context = DotToEmptyString(value);
        break;
      case kContextIndex:
        feat.context_index = std::stoul(value);
        break;
      case kMLScore:
        feat.ml_score = std::stod(value);
        break;
      case kFilterStatus:
        feat.filter_status = string::Split(value, ",");
        break;
      case kMapqLT60Ratio:
        feat.mapq_lt60_ratio = std::stod(value);
        break;
      case kMapqLT40Ratio:
        feat.mapq_lt40_ratio = std::stod(value);
        break;
      case kMapqLT30Ratio:
        feat.mapq_lt30_ratio = std::stod(value);
        break;
      case kMapqLT20Ratio:
        feat.mapq_lt20_ratio = std::stod(value);
        break;
      case kBaseqLT20Ratio:
        feat.baseq_lt20_ratio = std::stod(value);
        break;
      case kFamilysizeLT3Ratio:
        feat.familysize_lt3_ratio = std::stod(value);
        break;
      case kFamilysizeLT5Ratio:
        feat.familysize_lt5_ratio = std::stod(value);
        break;
      case kMQAF:
        feat.mq_af = std::stod(value);
        break;
      case kBQAF:
        feat.bq_af = std::stod(value);
        break;
      case kTumorSupport:
        feat.tumor_support = std::stoul(value);
        break;
      case kTumorMapqSum:
        feat.tumor_mapq_sum = std::stoul(value);
        break;
      case kTumorMapqMean:
        feat.tumor_mapq_mean = std::stod(value);
        break;
      case kTumorBaseqSum:
        feat.tumor_baseq_sum = std::stod(value);
        break;
      case kTumorBaseqMean:
        feat.tumor_mapq_mean = std::stod(value);
        break;
      case kTumorDistanceSum:
        feat.tumor_distance_sum = std::stoul(value);
        break;
      case kTumorDistanceMean:
        feat.tumor_distance_mean = std::stod(value);
        break;
      case kNormalSupport:
        feat.normal_support = std::stoul(value);
        break;
      case kNormalMapqSum:
        feat.normal_mapq_sum = std::stoul(value);
        break;
      case kNormalMapqMean:
        feat.normal_mapq_mean = std::stod(value);
        break;
      case kNormalBaseqSum:
        feat.normal_baseq_sum = std::stod(value);
        break;
      case kNormalBaseqMean:
        feat.normal_baseq_mean = std::stod(value);
        break;
      case kNormalDistanceSum:
        feat.normal_distance_sum = std::stoul(value);
        break;
      case kNormalDistanceMean:
        feat.normal_distance_mean = std::stod(value);
        break;
      case kBamTumorAF:
        feat.tumor_af = std::stod(value);
        break;
      case kBamNormalAF:
        feat.normal_af = std::stod(value);
        break;
      case kBamRAT:
        feat.rat = std::stod(value);
        break;
      case kSupportReverse:
        feat.support_reverse = std::stoul(value);
        break;
      case kTumorSupportReverse:
        feat.tumor_support_reverse = std::stoul(value);
        break;
      case kNormalSupportReverse:
        feat.normal_support_reverse = std::stoul(value);
        break;
      case kAlignmentBias:
        feat.alignmentbias = std::stod(value);
        break;
      case kTumorAlignmentBias:
        feat.tumor_alignmentbias = std::stod(value);
        break;
      case kNormalAlignmentBias:
        feat.normal_alignmentbias = std::stod(value);
        break;
      case kADT:
        feat.adt = std::stod(value);
        break;
      case kADTL:
        feat.adtl = std::stod(value);
        break;
      case kIndelAf:
        feat.indel_af = std::stod(value);
        break;
      default:
        break;
    }
  }
  return feat;
}

/**
 * @brief Deserialize a map of UnifiedFeatureCols to string into a UnifiedReferenceFeature struct
 * @param serialized a map of UnifiedFeatureCols enums to strings
 * @return a UnifiedReferenceFeature Struct
 */
UnifiedReferenceFeature VariantInfoSerializer::DeserializeRefFeature(
    const std::map<UnifiedFeatureCols, std::string>& serialized) {
  using enum UnifiedFeatureCols;
  UnifiedReferenceFeature feat{};
  for (const auto& [col, value] : serialized) {
    switch (col) {
      case kRefWeightedDepth:
        feat.weighted_depth = std::stod(value);
        break;
      case kRefNonhomopolymerWeightedDepth:
        feat.nonhomopolymer_weighted_depth = std::stod(value);
        break;
      case kRefSupport:
        feat.support = std::stoul(value);
        break;
      case kTumorRefSupport:
        feat.tumor_support = std::stoul(value);
        break;
      case kNormalRefSupport:
        feat.normal_support = std::stoul(value);
        break;
      case kRefNonhomopolymerSupport:
        feat.nonhomopolymer_support = std::stoul(value);
        break;
      case kRefMapqMin:
        feat.mapq_min = static_cast<u8>(std::stoul(value));
        break;
      case kRefMapqMax:
        feat.mapq_max = static_cast<u8>(std::stoul(value));
        break;
      case kRefMapqSum:
        feat.mapq_sum = std::stoul(value);
        break;
      case kRefMapqSumLowbq:
        feat.mapq_sum_lowbq = std::stoul(value);
        break;
      case kRefMapqSumSimplex:
        feat.mapq_sum_simplex = std::stoul(value);
        break;
      case kRefMapqMean:
        feat.mapq_mean = std::stod(value);
        break;
      case kRefMapqMeanLowbq:
        feat.mapq_mean_lowbq = std::stod(value);
        break;
      case kRefMapqMeanSimplex:
        feat.mapq_mean_simplex = std::stod(value);
        break;
      case kRefMapqLT60Count:
        feat.mapq_lt60_count = std::stoul(value);
        break;
      case kRefMapqLT40Count:
        feat.mapq_lt40_count = std::stoul(value);
        break;
      case kRefMapqLT30Count:
        feat.mapq_lt30_count = std::stoul(value);
        break;
      case kRefMapqLT20Count:
        feat.mapq_lt20_count = std::stoul(value);
        break;
      case kRefBaseqSum:
        feat.baseq_sum = std::stoul(value);
        break;
      case kRefBaseqLT20Count:
        feat.baseq_lt20_count = std::stoul(value);
        break;
      case kRefDistanceMin:
        feat.distance_min = std::stoul(value);
        break;
      case kRefDistanceMax:
        feat.distance_max = std::stoul(value);
        break;
      case kRefDistanceSum:
        feat.distance_sum = std::stoul(value);
        break;
      case kRefDistanceSumLowbq:
        feat.distance_sum_lowbq = std::stoul(value);
        break;
      case kRefDistanceSumSimplex:
        feat.distance_sum_simplex = std::stoul(value);
        break;
      case kRefDistanceMean:
        feat.distance_mean = std::stod(value);
        break;
      case kRefDistanceMeanLowbq:
        feat.distance_mean_lowbq = std::stod(value);
        break;
      case kRefDistanceMeanSimplex:
        feat.distance_mean_simplex = std::stod(value);
        break;
      case kRefDuplexLowbq:
        feat.duplex_lowbq = std::stod(value);
        break;
      case kRefSimplex:
        feat.simplex = std::stoul(value);
        break;
      case kRefDuplexAF:
        feat.duplex_af = std::stod(value);
        break;
      case kDuplexDP:
        feat.duplex_dp = std::stod(value);
        break;
      case kRefFamilysizeSum:
        feat.familysize_sum = std::stoul(value);
        break;
      case kRefFamilysizeLT3Count:
        feat.familysize_lt3_count = std::stoul(value);
        break;
      case kRefFamilysizeLT5Count:
        feat.familysize_lt5_count = std::stoul(value);
        break;
      case kRefNonhomopolymerMapqMin:
        feat.nonhomopolymer_mapq_min = static_cast<u8>(std::stoul(value));
        break;
      case kRefNonhomopolymerMapqMax:
        feat.nonhomopolymer_mapq_max = static_cast<u8>(std::stoul(value));
        break;
      case kRefNonhomopolymerMapqSum:
        feat.nonhomopolymer_mapq_sum = std::stoul(value);
        break;
      case kRefNonhomopolymerMapqLT60Count:
        feat.nonhomopolymer_mapq_lt60_count = std::stoul(value);
        break;
      case kRefNonhomopolymerMapqLT40Count:
        feat.nonhomopolymer_mapq_lt40_count = std::stoul(value);
        break;
      case kRefNonhomopolymerMapqLT30Count:
        feat.nonhomopolymer_mapq_lt30_count = std::stoul(value);
        break;
      case kRefNonhomopolymerMapqLT20Count:
        feat.nonhomopolymer_mapq_lt20_count = std::stoul(value);
        break;
      case kRefNonhomopolymerBaseqMin:
        feat.nonhomopolymer_baseq_min = static_cast<u8>(std::stoul(value));
        break;
      case kRefNonhomopolymerBaseqMax:
        feat.nonhomopolymer_baseq_max = static_cast<u8>(std::stoul(value));
        break;
      case kRefNonhomopolymerBaseqSum:
        feat.nonhomopolymer_baseq_sum = std::stoul(value);
        break;
      case kRefBaseqMean:
        feat.baseq_mean = std::stod(value);
        break;
      case kRefBaseqLT20Ratio:
        feat.baseq_lt20_ratio = std::stod(value);
        break;
      case kRefFamilysizeMean:
        feat.familysize_mean = std::stod(value);
        break;
      case kRefFamilysizeLT3Ratio:
        feat.familysize_lt3_ratio = std::stod(value);
        break;
      case kRefFamilysizeLT5Ratio:
        feat.familysize_lt5_ratio = std::stod(value);
        break;
      case kRefNonhomopolymerMapqMean:
        feat.nonhomopolymer_mapq_mean = std::stod(value);
        break;
      case kRefNonhomopolymerMapqLT60Ratio:
        feat.nonhomopolymer_mapq_lt60_ratio = std::stod(value);
        break;
      case kRefNonhomopolymerMapqLT40Ratio:
        feat.nonhomopolymer_mapq_lt40_ratio = std::stod(value);
        break;
      case kRefNonhomopolymerMapqLT30Ratio:
        feat.nonhomopolymer_mapq_lt30_ratio = std::stod(value);
        break;
      case kRefNonhomopolymerMapqLT20Ratio:
        feat.nonhomopolymer_mapq_lt20_ratio = std::stod(value);
        break;
      case kRefNonhomopolymerBaseqMean:
        feat.nonhomopolymer_baseq_mean = std::stod(value);
        break;
      case kRefMapqLT60Ratio:
        feat.mapq_lt60_ratio = std::stod(value);
        break;
      case kRefMapqLT40Ratio:
        feat.mapq_lt40_ratio = std::stod(value);
        break;
      case kRefMapqLT30Ratio:
        feat.mapq_lt30_ratio = std::stod(value);
        break;
      case kRefMapqLT20Ratio:
        feat.mapq_lt20_ratio = std::stod(value);
        break;
      case kRefMQAF:
        feat.mq_af = std::stod(value);
        break;
      case kRefBQAF:
        feat.bq_af = std::stod(value);
        break;
      case kNumAlt:
        feat.num_alt = std::stoul(value);
        break;
      default:
        break;
    }
  }
  return feat;
}

/**
 * @brief Deserialize a map of VcfFeatureCols to string into a VcfFeature struct
 * @param serialized a map of VcfFeatureCols to string
 * @return a VcfFeature struct
 */
VcfFeature VariantInfoSerializer::DeserializeVcfFeature(std::map<UnifiedFeatureCols, std::string>& serialized) {
  using enum UnifiedFeatureCols;
  VcfFeature feat{};
  for (const auto& [col, value] : serialized) {
    switch (col) {
      case kChrom:
        feat.chrom = value;
        break;
      case kPos:
        feat.pos = std::stoull(value) - 1;  // convert from 1-based to 0-based
        break;
      case kRef:
        feat.ref = value;
        break;
      case kAlt:
        feat.alt = value;
        break;
      case kVcfNalod:
        feat.nalod = std::stof(value);
        break;
      case kVcfNlod:
        feat.nlod = std::stof(value);
        break;
      case kVcfTlod:
        feat.tlod = std::stof(value);
        break;
      case kVcfMpos:
        feat.mpos = std::stoul(value);
        break;
      case kVcfMmqRef:
        feat.mmq_ref = static_cast<u8>(std::stoul(value));
        break;
      case kVcfMmqAlt:
        feat.mmq_alt = static_cast<u8>(std::stoul(value));
        break;
      case kVcfMbqRef:
        feat.mbq_ref = static_cast<u8>(std::stoul(value));
        break;
      case kVcfMbqAlt:
        feat.mbq_alt = static_cast<u8>(std::stoul(value));
        break;
      case kVcfPre2bContext:
        feat.pre_2bp_context = DotToEmptyString(value);
        break;
      case kVcfPost2bContext:
        feat.post_2bp_context = DotToEmptyString(value);
        break;
      case kVcfPost30bContext:
        feat.post_30bp_context = DotToEmptyString(value);
        break;
      case kVcfUnique3mers:
        feat.uniq_3mers = std::stoul(value);
        break;
      case kVcfUnique4mers:
        feat.uniq_4mers = std::stoul(value);
        break;
      case kVcfUnique5mers:
        feat.uniq_5mers = std::stoul(value);
        break;
      case kVcfUnique6mers:
        feat.uniq_6mers = std::stoul(value);
        break;
      case kVcfHomopolymer:
        feat.homopolymer = std::stoul(value);
        break;
      case kVcfDirepeat:
        feat.direpeat = std::stoul(value);
        break;
      case kVcfTrirepeat:
        feat.trirepeat = std::stoul(value);
        break;
      case kVcfQuadrepeat:
        feat.quadrepeat = std::stoul(value);
        break;
      case kVcfVariantDensity:
        feat.variant_density = std::stoul(value);
        break;
      case kVcfRefAd:
        feat.ref_ad = std::stoul(value);
        break;
      case kVcfAltAd:
        feat.alt_ad = std::stoul(value);
        break;
      case kVcfAltAd2:
        feat.alt_ad2 = std::stoul(value);
        break;
      case kVcfRefAdAf:
        feat.ref_ad_af = std::stod(value);
        break;
      case kVcfAltAdAf:
        feat.alt_ad_af = std::stod(value);
        break;
      case kVcfAltAd2Af:
        feat.alt_ad2_af = std::stod(value);
        break;
      case kVcfTumorAltAd:
        feat.tumor_alt_ad = std::stoul(value);
        break;
      case kVcfNormalAltAd:
        feat.normal_alt_ad = std::stoul(value);
        break;
      case kVcfTumorAf:
        feat.tumor_af = std::stof(value);
        break;
      case kVcfNormalAf:
        feat.normal_af = std::stof(value);
        break;
      case kVcfTumorNormalAfRatio:
        feat.tumor_normal_af_ratio = std::stof(value);
        break;
      case kVcfTumorDp:
        feat.tumor_dp = std::stoul(value);
        break;
      case kVcfNormalDp:
        feat.normal_dp = std::stoul(value);
        break;
      case kVcfPopAf:
        feat.popaf = std::stof(value);
        break;
      case kVcfGenotype:
        feat.genotype = StringToGenotype(value);
        break;
      case kVcfQual:
        feat.qual = std::stof(value);
        break;
      case kVcfGq:
        feat.gq = std::stoul(value);
        break;
      case kVcfHapcomp:
        feat.hapcomp = std::stoul(value);
        break;
      case kVcfHapdom:
        feat.hapdom = std::stof(value);
        break;
      case kVcfRu:
        feat.ru = DotToEmptyString(value);
        break;
      case kVcfRpaRef:
        feat.rpa_ref = std::stoul(value);
        break;
      case kVcfRpaAlt:
        feat.rpa_alt = std::stoul(value);
        break;
      case kVcfStr:
        feat.str = (std::stoul(value) != 0u);
        break;
      case kVcfAtInterest:
        feat.at_interest = (std::stoul(value) != 0u);
        break;
      default:
        break;
    }
  }
  return feat;
}

/**
 * @brief Serialize a VcfFeatures struct into a map of VcfFeatureCols to string
 * @param feature a VcfFeatures struct
 * @return a map of VcfFeatureCols to string
 */
std::map<UnifiedFeatureCols, std::string> VariantInfoSerializer::SerializeVcfFeature(const VcfFeature& feature) {
  using enum UnifiedFeatureCols;
  // serialize the VariantId part of the feature
  const VariantId vid(feature.chrom, feature.pos, feature.ref, feature.alt);
  std::map<UnifiedFeatureCols, std::string> result = SerializeVariantId(vid);
  // serialize the rest of the VcfFeature
  result[kVcfNalod] = std::to_string(feature.nalod);
  result[kVcfNlod] = std::to_string(feature.nlod);
  result[kVcfTlod] = std::to_string(feature.tlod);
  result[kVcfMpos] = std::to_string(feature.mpos);
  result[kVcfMmqRef] = std::to_string(feature.mmq_ref);
  result[kVcfMmqAlt] = std::to_string(feature.mmq_alt);
  result[kVcfMbqRef] = std::to_string(feature.mbq_ref);
  result[kVcfMbqAlt] = std::to_string(feature.mbq_alt);
  result[kVcfPre2bContext] = EmptyStringToDot(feature.pre_2bp_context);
  result[kVcfPost2bContext] = EmptyStringToDot(feature.post_2bp_context);
  result[kVcfPost30bContext] = EmptyStringToDot(feature.post_30bp_context);
  result[kVcfUnique3mers] = std::to_string(feature.uniq_3mers);
  result[kVcfUnique4mers] = std::to_string(feature.uniq_4mers);
  result[kVcfUnique5mers] = std::to_string(feature.uniq_5mers);
  result[kVcfUnique6mers] = std::to_string(feature.uniq_6mers);
  result[kVcfHomopolymer] = std::to_string(feature.homopolymer);
  result[kVcfDirepeat] = std::to_string(feature.direpeat);
  result[kVcfTrirepeat] = std::to_string(feature.trirepeat);
  result[kVcfQuadrepeat] = std::to_string(feature.quadrepeat);
  result[kVcfVariantDensity] = std::to_string(feature.variant_density);
  result[kVcfRefAd] = std::to_string(feature.ref_ad);
  result[kVcfAltAd] = std::to_string(feature.alt_ad);
  result[kVcfAltAd2] = std::to_string(feature.alt_ad2);
  result[kVcfRefAdAf] = std::to_string(feature.ref_ad_af);
  result[kVcfAltAdAf] = std::to_string(feature.alt_ad_af);
  result[kVcfAltAd2Af] = std::to_string(feature.alt_ad2_af);
  result[kVcfTumorAltAd] = std::to_string(feature.tumor_alt_ad);
  result[kVcfNormalAltAd] = std::to_string(feature.normal_alt_ad);
  result[kVcfTumorAf] = std::to_string(feature.tumor_af);
  result[kVcfNormalAf] = std::to_string(feature.normal_af);
  result[kVcfTumorNormalAfRatio] = std::to_string(feature.tumor_normal_af_ratio);
  result[kVcfTumorDp] = std::to_string(feature.tumor_dp);
  result[kVcfNormalDp] = std::to_string(feature.normal_dp);
  result[kVcfPopAf] = std::to_string(feature.popaf);
  result[kVcfGenotype] = GenotypeToString(feature.genotype);
  result[kVcfQual] = std::to_string(feature.qual);
  result[kVcfGq] = std::to_string(feature.gq);
  result[kVcfHapcomp] = std::to_string(feature.hapcomp);
  result[kVcfHapdom] = std::to_string(feature.hapdom);
  result[kVcfRu] = EmptyStringToDot(feature.ru);
  result[kVcfRpaRef] = std::to_string(feature.rpa_ref);
  result[kVcfRpaAlt] = std::to_string(feature.rpa_alt);
  result[kVcfStr] = feature.str ? "1" : "0";
  result[kVcfAtInterest] = feature.at_interest ? "1" : "0";
  return result;
}

/**
 * @brief Numericalize an attribute in a VariantId struct.
 * @param col Enum of the attribute
 * @param vid VariantId struct
 * @return double value of the attribute
 */
double VariantInfoSerializer::NumericalizeFeature(const UnifiedFeatureCols col, const VariantId& vid) {
  using enum UnifiedFeatureCols;
  switch (col) {
    case kChrom:
      // TODO: if possible, replace with htslib contig index
      return static_cast<double>(std::hash<std::string>{}(vid.chrom));
    case kPos:
      return static_cast<double>(vid.pos + 1);  // convert from 0-based to 1-based
    case kRef:
      return SeqToDouble(vid.ref);
    case kAlt:
      return SeqToDouble(vid.alt);
    case kSubTypeIndex:
      return SubstIndex(vid.ref, vid.alt);
    case kVariantType:
      return static_cast<double>(vid.type);
    case kAltLen:
      return static_cast<double>(vid.alt.size()) - static_cast<double>(vid.ref.size());
    default:
      return NAN;  // not an attribute of VariantId
  }
}

/**
 * @brief Numericalize an attribute in a UnifiedVariantFeature struct.
 * @param col Enum of the attribute
 * @param feat UnifiedVariantFeature struct
 * @return double value of the attribute
 */
double VariantInfoSerializer::NumericalizeFeature(const UnifiedFeatureCols col, const UnifiedVariantFeature& feat) {
  using enum UnifiedFeatureCols;
  switch (col) {
    case kWeightedDepth:
      return feat.weighted_depth;
    case kSupport:
      return feat.support;
    case kMapqMin:
      return feat.mapq_min;
    case kMapqMax:
      return feat.mapq_max;
    case kMapqSum:
      return feat.mapq_sum;
    case kMapqSumLowbq:
      return feat.mapq_sum_lowbq;
    case kMapqSumSimplex:
      return feat.mapq_sum_simplex;
    case kMapqMean:
      return feat.mapq_mean;
    case kMapqMeanLowbq:
      return feat.mapq_mean_lowbq;
    case kMapqMeanSimplex:
      return feat.mapq_mean_simplex;
    case kMapqLT60Count:
      return feat.mapq_lt60_count;
    case kMapqLT40Count:
      return feat.mapq_lt40_count;
    case kMapqLT30Count:
      return feat.mapq_lt30_count;
    case kMapqLT20Count:
      return feat.mapq_lt20_count;
    case kMapqLT60Ratio:
      return feat.mapq_lt60_ratio;
    case kMapqLT40Ratio:
      return feat.mapq_lt40_ratio;
    case kMapqLT30Ratio:
      return feat.mapq_lt30_ratio;
    case kMapqLT20Ratio:
      return feat.mapq_lt20_ratio;
    case kBaseqMin:
      return feat.baseq_min;
    case kBaseqMax:
      return feat.baseq_max;
    case kBaseqSum:
      return feat.baseq_sum;
    case kBaseqMean:
      return feat.baseq_mean;
    case kBaseqLT20Count:
      return feat.baseq_lt20_count;
    case kBaseqLT20Ratio:
      return feat.baseq_lt20_ratio;
    case kDistanceMin:
      return static_cast<double>(feat.distance_min);
    case kDistanceMax:
      return static_cast<double>(feat.distance_max);
    case kDistanceSum:
      return static_cast<double>(feat.distance_sum);
    case kDistanceSumLowbq:
      return static_cast<double>(feat.distance_sum_lowbq);
    case kDistanceSumSimplex:
      return static_cast<double>(feat.distance_sum_simplex);
    case kDistanceMean:
      return feat.distance_mean;
    case kDistanceMeanLowbq:
      return feat.distance_mean_lowbq;
    case kDistanceMeanSimplex:
      return feat.distance_mean_simplex;
    case kDuplexLowbq:
      return feat.duplex_lowbq;
    case kSimplex:
      return feat.simplex;
    case kDuplexAF:
      return feat.duplex_af;
    case kFamilysizeSum:
      return feat.familysize_sum;
    case kFamilysizeMean:
      return feat.familysize_mean;
    case kFamilysizeLT3Count:
      return feat.familysize_lt3_count;
    case kFamilysizeLT5Count:
      return feat.familysize_lt5_count;
    case kFamilysizeLT3Ratio:
      return feat.familysize_lt3_ratio;
    case kFamilysizeLT5Ratio:
      return feat.familysize_lt5_ratio;
    case kNonDuplex:
      return feat.nonduplex;
    case kDuplex:
      return feat.duplex;
    case kPlusOnly:
      return feat.plusonly;
    case kMinusOnly:
      return feat.minusonly;
    case kMQAF:
      return feat.mq_af;
    case kBQAF:
      return feat.bq_af;
    case kTumorSupport:
      return feat.tumor_support;
    case kTumorMapqSum:
      return feat.tumor_mapq_sum;
    case kTumorMapqMean:
      return feat.tumor_mapq_mean;
    case kTumorBaseqSum:
      return feat.tumor_baseq_sum;
    case kTumorBaseqMean:
      return feat.tumor_baseq_mean;
    case kTumorDistanceSum:
      return static_cast<double>(feat.tumor_distance_sum);
    case kTumorDistanceMean:
      return feat.tumor_distance_mean;
    case kNormalSupport:
      return feat.normal_support;
    case kNormalMapqSum:
      return feat.normal_mapq_sum;
    case kNormalMapqMean:
      return feat.normal_mapq_mean;
    case kNormalBaseqSum:
      return feat.normal_baseq_sum;
    case kNormalBaseqMean:
      return feat.normal_baseq_mean;
    case kNormalDistanceSum:
      return static_cast<double>(feat.normal_distance_sum);
    case kNormalDistanceMean:
      return feat.normal_distance_mean;
    case kBamTumorAF:
      return feat.tumor_af;
    case kBamNormalAF:
      return feat.normal_af;
    case kBamRAT:
      return feat.rat;
    case kSupportReverse:
      return feat.support_reverse;
    case kTumorSupportReverse:
      return feat.tumor_support_reverse;
    case kNormalSupportReverse:
      return feat.normal_support_reverse;
    case kAlignmentBias:
      return feat.alignmentbias;
    case kTumorAlignmentBias:
      return feat.tumor_alignmentbias;
    case kNormalAlignmentBias:
      return feat.normal_alignmentbias;
    case kWeightedScore:
      return feat.weighted_score;
    case kStrandBias:
      return feat.strandbias;
    case kContext:
      return UnifiedVariantFeature::ContextIndex(feat.context);
    case kContextIndex:
      return feat.context_index;
    case kMLScore:
      return feat.ml_score;
    case kFilterStatus:
      return static_cast<double>(util::hash::HashRange(feat.filter_status.begin(), feat.filter_status.end()));
    case kADT:
      return feat.adt;
    case kADTL:
      return feat.adtl;
    case kIndelAf:
      return feat.indel_af;
    default:
      return NAN;  // not an attribute of UnifiedVariantFeature
  }
}

/**
 * @brief Numericalize an attribute in a VcfFeature struct.
 * @param col Enum of the attribute
 * @param feat VcfFeature struct
 * @return double value of the attribute
 */
double VariantInfoSerializer::NumericalizeFeature(const UnifiedFeatureCols col, const VcfFeature& feat) {
  using enum UnifiedFeatureCols;
  switch (col) {
    case kChrom:
      // TODO: if possible, replace with htslib contig index
      return static_cast<double>(std::hash<std::string>{}(feat.chrom));
    case kPos:
      return static_cast<double>(feat.pos + 1);  // convert from 0-based to 1-based
    case kRef:
      return SeqToDouble(feat.ref);
    case kAlt:
      return SeqToDouble(feat.alt);
    case kVcfNalod:
      return feat.nalod;
    case kVcfNlod:
      return feat.nlod;
    case kVcfTlod:
      return feat.tlod;
    case kVcfMpos:
      return static_cast<double>(feat.mpos);
    case kVcfMmqRef:
      return feat.mmq_ref;
    case kVcfMmqAlt:
      return feat.mmq_alt;
    case kVcfMbqRef:
      return feat.mbq_ref;
    case kVcfMbqAlt:
      return feat.mbq_alt;
    case kVcfPre2bContext:
      if (IsAnyNotACTG(feat.pre_2bp_context)) {
        return UnifiedVariantFeature::ContextIndex("");
      } else {
        return UnifiedVariantFeature::ContextIndex(feat.pre_2bp_context);
      }
    case kVcfPost2bContext:
      if (IsAnyNotACTG(feat.post_2bp_context)) {
        return UnifiedVariantFeature::ContextIndex("");
      } else {
        return UnifiedVariantFeature::ContextIndex(feat.post_2bp_context);
      }
    case kVcfPost30bContext:
      return SeqToDouble(feat.post_30bp_context);
    case kVcfUnique3mers:
      return feat.uniq_3mers;
    case kVcfUnique4mers:
      return feat.uniq_4mers;
    case kVcfUnique5mers:
      return feat.uniq_5mers;
    case kVcfUnique6mers:
      return feat.uniq_6mers;
    case kVcfHomopolymer:
      return feat.homopolymer;
    case kVcfDirepeat:
      return feat.direpeat;
    case kVcfTrirepeat:
      return feat.trirepeat;
    case kVcfQuadrepeat:
      return feat.quadrepeat;
    case kVcfVariantDensity:
      return feat.variant_density;
    case kVcfRefAd:
      return feat.ref_ad;
    case kVcfAltAd:
      return feat.alt_ad;
    case kVcfAltAd2:
      return feat.alt_ad2;
    case kVcfRefAdAf:
      return feat.ref_ad_af;
    case kVcfAltAdAf:
      return feat.alt_ad_af;
    case kVcfAltAd2Af:
      return feat.alt_ad2_af;
    case kVcfTumorAltAd:
      return feat.tumor_alt_ad;
    case kVcfNormalAltAd:
      return feat.normal_alt_ad;
    case kVcfTumorAf:
      return feat.tumor_af;
    case kVcfNormalAf:
      return feat.normal_af;
    case kVcfTumorNormalAfRatio:
      return feat.tumor_normal_af_ratio;
    case kVcfTumorDp:
      return feat.tumor_dp;
    case kVcfNormalDp:
      return feat.normal_dp;
    case kVcfPopAf:
      return feat.popaf;
    case kVcfGenotype:
      return GenotypeToInt(feat.genotype);
    case kVcfQual:
      return feat.qual;
    case kVcfGq:
      return feat.gq;
    case kVcfHapcomp:
      return feat.hapcomp;
    case kVcfHapdom:
      return feat.hapdom;
    case kVcfRu:
      return SeqToDouble(feat.ru);
    case kVcfStr:
      return feat.str ? 1 : 0;
    case kVcfRpaRef:
      return feat.rpa_ref;
    case kVcfRpaAlt:
      return feat.rpa_alt;
    case kVcfAtInterest:
      return feat.at_interest ? 1 : 0;
    default:
      return NAN;  // not an attribute of VcfFeature
  }
}

/**
 * @brief Numericalize an attribute in a UnifiedReferenceFeature struct.
 * @param col Enum of the attribute
 * @param feat UnifiedReferenceFeature struct
 * @return double value of the attribute
 */
double VariantInfoSerializer::NumericalizeFeature(const UnifiedFeatureCols col, const UnifiedReferenceFeature& feat) {
  using enum UnifiedFeatureCols;
  switch (col) {
    case kRefSupport:
      return feat.support;
    case kRefWeightedDepth:
      return feat.weighted_depth;
    case kRefNonhomopolymerWeightedDepth:
      return feat.nonhomopolymer_weighted_depth;
    case kRefNonhomopolymerSupport:
      return feat.nonhomopolymer_support;
    case kRefMapqMin:
      return feat.mapq_min;
    case kRefMapqMax:
      return feat.mapq_max;
    case kRefMapqSum:
      return feat.mapq_sum;
    case kRefMapqSumLowbq:
      return feat.mapq_sum_lowbq;
    case kRefMapqSumSimplex:
      return feat.mapq_sum_simplex;
    case kRefMapqMean:
      return feat.mapq_mean;
    case kRefMapqMeanLowbq:
      return feat.mapq_mean_lowbq;
    case kRefMapqMeanSimplex:
      return feat.mapq_mean_simplex;
    case kRefMapqLT60Count:
      return feat.mapq_lt60_count;
    case kRefMapqLT40Count:
      return feat.mapq_lt40_count;
    case kRefMapqLT30Count:
      return feat.mapq_lt30_count;
    case kRefMapqLT20Count:
      return feat.mapq_lt20_count;
    case kRefBaseqSum:
      return feat.baseq_sum;
    case kRefBaseqMean:
      return feat.baseq_mean;
    case kRefBaseqLT20Count:
      return feat.baseq_lt20_count;
    case kRefBaseqLT20Ratio:
      return feat.baseq_lt20_ratio;
    case kRefDistanceMin:
      return static_cast<double>(feat.distance_min);
    case kRefDistanceMax:
      return static_cast<double>(feat.distance_max);
    case kRefDistanceSum:
      return static_cast<double>(feat.distance_sum);
    case kRefDistanceSumLowbq:
      return static_cast<double>(feat.distance_sum_lowbq);
    case kRefDistanceSumSimplex:
      return static_cast<double>(feat.distance_sum_simplex);
    case kRefDistanceMean:
      return feat.distance_mean;
    case kRefDistanceMeanLowbq:
      return feat.distance_mean_lowbq;
    case kRefDistanceMeanSimplex:
      return feat.distance_mean_simplex;
    case kRefDuplexLowbq:
      return feat.duplex_lowbq;
    case kRefSimplex:
      return feat.simplex;
    case kRefDuplexAF:
      return feat.duplex_af;
    case kDuplexDP:
      return feat.duplex_dp;
    case kRefFamilysizeSum:
      return feat.familysize_sum;
    case kRefFamilysizeMean:
      return feat.familysize_mean;
    case kRefFamilysizeLT3Count:
      return feat.familysize_lt3_count;
    case kRefFamilysizeLT5Count:
      return feat.familysize_lt5_count;
    case kRefFamilysizeLT3Ratio:
      return feat.familysize_lt3_ratio;
    case kRefFamilysizeLT5Ratio:
      return feat.familysize_lt5_ratio;
    case kRefNonhomopolymerMapqMin:
      return feat.nonhomopolymer_mapq_min;
    case kRefNonhomopolymerMapqMax:
      return feat.nonhomopolymer_mapq_max;
    case kRefNonhomopolymerMapqSum:
      return feat.nonhomopolymer_mapq_sum;
    case kRefNonhomopolymerMapqMean:
      return feat.nonhomopolymer_mapq_mean;
    case kRefNonhomopolymerMapqLT60Count:
      return feat.nonhomopolymer_mapq_lt60_count;
    case kRefNonhomopolymerMapqLT40Count:
      return feat.nonhomopolymer_mapq_lt40_count;
    case kRefNonhomopolymerMapqLT30Count:
      return feat.nonhomopolymer_mapq_lt30_count;
    case kRefNonhomopolymerMapqLT20Count:
      return feat.nonhomopolymer_mapq_lt20_count;
    case kRefNonhomopolymerMapqLT60Ratio:
      return feat.nonhomopolymer_mapq_lt60_ratio;
    case kRefNonhomopolymerMapqLT40Ratio:
      return feat.nonhomopolymer_mapq_lt40_ratio;
    case kRefNonhomopolymerMapqLT30Ratio:
      return feat.nonhomopolymer_mapq_lt30_ratio;
    case kRefNonhomopolymerMapqLT20Ratio:
      return feat.nonhomopolymer_mapq_lt20_ratio;
    case kRefNonhomopolymerBaseqMin:
      return feat.nonhomopolymer_baseq_min;
    case kRefNonhomopolymerBaseqMax:
      return feat.nonhomopolymer_baseq_max;
    case kRefNonhomopolymerBaseqSum:
      return feat.nonhomopolymer_baseq_sum;
    case kRefNonhomopolymerBaseqMean:
      return feat.nonhomopolymer_baseq_mean;
    case kRefMapqLT60Ratio:
      return feat.mapq_lt60_ratio;
    case kRefMapqLT40Ratio:
      return feat.mapq_lt40_ratio;
    case kRefMapqLT30Ratio:
      return feat.mapq_lt30_ratio;
    case kRefMapqLT20Ratio:
      return feat.mapq_lt20_ratio;
    case kRefMQAF:
      return feat.mq_af;
    case kRefBQAF:
      return feat.bq_af;
    case kNumAlt:
      return feat.num_alt;
    case kTumorRefSupport:
      return feat.tumor_support;
    case kNormalRefSupport:
      return feat.normal_support;
    default:
      return NAN;  // not an attribute of UnifiedReferenceFeature
  }
}

/**
 * @brief Numericalize an attribute in the feature structs.
 * @param col Enum of the attribute
 * @param vid VariantId struct
 * @param bam_feat UnifiedVariantFeature struct
 * @param vcf_feat VcfFeature struct
 * @param ref_feat UnifiedReferenceFeature struct
 * @return double value of the attribute
 */
double VariantInfoSerializer::NumericalizeFeature(const UnifiedFeatureCols col,
                                                  const VariantId& vid,
                                                  const UnifiedVariantFeature& bam_feat,
                                                  const VcfFeature& vcf_feat,
                                                  const UnifiedReferenceFeature& ref_feat) {
  auto val = NumericalizeFeature(col, vid);
  if (std::isnan(val)) {
    val = NumericalizeFeature(col, bam_feat);
    if (std::isnan(val)) {
      val = NumericalizeFeature(col, vcf_feat);
      if (std::isnan(val)) {
        val = NumericalizeFeature(col, ref_feat);
      }
    }
  }
  return val;
}

/**
 * Load BAM features from a features file.
 * @param features_file a text file containing BAM features
 * @return a pair of maps, one for UnifiedVariantFeatures, the other for UnifiedReferenceFeatures
 */
std::pair<ChromToVariantInfoMap, RefInfoMap> VariantInfoSerializer::LoadFeatures(const std::string& features_file) {
  ChromToVariantInfoMap features;
  RefInfoMap ref_features;
  std::ifstream input(features_file);
  std::string line;
  if (!std::getline(input, line)) {
    throw error::Error("Could not read header from features file " + features_file);
  }
  // Extract field names from header line.
  auto header_names = string::Split(line, "\t");
  if (header_names.empty()) {
    throw error::Error("Cannot find feature names in feature file: {}", features_file);
  }
  auto unsupported_fields = FindUnsupportedBamFeatureNames(header_names);
  if (!unsupported_fields.empty()) {
    throw error::Error("Unsupported BAM feature names: {}", string::Join(unsupported_fields, ", "));
  }
  // Convert header field names from strings to enums.
  vec<UnifiedFeatureCols> header_cols;
  header_cols.reserve(header_names.size());
  for (auto& f : header_names) {
    header_cols.emplace_back(GetCol(f));
  }
  while (std::getline(input, line)) {
    auto field_values = string::Split(line, "\t");
    if (field_values.size() != header_cols.size()) {
      throw error::Error("Number of fields in features file " + features_file + " does not match header");
    }
    std::map<UnifiedFeatureCols, std::string> row;
    for (size_t i = 0; i < header_cols.size(); i++) {
      row[header_cols[i]] = field_values[i];
    }
    const VariantId id{DeserializeVariantId(row)};
    features[id.chrom][id.pos][id] = DeserializeVariantFeature(row);
    ref_features[id.chrom][id.GetRefFeaturePos()] = DeserializeRefFeature(row);
  }
  return std::make_pair(features, ref_features);
}

/**
 * @brief Load BAM features from a features file for variants at VCF feature positions.
 * @param features_file Path to BAM features file
 * @param chr_to_vcf_feat_map Map of chromosomal positions to VCF features
 * @return Pair of maps, one for UnifiedVariantFeatures, the other for UnifiedReferenceFeatures
 */
std::pair<ChromToVariantInfoMap, RefInfoMap> VariantInfoSerializer::LoadFeatures(
    const std::string& features_file, const ChromToVcfFeaturesMap& chr_to_vcf_feat_map) {
  ChromToVariantInfoMap features;
  RefInfoMap ref_features;
  std::ifstream input(features_file);
  std::string line;
  if (!std::getline(input, line)) {
    throw error::Error("Could not read header from features file " + features_file);
  }
  // Extract field names from header line.
  auto header_names = string::Split(line, "\t");
  if (header_names.empty()) {
    throw error::Error("Cannot find feature names in feature file: {}", features_file);
  }
  auto unsupported_fields = FindUnsupportedBamFeatureNames(header_names);
  if (!unsupported_fields.empty()) {
    throw error::Error("Unsupported BAM feature names: {}", string::Join(unsupported_fields, ", "));
  }
  // Convert header field names from strings to enums.
  vec<UnifiedFeatureCols> header_cols;
  header_cols.reserve(header_names.size());
  for (auto& f : header_names) {
    header_cols.emplace_back(GetCol(f));
  }
  const auto num_cols = header_cols.size();
  while (std::getline(input, line)) {
    auto field_values = string::Split(line, "\t");
    if (field_values.size() != num_cols) {
      throw error::Error("Number of fields in features file " + features_file + " does not match header");
    }
    std::map<UnifiedFeatureCols, std::string> row;
    for (size_t i = 0; i < num_cols; i++) {
      row[header_cols[i]] = field_values[i];
    }
    const VariantId id{DeserializeVariantId(row)};
    // only keep BAM feature at positions containing VCF feature(s)
    auto itr = chr_to_vcf_feat_map.find(id.chrom);
    if (itr != chr_to_vcf_feat_map.end() && itr->second.contains(id.pos)) {
      // Store the BAM feature only if it corresponds to a VCF feature position.
      features[id.chrom][id.pos][id] = DeserializeVariantFeature(row);
      ref_features[id.chrom][id.GetRefFeaturePos()] = DeserializeRefFeature(row);
    }
  }
  return std::make_pair(features, ref_features);
}

ChromToVcfFeaturesMap VariantInfoSerializer::LoadVcfFeatures(const std::string& features_file) {
  ChromToVcfFeaturesMap features;
  std::ifstream input(features_file);
  std::string line;
  if (!std::getline(input, line)) {
    throw error::Error("Could not read header from VCF features file " + features_file);
  }
  // Extract field names from header line.
  auto header_names = string::Split(line, "\t");
  if (header_names.empty()) {
    throw error::Error("Cannot find feature names in feature file: {}", features_file);
  }
  auto unsupported_fields = FindUnsupportedVcfFeatureNames(header_names);
  if (!unsupported_fields.empty()) {
    throw error::Error("Unsupported VCF feature names: {}", string::Join(unsupported_fields, ", "));
  }
  // Convert header field names from strings to enums.
  vec<UnifiedFeatureCols> header_cols;
  header_cols.reserve(header_names.size());
  for (auto& f : header_names) {
    header_cols.emplace_back(GetCol(f));
  }
  while (std::getline(input, line)) {
    auto field_values = string::Split(line, "\t");
    if (field_values.size() != header_cols.size()) {
      throw error::Error("Number of fields in VCF features file " + features_file + " does not match header");
    }
    std::map<UnifiedFeatureCols, std::string> row;
    for (size_t i = 0; i < header_cols.size(); i++) {
      row[header_cols[i]] = field_values[i];
    }
    const VariantId id = DeserializeVariantId(row);
    features[id.chrom][id.pos][id] = DeserializeVcfFeature(row);
  }
  return features;
}

}  // namespace xoos::svc
