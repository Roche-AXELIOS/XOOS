#include "variant-info-serializer.h"

#include <cmath>
#include <fstream>
#include <limits>
#include <utility>
#include <vector>

#include <magic_enum/magic_enum.hpp>
#include <nlohmann/detail/input/parser.hpp>

#include <xoos/error/error.h>
#include <xoos/util/hash.h>
#include <xoos/util/parse-int.h>
#include <xoos/util/string-functions.h>

#include "core/column-names.h"
#include "core/genotype.h"
#include "util/seq-util.h"
#include "xoos/io/metadata-util.h"

namespace xoos::svc {

namespace {
/**
 * @brief Convert string to a dot if it is empty, otherwise return the string as is.
 * @param value Input string value
 * @return Converted string, empty strings become ".", otherwise the original string
 */
std::string EmptyStringToDot(const std::string& value) {
  return value.empty() ? "." : value;
}

/**
 * @brief Convert a dot to an empty string if the input is ".", otherwise return the string as is.
 * @param value Input string value
 * @return Converted string, "." becomes empty string, otherwise the original string
 */
std::string DotToEmptyString(const std::string& value) {
  return value == "." ? "" : value;
}

/**
 * @brief Parse an RPA (Repeat Per Allele) value with backward compatibility.
 *
 * Feature TSVs produced by older versions serialized RPA as u32, meaning negative values
 * were written as large unsigned integers (e.g., -1 became "4294967295"). This function
 * accepts both the current signed format (e.g., "-1") and the legacy unsigned format, mapping
 * values in the range (INT32_MAX, UINT32_MAX] back to their signed s32 equivalents.
 *
 * @param value String representation of an RPA value
 * @return Parsed s32 value
 * @throws xoos::error::Error if the value is out of u32 range or otherwise invalid
 */
s32 ParseRpaValue(const std::string_view value) {
  // ParseS64 accepts both negative strings (e.g., "-1") and large positive strings (e.g., "4294967295").
  s64 parsed{};
  try {
    parsed = util::ParseS64(value);
  } catch (const std::exception&) {
    throw error::Error("Cannot parse RPA value '{}': not a valid integer", value);
  }
  if (parsed >= std::numeric_limits<s32>::min() && parsed <= std::numeric_limits<s32>::max()) {
    return static_cast<s32>(parsed);
  }
  // Value is outside s32 range. If it fits in u32 it is a legacy underflowed value;
  // reinterpret the bit pattern as s32 to recover the original negative value.
  if (parsed >= 0 && std::cmp_less_equal(parsed, std::numeric_limits<u32>::max())) {
    return static_cast<s32>(static_cast<u32>(parsed));
  }
  throw error::Error("RPA value '{}' is out of range", value);
}
}  // namespace

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
  result[kDuplexConcordant] = std::to_string(variant_info.duplex_concordant);
  result[kDuplexSimplex] = std::to_string(variant_info.duplex_simplex);
  result[kDuplexDiscordant] = std::to_string(variant_info.duplex_discordant);
  result[kSimplex] = std::to_string(variant_info.simplex);
  result[kDuplexAF] = std::to_string(variant_info.duplex_af);
  result[kDuplexOther] = std::to_string(variant_info.duplex_other);
  result[kDuplexOtherAF] = std::to_string(variant_info.duplex_other_af);
  result[kDuplexGQ] = std::to_string(variant_info.duplex_gq);
  result[kAfqGT01] = std::to_string(variant_info.afq_gt01);
  result[kAfqGT11] = std::to_string(variant_info.afq_gt11);
  result[kAfqGT12] = std::to_string(variant_info.afq_gt12);
  result[kAfqGT] = GenotypeToString(variant_info.afq_gt);
  result[kPlusOnly] = std::to_string(variant_info.plusonly);
  result[kMinusOnly] = std::to_string(variant_info.minusonly);
  result[kMQAF] = std::to_string(variant_info.mq_af);
  result[kBQAF] = std::to_string(variant_info.bq_af);
  result[kBamTnAfRatio] = std::to_string(variant_info.tn_af_ratio);
  result[kSupportReverse] = std::to_string(variant_info.support_reverse);
  result[kAlignmentBias] = std::to_string(variant_info.alignmentbias);
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
  result[kDel1bpAF] = std::to_string(ref_info.del_1bp_af);
  result[kDel2bpAF] = std::to_string(ref_info.del_2bp_af);
  result[kDel3bpAF] = std::to_string(ref_info.del_3bp_af);
  result[kDelGt3bpAF] = std::to_string(ref_info.del_gt3bp_af);
  result[kIns1bpAF] = std::to_string(ref_info.ins_1bp_af);
  result[kIns2bpAF] = std::to_string(ref_info.ins_2bp_af);
  result[kIns3bpAF] = std::to_string(ref_info.ins_3bp_af);
  result[kInsGt3bpAF] = std::to_string(ref_info.ins_gt3bp_af);
  result[kDel1bp2bpAFLogRatio] = std::to_string(ref_info.del_1bp_2bp_af_logratio);
  result[kDel2bp3bpAFLogRatio] = std::to_string(ref_info.del_2bp_3bp_af_logratio);
  result[kDel3bpGt3bpAFLogRatio] = std::to_string(ref_info.del_3bp_gt3bp_af_logratio);
  result[kIns1bp2bpAFLogRatio] = std::to_string(ref_info.ins_1bp_2bp_af_logratio);
  result[kIns2bp3bpAFLogRatio] = std::to_string(ref_info.ins_2bp_3bp_af_logratio);
  result[kIns3bpGt3bpAFLogRatio] = std::to_string(ref_info.ins_3bp_gt3bp_af_logratio);
  result[kDel1bpIns1bpAFLogRatio] = std::to_string(ref_info.del_1bp_ins_1bp_af_logratio);
  result[kSnvAF] = std::to_string(ref_info.snv_af);
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
  result[kRefDuplexConcordant] = std::to_string(ref_info.duplex_concordant);
  result[kRefDuplexSimplex] = std::to_string(ref_info.duplex_simplex);
  result[kRefDuplexDiscordant] = std::to_string(ref_info.duplex_discordant);
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
  result[kRefMQAF] = std::to_string(ref_info.mq_af);
  result[kRefBQAF] = std::to_string(ref_info.bq_af);
  result[kNumAlt] = std::to_string(ref_info.num_alt);
  return result;
}

VariantId VariantInfoSerializer::DeserializeVariantId(const vec<FeatureColumn>& header,
                                                      const vec<std::string>& fields) {
  using enum UnifiedFeatureCols;
  std::string chrom;
  u64 pos = 0;
  std::string ref;
  std::string alt;
  for (size_t i = 0; i < header.size(); ++i) {
    const auto& col = header[i];
    const auto& value = fields[i];
    switch (col.enum_val) {
      case kChrom:
        chrom = value;
        break;
      case kPos:
        // convert from 1-based to 0-based
        pos = value.empty() ? 0 : static_cast<u64>(std::stol(value) - 1);
        break;
      case kRef:
        ref = value;
        break;
      case kAlt:
        alt = value;
        break;
      default:
        break;
    }
  }
  // Other fields in VariantId are derived from these four primary fields.
  // So, we construct VariantId after extracting these four fields.
  const auto& vid = VariantId(chrom, pos, ref, alt);
  return vid;
}

namespace {
bool UpdateVariantBamFeature(const FeatureColumn& col, const std::string& value, UnifiedVariantFeature& feat) {
  using enum UnifiedFeatureCols;
  switch (col.enum_val) {
    case kWeightedDepth:
      feat.weighted_depth = std::stod(value);
      break;
    case kSupport:
      feat.support = util::ParseU32(value);
      break;
    case kMapqMax:
      feat.mapq_max = util::ParseU8(value);
      break;
    case kMapqMin:
      feat.mapq_min = util::ParseU8(value);
      break;
    case kMapqSum:
      feat.mapq_sum = util::ParseU32(value);
      break;
    case kMapqSumLowbq:
      feat.mapq_sum_lowbq = util::ParseU32(value);
      break;
    case kMapqSumSimplex:
      feat.mapq_sum_simplex = util::ParseU32(value);
      break;
    case kMapqMean:
      feat.mapq_mean = std::stof(value);
      break;
    case kMapqMeanLowbq:
      feat.mapq_mean_lowbq = std::stof(value);
      break;
    case kMapqMeanSimplex:
      feat.mapq_mean_simplex = std::stof(value);
      break;
    case kMapqLT60Count:
      feat.mapq_lt60_count = util::ParseU32(value);
      break;
    case kMapqLT40Count:
      feat.mapq_lt40_count = util::ParseU32(value);
      break;
    case kMapqLT30Count:
      feat.mapq_lt30_count = util::ParseU32(value);
      break;
    case kMapqLT20Count:
      feat.mapq_lt20_count = util::ParseU32(value);
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
      feat.baseq_mean = std::stof(value);
      break;
    case kBaseqLT20Count:
      feat.baseq_lt20_count = util::ParseU32(value);
      break;
    case kDistanceMin:
      feat.distance_min = util::ParseU64(value);
      break;
    case kDistanceMax:
      feat.distance_max = util::ParseU64(value);
      break;
    case kDistanceSum:
      feat.distance_sum = util::ParseU64(value);
      break;
    case kDistanceSumLowbq:
      feat.distance_sum_lowbq = util::ParseU64(value);
      break;
    case kDistanceSumSimplex:
      feat.distance_sum_simplex = util::ParseU64(value);
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
      feat.familysize_sum = util::ParseU32(value);
      break;
    case kFamilysizeMean:
      feat.familysize_mean = std::stof(value);
      break;
    case kFamilysizeLT3Count:
      feat.familysize_lt3_count = util::ParseU32(value);
      break;
    case kFamilysizeLT5Count:
      feat.familysize_lt5_count = util::ParseU32(value);
      break;
    case kNonDuplex:
      feat.nonduplex = util::ParseU32(value);
      break;
    case kDuplex:
      feat.duplex = util::ParseU32(value);
      break;
    case kDuplexLowbq:
      feat.duplex_lowbq = std::stod(value);
      break;
    case kDuplexConcordant:
      feat.duplex_concordant = util::ParseU32(value);
      break;
    case kDuplexSimplex:
      feat.duplex_simplex = util::ParseU32(value);
      break;
    case kDuplexDiscordant:
      feat.duplex_discordant = util::ParseU32(value);
      break;
    case kSimplex:
      feat.simplex = util::ParseU32(value);
      break;
    case kDuplexAF:
      feat.duplex_af = std::stof(value);
      break;
    case kDuplexOther:
      feat.duplex_other = std::stod(value);
      break;
    case kDuplexOtherAF:
      feat.duplex_other_af = std::stof(value);
      break;
    case kDuplexGQ:
      feat.duplex_gq = std::stod(value);
      break;
    case kAfqGT01:
      feat.afq_gt01 = std::stod(value);
      break;
    case kAfqGT11:
      feat.afq_gt11 = std::stod(value);
      break;
    case kAfqGT12:
      feat.afq_gt12 = std::stod(value);
      break;
    case kAfqGT:
      feat.afq_gt = StringToGenotype(value);
      break;
    case kPlusOnly:
      feat.plusonly = util::ParseU32(value);
      break;
    case kMinusOnly:
      feat.minusonly = util::ParseU32(value);
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
      feat.context_index = util::ParseU32(value);
      break;
    case kMLScore:
      feat.ml_score = std::stod(value);
      break;
    case kFilterStatus:
      feat.filter_status = string::Split(value, ",");
      break;
    case kMapqLT60Ratio:
      feat.mapq_lt60_ratio = std::stof(value);
      break;
    case kMapqLT40Ratio:
      feat.mapq_lt40_ratio = std::stof(value);
      break;
    case kMapqLT30Ratio:
      feat.mapq_lt30_ratio = std::stof(value);
      break;
    case kMapqLT20Ratio:
      feat.mapq_lt20_ratio = std::stof(value);
      break;
    case kBaseqLT20Ratio:
      feat.baseq_lt20_ratio = std::stof(value);
      break;
    case kFamilysizeLT3Ratio:
      feat.familysize_lt3_ratio = std::stof(value);
      break;
    case kFamilysizeLT5Ratio:
      feat.familysize_lt5_ratio = std::stof(value);
      break;
    case kMQAF:
      feat.mq_af = std::stof(value);
      break;
    case kBQAF:
      feat.bq_af = std::stof(value);
      break;
    case kBamTnAfRatio:
      feat.tn_af_ratio = std::stof(value);
      break;
    case kSupportReverse:
      feat.support_reverse = util::ParseU32(value);
      break;
    case kAlignmentBias:
      feat.alignmentbias = std::stod(value);
      break;
    case kADT:
      feat.adt = std::stod(value);
      break;
    case kADTL:
      feat.adtl = std::stod(value);
      break;
    case kIndelAf:
      feat.indel_af = std::stof(value);
      break;
    default:
      return false;
  }
  return true;
}

bool UpdateReferenceBamFeature(const FeatureColumn& col, const std::string& value, UnifiedReferenceFeature& feat) {
  using enum UnifiedFeatureCols;
  switch (col.enum_val) {
    case kRefWeightedDepth:
      feat.weighted_depth = std::stod(value);
      break;
    case kRefNonhomopolymerWeightedDepth:
      feat.nonhomopolymer_weighted_depth = std::stod(value);
      break;
    case kRefSupport:
      feat.support = util::ParseU32(value);
      break;
    case kRefNonhomopolymerSupport:
      feat.nonhomopolymer_support = util::ParseU32(value);
      break;
    case kRefMapqMin:
      feat.mapq_min = util::ParseU8(value);
      break;
    case kRefMapqMax:
      feat.mapq_max = util::ParseU8(value);
      break;
    case kRefMapqSum:
      feat.mapq_sum = util::ParseU32(value);
      break;
    case kRefMapqSumLowbq:
      feat.mapq_sum_lowbq = util::ParseU32(value);
      break;
    case kRefMapqSumSimplex:
      feat.mapq_sum_simplex = util::ParseU32(value);
      break;
    case kRefMapqMean:
      feat.mapq_mean = std::stof(value);
      break;
    case kRefMapqMeanLowbq:
      feat.mapq_mean_lowbq = std::stof(value);
      break;
    case kRefMapqMeanSimplex:
      feat.mapq_mean_simplex = std::stof(value);
      break;
    case kRefMapqLT60Count:
      feat.mapq_lt60_count = util::ParseU32(value);
      break;
    case kRefMapqLT40Count:
      feat.mapq_lt40_count = util::ParseU32(value);
      break;
    case kRefMapqLT30Count:
      feat.mapq_lt30_count = util::ParseU32(value);
      break;
    case kRefMapqLT20Count:
      feat.mapq_lt20_count = util::ParseU32(value);
      break;
    case kRefBaseqSum:
      feat.baseq_sum = util::ParseU32(value);
      break;
    case kRefBaseqLT20Count:
      feat.baseq_lt20_count = util::ParseU32(value);
      break;
    case kRefDistanceMin:
      feat.distance_min = util::ParseU64(value);
      break;
    case kRefDistanceMax:
      feat.distance_max = util::ParseU64(value);
      break;
    case kRefDistanceSum:
      feat.distance_sum = util::ParseU64(value);
      break;
    case kRefDistanceSumLowbq:
      feat.distance_sum_lowbq = util::ParseU64(value);
      break;
    case kRefDistanceSumSimplex:
      feat.distance_sum_simplex = util::ParseU64(value);
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
    case kRefDuplexConcordant:
      feat.duplex_concordant = util::ParseU32(value);
      break;
    case kRefDuplexSimplex:
      feat.duplex_simplex = util::ParseU32(value);
      break;
    case kRefDuplexDiscordant:
      feat.duplex_discordant = util::ParseU32(value);
      break;
    case kRefSimplex:
      feat.simplex = util::ParseU32(value);
      break;
    case kRefDuplexAF:
      feat.duplex_af = std::stof(value);
      break;
    case kDuplexDP:
      feat.duplex_dp = std::stod(value);
      break;
    case kRefFamilysizeSum:
      feat.familysize_sum = util::ParseU32(value);
      break;
    case kRefFamilysizeLT3Count:
      feat.familysize_lt3_count = util::ParseU32(value);
      break;
    case kRefFamilysizeLT5Count:
      feat.familysize_lt5_count = util::ParseU32(value);
      break;
    case kRefNonhomopolymerMapqMin:
      feat.nonhomopolymer_mapq_min = util::ParseU8(value);
      break;
    case kRefNonhomopolymerMapqMax:
      feat.nonhomopolymer_mapq_max = util::ParseU8(value);
      break;
    case kRefNonhomopolymerMapqSum:
      feat.nonhomopolymer_mapq_sum = util::ParseU32(value);
      break;
    case kRefNonhomopolymerMapqLT60Count:
      feat.nonhomopolymer_mapq_lt60_count = util::ParseU32(value);
      break;
    case kRefNonhomopolymerMapqLT40Count:
      feat.nonhomopolymer_mapq_lt40_count = util::ParseU32(value);
      break;
    case kRefNonhomopolymerMapqLT30Count:
      feat.nonhomopolymer_mapq_lt30_count = util::ParseU32(value);
      break;
    case kRefNonhomopolymerMapqLT20Count:
      feat.nonhomopolymer_mapq_lt20_count = util::ParseU32(value);
      break;
    case kRefNonhomopolymerBaseqMin:
      feat.nonhomopolymer_baseq_min = util::ParseU8(value);
      break;
    case kRefNonhomopolymerBaseqMax:
      feat.nonhomopolymer_baseq_max = util::ParseU8(value);
      break;
    case kRefNonhomopolymerBaseqSum:
      feat.nonhomopolymer_baseq_sum = util::ParseU32(value);
      break;
    case kRefBaseqMean:
      feat.baseq_mean = std::stof(value);
      break;
    // AF and log-ratio fields are stored as f32; parse with std::stof to avoid a narrowing f64->f32 conversion.
    case kDel1bpAF:
      feat.del_1bp_af = std::stof(value);
      break;
    case kDel2bpAF:
      feat.del_2bp_af = std::stof(value);
      break;
    case kDel3bpAF:
      feat.del_3bp_af = std::stof(value);
      break;
    case kDelGt3bpAF:
      feat.del_gt3bp_af = std::stof(value);
      break;
    case kIns1bpAF:
      feat.ins_1bp_af = std::stof(value);
      break;
    case kIns2bpAF:
      feat.ins_2bp_af = std::stof(value);
      break;
    case kIns3bpAF:
      feat.ins_3bp_af = std::stof(value);
      break;
    case kInsGt3bpAF:
      feat.ins_gt3bp_af = std::stof(value);
      break;
    case kDel1bp2bpAFLogRatio:
      feat.del_1bp_2bp_af_logratio = std::stof(value);
      break;
    case kDel2bp3bpAFLogRatio:
      feat.del_2bp_3bp_af_logratio = std::stof(value);
      break;
    case kDel3bpGt3bpAFLogRatio:
      feat.del_3bp_gt3bp_af_logratio = std::stof(value);
      break;
    case kIns1bp2bpAFLogRatio:
      feat.ins_1bp_2bp_af_logratio = std::stof(value);
      break;
    case kIns2bp3bpAFLogRatio:
      feat.ins_2bp_3bp_af_logratio = std::stof(value);
      break;
    case kIns3bpGt3bpAFLogRatio:
      feat.ins_3bp_gt3bp_af_logratio = std::stof(value);
      break;
    case kDel1bpIns1bpAFLogRatio:
      feat.del_1bp_ins_1bp_af_logratio = std::stof(value);
      break;
    case kSnvAF:
      feat.snv_af = std::stof(value);
      break;
    case kRefBaseqLT20Ratio:
      feat.baseq_lt20_ratio = std::stof(value);
      break;
    case kRefFamilysizeMean:
      feat.familysize_mean = std::stof(value);
      break;
    case kRefFamilysizeLT3Ratio:
      feat.familysize_lt3_ratio = std::stof(value);
      break;
    case kRefFamilysizeLT5Ratio:
      feat.familysize_lt5_ratio = std::stof(value);
      break;
    case kRefNonhomopolymerMapqMean:
      feat.nonhomopolymer_mapq_mean = std::stof(value);
      break;
    case kRefNonhomopolymerMapqLT60Ratio:
      feat.nonhomopolymer_mapq_lt60_ratio = std::stof(value);
      break;
    case kRefNonhomopolymerMapqLT40Ratio:
      feat.nonhomopolymer_mapq_lt40_ratio = std::stof(value);
      break;
    case kRefNonhomopolymerMapqLT30Ratio:
      feat.nonhomopolymer_mapq_lt30_ratio = std::stof(value);
      break;
    case kRefNonhomopolymerMapqLT20Ratio:
      feat.nonhomopolymer_mapq_lt20_ratio = std::stof(value);
      break;
    case kRefNonhomopolymerBaseqMean:
      feat.nonhomopolymer_baseq_mean = std::stof(value);
      break;
    case kRefMapqLT60Ratio:
      feat.mapq_lt60_ratio = std::stof(value);
      break;
    case kRefMapqLT40Ratio:
      feat.mapq_lt40_ratio = std::stof(value);
      break;
    case kRefMapqLT30Ratio:
      feat.mapq_lt30_ratio = std::stof(value);
      break;
    case kRefMapqLT20Ratio:
      feat.mapq_lt20_ratio = std::stof(value);
      break;
    case kRefMQAF:
      feat.mq_af = std::stof(value);
      break;
    case kRefBQAF:
      feat.bq_af = std::stof(value);
      break;
    case kNumAlt:
      feat.num_alt = util::ParseU32(value);
      break;
    default:
      return false;
  }
  return true;
}
}  // namespace

vec<std::string> VariantInfoSerializer::SerializeBamFeatureRow(const vec<FeatureColumn>& header,
                                                               const VariantId& vid,
                                                               const BamFeatureTuple& feat_row) {
  vec<std::string> result(header.size());
  // Serialize the variant ID and features into maps for easy lookup
  const auto vid_map = SerializeVariantId(vid);
  const auto var_feat_map = SerializeVariantFeature(feat_row.var_feat);
  const auto ref_feat_map = SerializeReferenceFeature(feat_row.ref_feat);

  // Fill in the result based on the feature column
  for (size_t i = 0; i < header.size(); ++i) {
    const auto& col = header.at(i);
    if (vid_map.contains(col.enum_val)) {
      result[i] = vid_map.at(col.enum_val);
    } else if (var_feat_map.contains(col.enum_val)) {
      result[i] = var_feat_map.at(col.enum_val);
    } else if (ref_feat_map.contains(col.enum_val)) {
      result[i] = ref_feat_map.at(col.enum_val);
    } else {
      throw error::Error("Cannot serialize feature column {} for variant {}", GetFeatureName(col), vid.ToString());
    }
  }
  return result;
}

vec<std::string> VariantInfoSerializer::SerializeTumorNormalBamFeatureRow(const vec<FeatureColumn>& header,
                                                                          const VariantId& vid,
                                                                          const TumorNormalBamFeatureTuple& feat_row,
                                                                          const bool has_tumor_feat,
                                                                          const bool has_normal_feat) {
  using enum SampleContext;
  vec<std::string> result(header.size());
  // Serialize the variant ID and BAM features into maps for easy lookup
  const auto vid_map = SerializeVariantId(vid);
  const auto var_feat_map = SerializeVariantFeature(feat_row.var_feat);
  const auto ref_feat_map = SerializeReferenceFeature(feat_row.ref_feat);
  const auto tumor_var_feat_map =
      has_tumor_feat ? SerializeVariantFeature(feat_row.tumor_var_feat) : std::map<UnifiedFeatureCols, std::string>();
  const auto tumor_ref_feat_map =
      has_tumor_feat ? SerializeReferenceFeature(feat_row.tumor_ref_feat) : std::map<UnifiedFeatureCols, std::string>();
  const auto normal_var_feat_map =
      has_normal_feat ? SerializeVariantFeature(feat_row.normal_var_feat) : std::map<UnifiedFeatureCols, std::string>();
  const auto normal_ref_feat_map = has_normal_feat ? SerializeReferenceFeature(feat_row.normal_ref_feat)
                                                   : std::map<UnifiedFeatureCols, std::string>();

  // Fill in the result based on the feature column
  for (size_t i = 0; i < header.size(); ++i) {
    const auto& col = header.at(i);
    if (vid_map.contains(col.enum_val)) {
      result[i] = vid_map.at(col.enum_val);
    } else {
      // select the maps to use based on the column prefix
      const std::map<UnifiedFeatureCols, std::string>* var_map = &var_feat_map;
      const std::map<UnifiedFeatureCols, std::string>* ref_map = &ref_feat_map;
      if (col.sample_context == kTumor) {
        var_map = &tumor_var_feat_map;
        ref_map = &tumor_ref_feat_map;
      } else if (col.sample_context == kNormal) {
        ref_map = &normal_ref_feat_map;
        var_map = &normal_var_feat_map;
      } else {
        // default sample context; use the pre-initialized maps
      }
      // lookup in the selected maps
      if (var_map->contains(col.enum_val)) {
        result[i] = var_map->at(col.enum_val);
      } else if (ref_map->contains(col.enum_val)) {
        result[i] = ref_map->at(col.enum_val);
      } else {
        throw error::Error("Cannot serialize feature column {} for variant {}", GetFeatureName(col), vid.ToString());
      }
    }
  }
  return result;
}

BamFeatureTuple VariantInfoSerializer::DeserializeBamFeatureRow(const vec<FeatureColumn>& header,
                                                                const vec<std::string>& fields) {
  BamFeatureTuple result;
  for (size_t i = 0; i < header.size(); ++i) {
    const auto& column = header.at(i);
    if (column.sample_context == SampleContext::kNone) {
      const auto& value = fields.at(i);
      if (!UpdateVariantBamFeature(column, value, result.var_feat)) {
        UpdateReferenceBamFeature(column, value, result.ref_feat);
      }
    }
  }
  return result;
}

TumorNormalBamFeatureTuple VariantInfoSerializer::DeserializeTumorNormalBamFeatureRow(const vec<FeatureColumn>& header,
                                                                                      const vec<std::string>& fields) {
  using enum SampleContext;
  TumorNormalBamFeatureTuple result;
  for (size_t i = 0; i < header.size(); ++i) {
    const auto& column = header.at(i);
    const auto& value = fields.at(i);
    auto select_var_feat = [&column, &result]() -> UnifiedVariantFeature& {
      if (column.sample_context == kTumor) {
        return result.tumor_var_feat;
      }
      if (column.sample_context == kNormal) {
        return result.normal_var_feat;
      }
      return result.var_feat;
    };
    auto select_ref_feat = [&column, &result]() -> UnifiedReferenceFeature& {
      if (column.sample_context == kTumor) {
        return result.tumor_ref_feat;
      }
      if (column.sample_context == kNormal) {
        return result.normal_ref_feat;
      }
      return result.ref_feat;
    };
    if (!UpdateVariantBamFeature(column, value, select_var_feat())) {
      UpdateReferenceBamFeature(column, value, select_ref_feat());
    }
  }
  return result;
}

VcfFeature VariantInfoSerializer::DeserializeVcfFeatureRow(const vec<FeatureColumn>& header,
                                                           const vec<std::string>& fields) {
  using enum UnifiedFeatureCols;
  VcfFeature feat{};
  for (size_t i = 0; i < header.size(); ++i) {
    const auto col = header[i].enum_val;
    const auto& value = fields[i];
    switch (col) {
      case kChrom:
        feat.chrom = value;
        break;
      case kPos:
        // convert from 1-based to 0-based
        feat.pos = util::ParseU64(value) - 1;
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
        feat.mpos = util::ParseU64(value);
        break;
      case kVcfMmqRef:
        feat.mmq_ref = util::ParseU8(value);
        break;
      case kVcfMmqAlt:
        feat.mmq_alt = util::ParseU8(value);
        break;
      case kVcfMbqRef:
        feat.mbq_ref = util::ParseU8(value);
        break;
      case kVcfMbqAlt:
        feat.mbq_alt = util::ParseU8(value);
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
        feat.uniq_3mers = util::ParseU32(value);
        break;
      case kVcfUnique4mers:
        feat.uniq_4mers = util::ParseU32(value);
        break;
      case kVcfUnique5mers:
        feat.uniq_5mers = util::ParseU32(value);
        break;
      case kVcfUnique6mers:
        feat.uniq_6mers = util::ParseU32(value);
        break;
      case kVcfHomopolymer:
        feat.homopolymer = util::ParseU32(value);
        break;
      case kVcfDirepeat:
        feat.direpeat = util::ParseU32(value);
        break;
      case kVcfTrirepeat:
        feat.trirepeat = util::ParseU32(value);
        break;
      case kVcfQuadrepeat:
        feat.quadrepeat = util::ParseU32(value);
        break;
      case kVcfVariantDensity:
        feat.variant_density = util::ParseU32(value);
        break;
      case kVcfRefAd:
        feat.ref_ad = util::ParseU32(value);
        break;
      case kVcfAltAd:
        feat.alt_ad = util::ParseU32(value);
        break;
      case kVcfAltAd2:
        feat.alt_ad2 = util::ParseU32(value);
        break;
      case kVcfRefAdAf:
        feat.ref_ad_af = std::stof(value);
        break;
      case kVcfAltAdAf:
        feat.alt_ad_af = std::stof(value);
        break;
      case kVcfAltAd2Af:
        feat.alt_ad2_af = std::stof(value);
        break;
      case kVcfTumorAltAd:
        feat.tumor_alt_ad = util::ParseU32(value);
        break;
      case kVcfNormalAltAd:
        feat.normal_alt_ad = util::ParseU32(value);
        break;
      case kVcfTumorAf:
        feat.tumor_af = std::stof(value);
        break;
      case kVcfNormalAf:
        feat.normal_af = std::stof(value);
        break;
      case kVcfTnAfRatio:
        feat.tn_af_ratio = std::stof(value);
        break;
      case kVcfTumorDp:
        feat.tumor_dp = util::ParseU32(value);
        break;
      case kVcfNormalDp:
        feat.normal_dp = util::ParseU32(value);
        break;
      case kVcfPopAf:
        feat.popaf = std::stof(value);
        break;
      case kVcfPonOcc:
        feat.pon_occ = static_cast<u32>(util::ParseU64(value));
        break;
      case kVcfPonOccRatio:
        feat.pon_occ_ratio = std::stof(value);
        break;
      case kVcfPonMeanDuplexAf:
        feat.pon_mean_duplex_af = std::stof(value);
        break;
      case kVcfPonMinDuplexAf:
        feat.pon_min_duplex_af = std::stof(value);
        break;
      case kVcfPonMaxDuplexAf:
        feat.pon_max_duplex_af = std::stof(value);
        break;
      case kVcfGenotype:
        feat.genotype = StringToGenotype(value);
        break;
      case kVcfQual:
        feat.qual = std::stof(value);
        break;
      case kVcfGq:
        feat.gq = util::ParseU32(value);
        break;
      case kVcfHapcomp:
        feat.hapcomp = util::ParseU32(value);
        break;
      case kVcfHapdom:
        feat.hapdom = std::stof(value);
        break;
      case kVcfMqRankSum:
        feat.mq_rank_sum = std::stof(value);
        break;
      case kVcfSor:
        feat.sor = std::stof(value);
        break;
      case kVcfRu:
        feat.ru = DotToEmptyString(value);
        break;
      case kVcfRpaRef:
        feat.rpa_ref = ParseRpaValue(value);
        break;
      case kVcfRpaAlt:
        feat.rpa_alt = ParseRpaValue(value);
        break;
      case kVcfStr:
        feat.str = (util::ParseU64(value) != 0U);
        break;
      case kVcfAtInterest:
        feat.at_interest = (util::ParseU64(value) != 0U);
        break;
      case kVcfInCatalog:
        feat.in_catalog = (util::ParseU64(value) != 0U);
        break;
      case kVcfIsInRepeat:
        feat.is_in_repeat = (util::ParseU64(value) != 0U);
        break;
      case kVcfRepeatLen:
        feat.repeat_len = util::ParseU32(value);
        break;
      case kVcfRepeatUnitLen:
        feat.repeat_unit_len = util::ParseU32(value);
        break;
      case kVcfIsUnitSizeMatch:
        feat.is_unit_size_match = (util::ParseU64(value) != 0U);
        break;
      case kVcfGcContent:
        feat.gc_content = std::stof(value);
        break;
      default:
        break;
    }
  }
  return feat;
}

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
  result[kVcfTnAfRatio] = std::to_string(feature.tn_af_ratio);
  result[kVcfTumorDp] = std::to_string(feature.tumor_dp);
  result[kVcfNormalDp] = std::to_string(feature.normal_dp);
  result[kVcfPopAf] = std::to_string(feature.popaf);
  result[kVcfPonOcc] = std::to_string(feature.pon_occ);
  result[kVcfPonOccRatio] = std::to_string(feature.pon_occ_ratio);
  result[kVcfPonMeanDuplexAf] = std::to_string(feature.pon_mean_duplex_af);
  result[kVcfPonMinDuplexAf] = std::to_string(feature.pon_min_duplex_af);
  result[kVcfPonMaxDuplexAf] = std::to_string(feature.pon_max_duplex_af);
  result[kVcfGenotype] = GenotypeToString(feature.genotype);
  result[kVcfQual] = std::to_string(feature.qual);
  result[kVcfGq] = std::to_string(feature.gq);
  result[kVcfHapcomp] = std::to_string(feature.hapcomp);
  result[kVcfHapdom] = std::to_string(feature.hapdom);
  result[kVcfMqRankSum] = std::to_string(feature.mq_rank_sum);
  result[kVcfSor] = std::to_string(feature.sor);
  result[kVcfRu] = EmptyStringToDot(feature.ru);
  result[kVcfRpaRef] = std::to_string(feature.rpa_ref);
  result[kVcfRpaAlt] = std::to_string(feature.rpa_alt);
  result[kVcfStr] = feature.str ? "1" : "0";
  result[kVcfAtInterest] = feature.at_interest ? "1" : "0";
  result[kVcfIsInRepeat] = feature.is_in_repeat ? "1" : "0";
  result[kVcfRepeatLen] = std::to_string(feature.repeat_len);
  result[kVcfRepeatUnitLen] = std::to_string(feature.repeat_unit_len);
  result[kVcfIsUnitSizeMatch] = feature.is_unit_size_match ? "1" : "0";
  result[kVcfInCatalog] = feature.in_catalog ? "1" : "0";
  result[kVcfGcContent] = std::to_string(feature.gc_content);
  return result;
}

f64 VariantInfoSerializer::NumericalizeFeature(const UnifiedFeatureCols col, const VariantId& vid) {
  using enum UnifiedFeatureCols;
  switch (col) {
    case kChrom:
      return static_cast<f64>(std::hash<std::string>{}(vid.chrom));
    case kPos:
      // convert from 0-based to 1-based
      return static_cast<f64>(vid.pos + 1);
    case kRef:
      return SeqToDouble(vid.ref);
    case kAlt:
      return SeqToDouble(vid.alt);
    case kSubTypeIndex:
      return SubstIndex(vid.ref, vid.alt);
    case kVariantType:
      return static_cast<f64>(vid.type);
    case kAltLen:
      return static_cast<f64>(vid.alt.size()) - static_cast<f64>(vid.ref.size());
    default:
      // not an attribute of VariantId
      return NAN;
  }
}

f64 VariantInfoSerializer::NumericalizeFeature(const UnifiedFeatureCols col, const UnifiedVariantFeature& feat) {
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
      return static_cast<f64>(feat.distance_min);
    case kDistanceMax:
      return static_cast<f64>(feat.distance_max);
    case kDistanceSum:
      return static_cast<f64>(feat.distance_sum);
    case kDistanceSumLowbq:
      return static_cast<f64>(feat.distance_sum_lowbq);
    case kDistanceSumSimplex:
      return static_cast<f64>(feat.distance_sum_simplex);
    case kDistanceMean:
      return feat.distance_mean;
    case kDistanceMeanLowbq:
      return feat.distance_mean_lowbq;
    case kDistanceMeanSimplex:
      return feat.distance_mean_simplex;
    case kDuplexLowbq:
      return feat.duplex_lowbq;
    case kDuplexConcordant:
      return feat.duplex_concordant;
    case kDuplexSimplex:
      return feat.duplex_simplex;
    case kDuplexDiscordant:
      return feat.duplex_discordant;
    case kSimplex:
      return feat.simplex;
    case kDuplexAF:
      return feat.duplex_af;
    case kDuplexOther:
      return feat.duplex_other;
    case kDuplexOtherAF:
      return feat.duplex_other_af;
    case kDuplexGQ:
      return feat.duplex_gq;
    case kAfqGT01:
      return feat.afq_gt01;
    case kAfqGT11:
      return feat.afq_gt11;
    case kAfqGT12:
      return feat.afq_gt12;
    case kAfqGT:
      return GenotypeToInt(feat.afq_gt);
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
    case kBamTnAfRatio:
      return feat.tn_af_ratio;
    case kSupportReverse:
      return feat.support_reverse;
    case kAlignmentBias:
      return feat.alignmentbias;
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
      return static_cast<f64>(util::hash::HashRange(feat.filter_status.begin(), feat.filter_status.end()));
    case kADT:
      return feat.adt;
    case kADTL:
      return feat.adtl;
    case kIndelAf:
      return feat.indel_af;
    default:
      // not an attribute of UnifiedVariantFeature
      return NAN;
  }
}

f64 VariantInfoSerializer::NumericalizeFeature(const UnifiedFeatureCols col, const VcfFeature& feat) {
  using enum UnifiedFeatureCols;
  switch (col) {
    case kChrom:
      return static_cast<f64>(std::hash<std::string>{}(feat.chrom));
    case kPos:
      // convert from 0-based to 1-based
      return static_cast<f64>(feat.pos + 1);
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
      return static_cast<f64>(feat.mpos);
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
      }
      return UnifiedVariantFeature::ContextIndex(feat.pre_2bp_context);
    case kVcfPost2bContext:
      if (IsAnyNotACTG(feat.post_2bp_context)) {
        return UnifiedVariantFeature::ContextIndex("");
      }
      return UnifiedVariantFeature::ContextIndex(feat.post_2bp_context);
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
    case kVcfTnAfRatio:
      return feat.tn_af_ratio;
    case kVcfTumorDp:
      return feat.tumor_dp;
    case kVcfNormalDp:
      return feat.normal_dp;
    case kVcfPopAf:
      return feat.popaf;
    case kVcfPonOcc:
      return feat.pon_occ;
    case kVcfPonOccRatio:
      return feat.pon_occ_ratio;
    case kVcfPonMeanDuplexAf:
      return feat.pon_mean_duplex_af;
    case kVcfPonMinDuplexAf:
      return feat.pon_min_duplex_af;
    case kVcfPonMaxDuplexAf:
      return feat.pon_max_duplex_af;
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
    case kVcfMqRankSum:
      return feat.mq_rank_sum;
    case kVcfSor:
      return feat.sor;
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
    case kVcfIsInRepeat:
      return feat.is_in_repeat ? 1 : 0;
    case kVcfRepeatLen:
      return feat.repeat_len;
    case kVcfRepeatUnitLen:
      return feat.repeat_unit_len;
    case kVcfIsUnitSizeMatch:
      return feat.is_unit_size_match ? 1 : 0;
    case kVcfInCatalog:
      return feat.in_catalog ? 1 : 0;
    case kVcfGcContent:
      return feat.gc_content;
    default:
      // not an attribute of VcfFeature
      return NAN;
  }
}

f64 VariantInfoSerializer::NumericalizeFeature(const UnifiedFeatureCols col, const UnifiedReferenceFeature& feat) {
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
    case kDel1bpAF:
      return feat.del_1bp_af;
    case kDel2bpAF:
      return feat.del_2bp_af;
    case kDel3bpAF:
      return feat.del_3bp_af;
    case kDelGt3bpAF:
      return feat.del_gt3bp_af;
    case kIns1bpAF:
      return feat.ins_1bp_af;
    case kIns2bpAF:
      return feat.ins_2bp_af;
    case kIns3bpAF:
      return feat.ins_3bp_af;
    case kInsGt3bpAF:
      return feat.ins_gt3bp_af;
    case kDel1bp2bpAFLogRatio:
      return feat.del_1bp_2bp_af_logratio;
    case kDel2bp3bpAFLogRatio:
      return feat.del_2bp_3bp_af_logratio;
    case kDel3bpGt3bpAFLogRatio:
      return feat.del_3bp_gt3bp_af_logratio;
    case kIns1bp2bpAFLogRatio:
      return feat.ins_1bp_2bp_af_logratio;
    case kIns2bp3bpAFLogRatio:
      return feat.ins_2bp_3bp_af_logratio;
    case kIns3bpGt3bpAFLogRatio:
      return feat.ins_3bp_gt3bp_af_logratio;
    case kDel1bpIns1bpAFLogRatio:
      return feat.del_1bp_ins_1bp_af_logratio;
    case kSnvAF:
      return feat.snv_af;
    case kRefBaseqLT20Count:
      return feat.baseq_lt20_count;
    case kRefBaseqLT20Ratio:
      return feat.baseq_lt20_ratio;
    case kRefDistanceMin:
      return static_cast<f64>(feat.distance_min);
    case kRefDistanceMax:
      return static_cast<f64>(feat.distance_max);
    case kRefDistanceSum:
      return static_cast<f64>(feat.distance_sum);
    case kRefDistanceSumLowbq:
      return static_cast<f64>(feat.distance_sum_lowbq);
    case kRefDistanceSumSimplex:
      return static_cast<f64>(feat.distance_sum_simplex);
    case kRefDistanceMean:
      return feat.distance_mean;
    case kRefDistanceMeanLowbq:
      return feat.distance_mean_lowbq;
    case kRefDistanceMeanSimplex:
      return feat.distance_mean_simplex;
    case kRefDuplexLowbq:
      return feat.duplex_lowbq;
    case kRefDuplexConcordant:
      return feat.duplex_concordant;
    case kRefDuplexSimplex:
      return feat.duplex_simplex;
    case kRefDuplexDiscordant:
      return feat.duplex_discordant;
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
    default:
      // not an attribute of UnifiedReferenceFeature
      return NAN;
  }
}

f64 VariantInfoSerializer::NumericalizeFeature(const UnifiedFeatureCols col,
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

namespace {
/**
 * @brief Parse the header line of a TSV features file and return a vector of UnifiedFeatureCols
 * @param input Input file stream
 * @param tsv_file Path to the features file (for error messages)
 * @return Vector of header column names
 */
vec<std::string> ParseTsvHeader(std::ifstream& input, const fs::path& tsv_file) {
  std::string line;
  bool header_found = false;
  while (std::getline(input, line)) {
    if (line.empty() || line.starts_with(io::kTsvCommentLinePrefix)) {
      continue;
    }
    // first non-empty, non-comment line
    header_found = true;
    break;
  }
  if (!header_found) {
    throw error::Error("Cannot find header from TSV file: {}", tsv_file);
  }
  // trim leading and trailing whitespace from the header line
  string::Trim(line);
  return string::Split(line, "\t");
}

/**
 * @brief Convert a vector of header column names to a vector of UnifiedFeatureCols
 * @param header_names Vector of header column names
 * @param infer_sample_context Whether to infer sample context from column name prefixes (e.g. "tumor_", "normal_")
 * @return Vector of UnifiedFeatureCols
 */
vec<FeatureColumn> GetFeatureColumns(const vec<std::string>& header_names, const bool infer_sample_context) {
  vec<FeatureColumn> header_cols;
  header_cols.reserve(header_names.size());
  for (const auto& f : header_names) {
    header_cols.emplace_back(GetFeatureColumn(f, infer_sample_context));
  }
  return header_cols;
}
}  // namespace

vec<FeatureColumn> VariantInfoSerializer::ParseBamFeaturesFileHeader(std::ifstream& input,
                                                                     const fs::path& features_file,
                                                                     const bool infer_sample_context) {
  const vec<std::string>& header_names = ParseTsvHeader(input, features_file);
  const auto& unsupported_fields = FindUnsupportedBamFeatureNames(header_names, infer_sample_context);
  if (!unsupported_fields.empty()) {
    throw error::Error("Unsupported BAM feature names: {}", string::Join(unsupported_fields, ", "));
  }
  return GetFeatureColumns(header_names, infer_sample_context);
}

/**
 * @brief Parse the header of a VCF features file and return a vector of UnifiedFeatureCols
 * @param input Input file stream
 * @param features_file Path to VCF features file
 * @param infer_sample_context Whether to infer sample context from column name prefixes (e.g. "tumor_", "normal_")
 * @return Vector of UnifiedFeatureCols
 */
vec<FeatureColumn> VariantInfoSerializer::ParseVcfFeaturesFileHeader(std::ifstream& input,
                                                                     const fs::path& features_file,
                                                                     const bool infer_sample_context) {
  const vec<std::string>& header_names = ParseTsvHeader(input, features_file);
  // VCF feature names already have hard-coded sample context prefixes.
  // So, we check for unsupported feature names without inferring sample context, and throw an error if any unsupported
  // feature names are found.
  const auto& unsupported_fields = FindUnsupportedVcfFeatureNames(header_names);
  if (!unsupported_fields.empty()) {
    throw error::Error("Unsupported VCF feature names: {}", string::Join(unsupported_fields, ", "));
  }
  return GetFeatureColumns(header_names, infer_sample_context);
}

void VariantInfoSerializer::ValidateFeatureFileHeaderForNormalization(const fs::path& features_file,
                                                                      const bool has_sample_context) {
  using enum UnifiedFeatureCols;
  std::ifstream input(features_file);
  const auto columns = ParseVcfFeaturesFileHeader(input, features_file, false);
  bool tumor_dp_found = false;
  bool normal_dp_found = false;
  for (const auto& [enum_val, sample_context] : columns) {
    if (enum_val == kVcfTumorDp) {
      tumor_dp_found = true;
    } else if (enum_val == kVcfNormalDp) {
      normal_dp_found = true;
    } else {
      // not a DP column; skip
    }
    if (tumor_dp_found && normal_dp_found) {
      break;
    }
  }
  if (has_sample_context && !tumor_dp_found) {
    throw error::Error(
        "VCF feature '{}' is required for feature normalization in the tumor sample, but it is not found in features "
        "file: {}",
        kNameTumorDP,
        features_file);
  }
  if (!normal_dp_found) {
    throw error::Error(
        "VCF feature '{}' is required for feature normalization, but it is not found in features file: {}",
        kNameNormalDP,
        features_file);
  }
}

vec<std::string> VariantInfoSerializer::ParseNextTsvRow(std::ifstream& input, const size_t num_cols) {
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line.starts_with(io::kTsvCommentLinePrefix)) {
      continue;
    }
    auto field_values = string::Split(line, "\t");
    if (field_values.size() != num_cols) {
      throw error::Error("Number of TSV fields ({}) does not match TSV header ({})", field_values.size(), num_cols);
    }
    // trim leading and trailing whitespace from each field
    for (auto& field : field_values) {
      string::Trim(field);
    }
    // valid line read
    return field_values;
  }
  // end of file reached
  return {};
}

VarIdToVcfFeatures VariantInfoSerializer::LoadVcfFeatures(const fs::path& features_file,
                                                          const bool infer_sample_context) {
  auto empty_set = std::unordered_set<VariantId>();
  return LoadVcfFeatures<std::unordered_set<VariantId>>(features_file, empty_set, infer_sample_context);
}

auto VariantInfoSerializer::BamFeatureTupleGenerator(const fs::path& features_file)
    -> std::function<std::optional<std::pair<VariantId, BamFeatureTuple>>()> {
  auto input = std::make_shared<std::ifstream>(features_file);
  if (!input->is_open()) {
    throw error::Error("Failed to open BAM feature file: " + features_file.string());
  }
  // Parse header columns
  const vec<FeatureColumn> header = ParseBamFeaturesFileHeader(*input, features_file, false);
  const auto num_cols = header.size();

  // Lambda state for generator
  return [input_stream = std::move(input),
          header,
          num_cols]() mutable -> std::optional<std::pair<VariantId, BamFeatureTuple>> {
    // Read and return the next row as a (VariantId, BamFeatureTuple) pair, or nullopt at EOF
    const vec<std::string> fields = ParseNextTsvRow(*input_stream, num_cols);
    if (fields.empty()) {
      return std::nullopt;
    }
    const VariantId vid = DeserializeVariantId(header, fields);
    const BamFeatureTuple tuple = DeserializeBamFeatureRow(header, fields);
    return std::make_pair(vid, tuple);
  };
}

auto VariantInfoSerializer::TumorNormalBamFeatureTupleGenerator(const fs::path& features_file)
    -> std::function<std::optional<std::pair<VariantId, TumorNormalBamFeatureTuple>>()> {
  auto input = std::make_shared<std::ifstream>(features_file);
  if (!input->is_open()) {
    throw error::Error("Failed to open BAM feature file: " + features_file.string());
  }
  // Parse header columns
  const vec<FeatureColumn> header = ParseBamFeaturesFileHeader(*input, features_file, true);
  const auto num_cols = header.size();

  // Lambda state for generator
  return [input_stream = std::move(input),
          header,
          num_cols]() mutable -> std::optional<std::pair<VariantId, TumorNormalBamFeatureTuple>> {
    // Read and return the next row as a (VariantId, TumorNormalBamFeatureTuple) pair, or nullopt at EOF
    const vec<std::string> fields = ParseNextTsvRow(*input_stream, num_cols);
    if (fields.empty()) {
      return std::nullopt;
    }
    const VariantId vid = DeserializeVariantId(header, fields);
    const TumorNormalBamFeatureTuple tuple = DeserializeTumorNormalBamFeatureRow(header, fields);
    return std::make_pair(vid, tuple);
  };
}

}  // namespace xoos::svc
