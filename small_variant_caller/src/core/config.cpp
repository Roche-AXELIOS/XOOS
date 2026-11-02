#include "config.h"

#include <algorithm>
#include <fstream>

#include <fmt/format.h>

#include <xoos/enum/enum-util.h>
#include <xoos/error/error.h>
#include <xoos/util/string-functions.h>

#include "column-names.h"
#include "util/model-resolver.h"
#include "variant-info.h"

namespace xoos::svc {

namespace {
constexpr u32 kGermlineTaggingModelIterations{1500};
// Training iterations for the large germline WGS model family.
constexpr u32 kGermlineWgsLargeSnvIterations{1000};
constexpr u32 kGermlineWgsLargeIndelIterations{8000};
// Training iterations for the small germline WGS model family.
constexpr u32 kGermlineWgsSmallSnvIterations{1500};
constexpr u32 kGermlineWgsSmallIndelIterations{1500};
constexpr auto* kGermlineWgsLargeSnvModelLgbmParams =
    "objective=multiclass boosting=gbdt metric=multi_logloss seed=1238845 learning_rate=0.0286803679209628 "
    "bagging_fraction=0.683997122877158 bagging_freq=1 feature_fraction=0.753244391002605 num_leaves=58 "
    "min_data_in_leaf=435 num_classes=4 ";
constexpr auto* kGermlineWgsLargeIndelModelLgbmParams =
    "objective=multiclass boosting=gbdt metric=multi_logloss seed=1238845 learning_rate=0.0420473270669401 "
    "bagging_fraction=0.88062343020432 bagging_freq=1 num_leaves=59 feature_fraction=0.902813236067906 "
    "min_data_in_leaf=1325 num_classes=4 ";
constexpr auto* kGermlineWgsSmallSnvModelLgbmParams =
    "objective=multiclass boosting=gbdt metric=multi_logloss seed=1238845 learning_rate=0.01 "
    "bagging_fraction=0.9 bagging_freq=1 feature_fraction=0.5 num_leaves=64 num_classes=4 ";
constexpr auto* kGermlineWgsSmallIndelModelLgbmParams =
    "objective=multiclass boosting=gbdt metric=multi_logloss seed=1238845 learning_rate=0.02 "
    "bagging_fraction=0.9 bagging_freq=1 num_leaves=64 feature_fraction=0.5 num_classes=4 ";
// Minimum homopolymer length used by the non-germline workflows.
constexpr u32 kNonGermlineMinHomopolymerLength{7};

// Tumor-only TE workflow tuning parameters.
constexpr u8 kTumorOnlyTeMinMapq{9};
constexpr u32 kTumorOnlyTeMaxVariantsPerRead{10};
constexpr u32 kTumorOnlyTeIterations{1000};
constexpr f32 kTumorOnlyTeMinAf{0.01F};

// Phased-variant allele-frequency bounds shared by the somatic workflows.
constexpr f32 kSomaticMinPhasedAf{0.001F};
constexpr f32 kSomaticMaxPhasedAf{0.5F};
// Default ML score threshold and hotspot ML score threshold shared by the somatic workflows.
constexpr f32 kSomaticDefaultMlScore{0.3F};
constexpr f32 kSomaticHotspotMlScore{0.01F};

// Tumor-normal WGS workflow tuning parameters.
constexpr u32 kTumorNormalWgsIterations{6000};
constexpr f32 kTumorNormalWgsSnvMinMlScore{0.03F};
constexpr f32 kTumorNormalWgsIndelMinMlScore{0.008F};
constexpr f32 kTumorNormalWgsFfpeMinTumorAf{0.05F};
constexpr f32 kTumorNormalWgsCellLineMinTumorAf{0.03F};
constexpr u32 kTumorNormalWgsMaxIndelSize{30};
constexpr f32 kTumorNormalWgsMinDpRatio{0.15F};
// Per-sample-type model score thresholds for the tumor-normal WGS workflow.
// Converted from the validated FFPE model QUAL cutoffs: SNV 12.29 and indel 2.46.
constexpr f32 kTumorNormalWgsFfpeSnvMlScore{0.940979891982816F};
constexpr f32 kTumorNormalWgsFfpeIndelMlScore{0.432455394591453F};
constexpr f32 kTumorNormalWgsCellLineSnvMlScore{0.5257580147397554F};
constexpr f32 kTumorNormalWgsCellLineIndelMlScore{0.609159104207598F};

// Tumor-only WGS workflow tuning parameters.
constexpr u32 kTumorOnlyWgsIterations{6000};
constexpr f32 kTumorOnlyWgsSnvMinMlScore{0.001F};
constexpr f32 kTumorOnlyWgsIndelMinMlScore{0.001F};
constexpr f32 kTumorOnlyWgsMinAf{0.05F};
constexpr u32 kTumorOnlyWgsMaxIndelSize{30};
// Per-sample-type score thresholds for the auto-resolved tumor-only WGS BWA models. `cell-line`
// selects the multi-cell-line model (HCC1395 + H2009 + HG008); `ffpe` selects the FFPE model.
// The HCC1395-only model remains available through an explicit `--model`.
constexpr f32 kTumorOnlyWgsCellLineSnvMlScore{0.07F};
constexpr f32 kTumorOnlyWgsCellLineIndelMlScore{0.09F};
constexpr f32 kTumorOnlyWgsCellLineMinAf{0.03F};
constexpr f32 kTumorOnlyWgsFfpeSnvMlScore{0.929F};
constexpr f32 kTumorOnlyWgsFfpeIndelMlScore{0.818F};
constexpr f32 kTumorOnlyWgsFfpeMinAf{0.05F};

constexpr f32 kGermlineMitochondriaMinDist{5.0F};
constexpr f32 kGermlineMitochondriaMinAf{0.015F};
constexpr u32 kGermlineMitochondriaMinDp{400};
constexpr f32 kGermlineMitochondriaMinHomAltAf{0.95F};

/**
 * @brief Convert a germline model type to the JSON map key used in profile configuration.
 * @param model_type Model type enum value.
 * @return Profile JSON key for the model type.
 */
std::string GermlineModelTypeConfigKey(const GermlineModelType model_type) {
  return enum_util::FormatEnumName(model_type);
}

/**
 * @brief Build the current large germline WGS training defaults.
 * @return Training defaults for the large bundled germline WGS model family.
 */
GermlineModelTrainingDefaults LargeGermlineModelTrainingDefaults() {
  return {.snv_iterations = kGermlineWgsLargeSnvIterations,
          .indel_iterations = kGermlineWgsLargeIndelIterations,
          .snv_model_lgbm_params = kGermlineWgsLargeSnvModelLgbmParams,
          .indel_model_lgbm_params = kGermlineWgsLargeIndelModelLgbmParams};
}

/**
 * @brief Build the small germline WGS training defaults from the legacy base-branch profile.
 * @return Training defaults for the small bundled germline WGS model family.
 */
GermlineModelTrainingDefaults SmallGermlineModelTrainingDefaults() {
  return {.snv_iterations = kGermlineWgsSmallSnvIterations,
          .indel_iterations = kGermlineWgsSmallIndelIterations,
          .snv_model_lgbm_params = kGermlineWgsSmallSnvModelLgbmParams,
          .indel_model_lgbm_params = kGermlineWgsSmallIndelModelLgbmParams};
}

/**
 * @brief Populate germline WGS model-type training defaults from built-in values.
 * @details Both germline library profiles share model-family training hyperparameters today. Keeping this setup in one
 * helper prevents the simplex and duplex mirrors from drifting as the model families evolve.
 * @param config Config profile whose defaults are updated.
 */
void PopulateGermlineModelTrainingDefaults(SVCConfig& config) {
  config.germline_model_type_training_defaults = {
      {GermlineModelTypeConfigKey(GermlineModelType::kLarge), LargeGermlineModelTrainingDefaults()},
      {GermlineModelTypeConfigKey(GermlineModelType::kSmall), SmallGermlineModelTrainingDefaults()}};
}

/**
 * @brief Normalize germline WGS training defaults after JSON parsing.
 * @details New profile JSON stores model-family training settings in `germline_model_type_training_defaults`, so the
 * top-level active training fields are derived from the `large` entry. Older custom JSON files without the nested map
 * keep working by treating their top-level training fields as the large-family defaults. Small-model training remains
 * an explicit config requirement so stale custom configs fail clearly instead of silently using large settings.
 * @param config Config parsed from JSON.
 */
void NormalizeGermlineModelTrainingDefaults(SVCConfig& config) {
  if (config.workflow != Workflow::kGermlineWgs) {
    return;
  }
  const auto large_key = GermlineModelTypeConfigKey(GermlineModelType::kLarge);
  const auto large_defaults_it = config.germline_model_type_training_defaults.find(large_key);
  if (large_defaults_it == config.germline_model_type_training_defaults.end()) {
    config.germline_model_type_training_defaults.try_emplace(
        large_key,
        GermlineModelTrainingDefaults{.snv_iterations = config.snv_iterations,
                                      .indel_iterations = config.indel_iterations,
                                      .snv_model_lgbm_params = config.snv_model_lgbm_params,
                                      .indel_model_lgbm_params = config.indel_model_lgbm_params});
    return;
  }

  const auto& large_defaults = large_defaults_it->second;
  config.snv_iterations = large_defaults.snv_iterations;
  config.indel_iterations = large_defaults.indel_iterations;
  config.snv_model_lgbm_params = large_defaults.snv_model_lgbm_params;
  config.indel_model_lgbm_params = large_defaults.indel_model_lgbm_params;
}
}  // namespace

std::string AddTumorPrefix(std::string_view feature_name) {
  return fmt::format("{}{}", kTumorPrefix, feature_name);
}

std::string AddNormalPrefix(std::string_view feature_name) {
  return fmt::format("{}{}", kNormalPrefix, feature_name);
}

vec<std::string> DefaultBamFeatures() {
  return ToStringVec(
      kNameChrom, kNamePos, kNameRef, kNameAlt, kNameDuplex, kNameDuplexLowbq, kNameDuplexDP, kNameDuplexAF);
}

vec<std::string> GermlineMitochondriaBamFeatures() {
  return ToStringVec(kNameChrom, kNamePos, kNameRef, kNameAlt, kNameDuplex, kNameDuplexDP, kNameDuplexAF);
}

vec<std::string> GermlineMitochondriaVcfFeatures() {
  return ToStringVec(kNameChrom, kNamePos, kNameRef, kNameAlt, kNameAltAD, kNameNormalDP, kNameAltAdAF);
}

vec<std::string> DefaultVcfFeatures() {
  return ToStringVec(kNameChrom,
                     kNamePos,
                     kNameRef,
                     kNameAlt,
                     kNameVariantType,
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
                     kNameNormalDP,
                     kNamePopAF);
}

vec<std::string> TumorOnlyTeBamFeatures() {
  return ToStringVec(kNameChrom,
                     kNamePos,
                     kNameRef,
                     kNameAlt,
                     kNameMapqSum,
                     kNameBaseqSum,
                     kNameNonDuplex,
                     kNameDuplex,
                     kNameFamilysizeSum,
                     kNameDistanceSum,
                     kNamePlusOnly,
                     kNameMinusOnly,
                     kNameRefSupport,
                     kNameSupport,
                     kNameContext,
                     kNameWeightedScore,
                     kNameBaseqMean,
                     kNameDistanceMean,
                     kNameMapqMean,
                     kNameFamilysizeMean,
                     kNameStrandBias,
                     kNameSubtypeIndex,
                     kNameContextIndex);
}

vec<std::string> TumorOnlyTeScoringNames() {
  return ToStringVec(kNameWeightedScore,
                     kNameBaseqMean,
                     kNameDistanceMean,
                     kNameMapqMean,
                     kNameDuplex,
                     kNameSubtypeIndex,
                     kNameContextIndex,
                     kNameStrandBias,
                     kNamePlusOnly,
                     kNameMinusOnly,
                     kNameFamilysizeSum);
}

vec<std::string> TumorOnlyTeCategoricalNames() {
  return ToStringVec(kNameSubtypeIndex, kNameContextIndex);
}

vec<std::string> GermlineWgsDuplexIndelScoringNames() {
  return ToStringVec(kNameRef,
                     kNameAlt,
                     kNameDuplex,
                     kNameRefNonHpSupport,
                     kNameAltLen,
                     kNameRefSupport,
                     kNameDistanceMin,
                     kNameDistanceMax,
                     kNameDistanceMean,
                     kNameDistanceSum,
                     kNameMapqAF,
                     kNameRefMapqAF,
                     kNameWeightedDepth,
                     kNamePopAF,
                     kNameIndelAF,
                     kNamePre2BpContext,
                     kNamePost2BpContext,
                     kNameHomopolymer,
                     kNameDirepeat,
                     kNameTrirepeat,
                     kNameQuadrepeat,
                     kNameHapcomp,
                     kNameHapdom,
                     kNameRu,
                     kNameRpaRef,
                     kNameRpaAlt,
                     kNameStr,
                     kNameUniq3mers,
                     kNameUniq4mers,
                     kNameUniq5mers,
                     kNameUniq6mers,
                     kNameMapqMin,
                     kNameMapqMean,
                     kNameMapqLT60Ratio,
                     kNameMapqLT30Ratio,
                     kNameRefNonHpMapqLT60Ratio,
                     kNameADTL,
                     kNameNumAlt,
                     kNameNormalDP,
                     kNameRefAD,
                     kNameAltAD,
                     kNameAltAD2,
                     kNameRefAdAF,
                     kNameAltAdAF,
                     kNameAltAd2AF,
                     kNameQual,
                     kNameGq,
                     kNameGenotype,
                     kNameDuplexLowbq,
                     kNameDuplexAF,
                     kNameMapqSumLowbq,
                     kNameMapqMeanLowbq,
                     kNameDistanceSumLowbq,
                     kNameDistanceMeanLowbq,
                     kNameRefDuplexLowbq,
                     kNameRefDuplexAF,
                     kNameRefMapqSumLowbq,
                     kNameRefMapqMeanLowbq,
                     kNameRefDistanceSumLowbq,
                     kNameRefDistanceMeanLowbq,
                     kNameAtInterest);
}

vec<std::string> GermlineWgsDuplexSnvCategoricalScoringNames() {
  return ToStringVec(kNameRef, kNameAlt, kNameGenotype, kNameAtInterest);
}

vec<std::string> GermlineWgsDuplexIndelCategoricalScoringNames() {
  return ToStringVec(
      kNameRef, kNameAlt, kNamePre2BpContext, kNamePost2BpContext, kNameGenotype, kNameRu, kNameStr, kNameAtInterest);
}

vec<std::string> GermlineWgsDuplexBamFeatures() {
  return ToStringVec(kNameChrom,
                     kNamePos,
                     kNameRef,
                     kNameAlt,
                     kNameDuplex,
                     kNameRefNonHpSupport,
                     kNameRefSupport,
                     kNameDistanceMin,
                     kNameDistanceMax,
                     kNameDistanceMean,
                     kNameDistanceSum,
                     kNameMapqAF,
                     kNameRefMapqAF,
                     kNameSupport,
                     kNameWeightedDepth,
                     kNameMapqMin,
                     kNameMapqSum,
                     kNameMapqMean,
                     kNameMapqLT60Ratio,
                     kNameMapqLT30Ratio,
                     kNameRefNonHpMapqLT60Ratio,
                     kNameAltLen,
                     kNameRefDistanceMean,
                     kNameRefMapqMean,
                     kNameMapqLT40Ratio,
                     kNameRefMapqLT30Ratio,
                     kNameRefMapqLT40Ratio,
                     kNameDuplexLowbq,
                     kNameDuplexDP,
                     kNameDuplexAF,
                     kNameMapqSumLowbq,
                     kNameMapqMeanLowbq,
                     kNameDistanceSumLowbq,
                     kNameDistanceMeanLowbq,
                     kNameRefDuplexLowbq,
                     kNameRefDuplexAF,
                     kNameRefMapqSumLowbq,
                     kNameRefMapqMeanLowbq,
                     kNameRefDistanceSumLowbq,
                     kNameRefDistanceMeanLowbq,
                     kNameNumAlt,
                     kNameADTL,
                     kNameIndelAF);
}

vec<std::string> GermlineWgsDuplexVcfFeatures() {
  return ToStringVec(kNameChrom,
                     kNamePos,
                     kNameRef,
                     kNameAlt,
                     kNamePopAF,
                     kNameVariantDensity,
                     kNameGenotype,
                     kNameQual,
                     kNameGq,
                     kNameNormalAF,
                     kNameNormalDP,
                     kNamePre2BpContext,
                     kNamePost2BpContext,
                     kNamePost30BpContext,
                     kNameRefAD,
                     kNameAltAD,
                     kNameAltAD2,
                     kNameRefAdAF,
                     kNameAltAdAF,
                     kNameAltAd2AF,
                     kNameUniq3mers,
                     kNameUniq4mers,
                     kNameUniq5mers,
                     kNameUniq6mers,
                     kNameVariantType,
                     kNameHomopolymer,
                     kNameDirepeat,
                     kNameTrirepeat,
                     kNameQuadrepeat,
                     kNameHapcomp,
                     kNameHapdom,
                     kNameRu,
                     kNameRpaRef,
                     kNameRpaAlt,
                     kNameStr,
                     kNameAtInterest);
}

vec<std::string> GermlineWgsDuplexSnvScoringNames() {
  return ToStringVec(kNameRef,
                     kNameAlt,
                     kNameRefSupport,
                     kNameSupport,
                     kNameRefMapqMean,
                     kNameMapqMean,
                     kNameDistanceMean,
                     kNameRefDistanceMean,
                     kNamePopAF,
                     kNameVariantDensity,
                     kNameQual,
                     kNameGq,
                     kNameNormalDP,
                     kNameRefAD,
                     kNameAltAD,
                     kNameAltAD2,
                     kNameRefAdAF,
                     kNameAltAdAF,
                     kNameAltAd2AF,
                     kNameHapcomp,
                     kNameHapdom,
                     kNameGenotype,
                     kNameMapqLT30Ratio,
                     kNameRefMapqLT30Ratio,
                     kNameMapqLT40Ratio,
                     kNameRefMapqLT40Ratio,
                     kNameDuplexLowbq,
                     kNameDuplexAF,
                     kNameMapqAF,
                     kNameMapqSumLowbq,
                     kNameMapqMeanLowbq,
                     kNameDistanceSumLowbq,
                     kNameDistanceMeanLowbq,
                     kNameRefDuplexLowbq,
                     kNameRefDuplexAF,
                     kNameRefMapqSumLowbq,
                     kNameRefMapqMeanLowbq,
                     kNameRefDistanceSumLowbq,
                     kNameRefDistanceMeanLowbq,
                     kNameAtInterest);
}

vec<std::string> GermlineWgsSimplexBamFeatures() {
  return ToStringVec(kNameChrom,
                     kNamePos,
                     kNameRef,
                     kNameAlt,
                     kNameRefNonHpSupport,
                     kNameRefSupport,
                     kNameDistanceMin,
                     kNameDistanceMax,
                     kNameDistanceMean,
                     kNameDistanceSum,
                     kNameMapqAF,
                     kNameRefMapqAF,
                     kNameSupport,
                     kNameWeightedDepth,
                     kNameMapqMin,
                     kNameMapqMean,
                     kNameMapqLT60Ratio,
                     kNameMapqLT30Ratio,
                     kNameRefNonHpMapqLT60Ratio,
                     kNameAltLen,
                     kNameRefDistanceMean,
                     kNameRefMapqMean,
                     kNameMapqLT40Ratio,
                     kNameRefMapqLT30Ratio,
                     kNameRefMapqLT40Ratio,
                     kNameNumAlt,
                     kNameADTL);
}

vec<std::string> GermlineWgsSimplexVcfFeatures() {
  return ToStringVec(kNameChrom,
                     kNamePos,
                     kNameRef,
                     kNameAlt,
                     kNamePopAF,
                     kNameVariantDensity,
                     kNameGenotype,
                     kNameQual,
                     kNameGq,
                     kNameNormalDP,
                     kNamePre2BpContext,
                     kNamePost2BpContext,
                     kNameRefAD,
                     kNameAltAD,
                     kNameAltAD2,
                     kNameRefAdAF,
                     kNameAltAdAF,
                     kNameAltAd2AF,
                     kNameUniq3mers,
                     kNameUniq4mers,
                     kNameUniq5mers,
                     kNameUniq6mers,
                     kNameHomopolymer,
                     kNameDirepeat,
                     kNameTrirepeat,
                     kNameQuadrepeat,
                     kNameHapcomp,
                     kNameHapdom,
                     kNameMqRankSum,
                     kNameSor,
                     kNameRu,
                     kNameRpaRef,
                     kNameRpaAlt,
                     kNameStr,
                     kNameAtInterest);
}

vec<std::string> GermlineWgsSimplexSnvScoringNames() {
  return ToStringVec(kNameRef,
                     kNameAlt,
                     kNameRefSupport,
                     kNameSupport,
                     kNameRefMapqMean,
                     kNameMapqMean,
                     kNameDistanceMean,
                     kNameRefDistanceMean,
                     kNamePopAF,
                     kNameVariantDensity,
                     kNameQual,
                     kNameGq,
                     kNameNormalDP,
                     kNameRefAD,
                     kNameAltAD,
                     kNameAltAD2,
                     kNameRefAdAF,
                     kNameAltAdAF,
                     kNameAltAd2AF,
                     kNameHapcomp,
                     kNameHapdom,
                     kNameMqRankSum,
                     kNameSor,
                     kNameGenotype,
                     kNameMapqLT30Ratio,
                     kNameRefMapqLT30Ratio,
                     kNameMapqLT40Ratio,
                     kNameRefMapqLT40Ratio,
                     kNameMapqAF);
}

vec<std::string> GermlineWgsSimplexIndelScoringNames() {
  return ToStringVec(kNameRef,
                     kNameAlt,
                     kNameRefNonHpSupport,
                     kNameAltLen,
                     kNameRefSupport,
                     kNameDistanceMin,
                     kNameDistanceMax,
                     kNameDistanceMean,
                     kNameDistanceSum,
                     kNameMapqAF,
                     kNameRefMapqAF,
                     kNameWeightedDepth,
                     kNamePopAF,
                     kNamePre2BpContext,
                     kNamePost2BpContext,
                     kNameHomopolymer,
                     kNameDirepeat,
                     kNameTrirepeat,
                     kNameQuadrepeat,
                     kNameHapcomp,
                     kNameHapdom,
                     kNameMqRankSum,
                     kNameSor,
                     kNameRu,
                     kNameRpaRef,
                     kNameRpaAlt,
                     kNameStr,
                     kNameUniq3mers,
                     kNameUniq4mers,
                     kNameUniq5mers,
                     kNameUniq6mers,
                     kNameMapqMin,
                     kNameMapqMean,
                     kNameMapqLT60Ratio,
                     kNameMapqLT30Ratio,
                     kNameRefNonHpMapqLT60Ratio,
                     kNameADTL,
                     kNameNumAlt,
                     kNameNormalDP,
                     kNameRefAD,
                     kNameAltAD,
                     kNameAltAD2,
                     kNameRefAdAF,
                     kNameAltAdAF,
                     kNameAltAd2AF,
                     kNameQual,
                     kNameGq,
                     kNameGenotype);
}

vec<std::string> GermlineWgsSimplexSnvCategoricalScoringNames() {
  return ToStringVec(kNameRef, kNameAlt, kNameGenotype);
}

vec<std::string> GermlineWgsSimplexIndelCategoricalScoringNames() {
  return ToStringVec(kNameRef, kNameAlt, kNamePre2BpContext, kNamePost2BpContext, kNameGenotype, kNameRu, kNameStr);
}

vec<std::string> TumorNormalWgsBamFeatures() {
  return ToStringVec(kNameChrom,
                     kNamePos,
                     kNameRef,
                     kNameAlt,
                     AddTumorPrefix(kNameSupport),
                     AddTumorPrefix(kNameRefSupport),
                     AddTumorPrefix(kNameDistanceMean),
                     AddTumorPrefix(kNameMapqMean),
                     AddTumorPrefix(kNameBaseqMean),
                     AddNormalPrefix(kNameSupport),
                     AddNormalPrefix(kNameRefSupport),
                     AddNormalPrefix(kNameDistanceMean),
                     AddNormalPrefix(kNameMapqMean),
                     AddNormalPrefix(kNameBaseqMean),
                     kNameMapqMean,
                     kNameMapqMin,
                     kNameMapqLT60Ratio,
                     kNameMapqLT40Ratio,
                     kNameMapqLT30Ratio,
                     kNameMapqLT20Ratio,
                     kNameRefMapqMean,
                     kNameRefMapqMin,
                     kNameRefMapqLT60Ratio,
                     kNameRefMapqLT40Ratio,
                     kNameRefMapqLT30Ratio,
                     kNameRefMapqLT20Ratio,
                     kNameBaseqMean,
                     kNameBaseqMin,
                     kNameRefBaseqMean,
                     AddTumorPrefix(kNameDel1bpAF),
                     AddTumorPrefix(kNameDel2bpAF),
                     AddTumorPrefix(kNameDel3bpAF),
                     AddTumorPrefix(kNameDelGt3bpAF),
                     AddTumorPrefix(kNameIns1bpAF),
                     AddTumorPrefix(kNameIns2bpAF),
                     AddTumorPrefix(kNameIns3bpAF),
                     AddTumorPrefix(kNameInsGt3bpAF),
                     AddNormalPrefix(kNameDel1bpAF),
                     AddNormalPrefix(kNameDel2bpAF),
                     AddNormalPrefix(kNameDel3bpAF),
                     AddNormalPrefix(kNameDelGt3bpAF),
                     AddNormalPrefix(kNameIns1bpAF),
                     AddNormalPrefix(kNameIns2bpAF),
                     AddNormalPrefix(kNameIns3bpAF),
                     AddNormalPrefix(kNameInsGt3bpAF),
                     AddTumorPrefix(kNameDel1bpIns1bpAFLogRatio),
                     AddTumorPrefix(kNameDel1bp2bpAFLogRatio),
                     AddTumorPrefix(kNameDel2bp3bpAFLogRatio),
                     AddTumorPrefix(kNameDel3bpGt3bpAFLogRatio),
                     AddTumorPrefix(kNameIns1bp2bpAFLogRatio),
                     AddTumorPrefix(kNameIns2bp3bpAFLogRatio),
                     AddTumorPrefix(kNameIns3bpGt3bpAFLogRatio),
                     AddTumorPrefix(kNameSnvAF),
                     kNameDistanceMean,
                     kNameRefDistanceMean,
                     AddTumorPrefix(kNameDuplexAF),
                     AddNormalPrefix(kNameDuplexAF),
                     kNameBamTnAfRatio,
                     AddTumorPrefix(kNameSupportReverse),
                     AddTumorPrefix(kNameAlignmentBias),
                     kNameContext);
}

vec<std::string> TumorNormalWgsVcfFeatures() {
  return ToStringVec(kNameChrom,
                     kNamePos,
                     kNameRef,
                     kNameAlt,
                     kNameTumorDP,
                     kNameNormalDP,
                     kNameNalod,
                     kNameNlod,
                     kNameTlod,
                     kNameMpos,
                     kNamePopAF,
                     kNameHapcomp,
                     kNameHapdom,
                     kNameStr,
                     kNameRu,
                     kNameRpaRef,
                     kNameRpaAlt,
                     kNameSubtypeIndex,
                     kNameVariantType,
                     kNameUniq3mers,
                     kNameUniq4mers,
                     kNameUniq5mers,
                     kNameUniq6mers,
                     kNamePre2BpContext,
                     kNamePost2BpContext,
                     kNameHomopolymer,
                     kNameDirepeat,
                     kNameTrirepeat,
                     kNameQuadrepeat,
                     kNameIsInRepeat,
                     kNameRepeatLen,
                     kNameRepeatUnitLen,
                     kNameIsUnitSizeMatch,
                     kNamePonOcc,
                     kNamePonOccRatio,
                     kNamePonMeanDuplexAf,
                     kNamePonMinDuplexAf,
                     kNamePonMaxDuplexAf);
}

vec<std::string> TumorNormalWgsScoringNames() {
  return ToStringVec(kNameNalod,
                     kNameNlod,
                     kNameTlod,
                     kNameMpos,
                     kNamePopAF,
                     kNameHapcomp,
                     kNameHapdom,
                     kNameStr,
                     kNameRu,
                     kNameRpaRef,
                     kNameRpaAlt,
                     kNameSubtypeIndex,
                     kNameVariantType,
                     kNameUniq3mers,
                     kNameUniq4mers,
                     kNameUniq5mers,
                     kNameUniq6mers,
                     kNamePre2BpContext,
                     kNamePost2BpContext,
                     kNameHomopolymer,
                     kNameDirepeat,
                     kNameTrirepeat,
                     kNameQuadrepeat,
                     kNameIsInRepeat,
                     kNameRepeatLen,
                     kNameRepeatUnitLen,
                     kNameIsUnitSizeMatch,
                     AddTumorPrefix(kNameSupport),
                     AddTumorPrefix(kNameRefSupport),
                     AddTumorPrefix(kNameDistanceMean),
                     AddTumorPrefix(kNameMapqMean),
                     AddTumorPrefix(kNameBaseqMean),
                     AddNormalPrefix(kNameSupport),
                     AddNormalPrefix(kNameRefSupport),
                     AddNormalPrefix(kNameDistanceMean),
                     AddNormalPrefix(kNameMapqMean),
                     AddNormalPrefix(kNameBaseqMean),
                     kNameMapqMean,
                     kNameMapqMin,
                     kNameMapqLT60Ratio,
                     kNameMapqLT40Ratio,
                     kNameMapqLT30Ratio,
                     kNameMapqLT20Ratio,
                     kNameRefMapqMean,
                     kNameRefMapqMin,
                     kNameRefMapqLT60Ratio,
                     kNameRefMapqLT40Ratio,
                     kNameRefMapqLT30Ratio,
                     kNameRefMapqLT20Ratio,
                     kNameBaseqMean,
                     kNameBaseqMin,
                     kNameRefBaseqMean,
                     AddTumorPrefix(kNameDel1bpAF),
                     AddTumorPrefix(kNameDel2bpAF),
                     AddTumorPrefix(kNameDel3bpAF),
                     AddTumorPrefix(kNameDelGt3bpAF),
                     AddTumorPrefix(kNameIns1bpAF),
                     AddTumorPrefix(kNameIns2bpAF),
                     AddTumorPrefix(kNameIns3bpAF),
                     AddTumorPrefix(kNameInsGt3bpAF),
                     AddNormalPrefix(kNameDel1bpAF),
                     AddNormalPrefix(kNameDel2bpAF),
                     AddNormalPrefix(kNameDel3bpAF),
                     AddNormalPrefix(kNameDelGt3bpAF),
                     AddNormalPrefix(kNameIns1bpAF),
                     AddNormalPrefix(kNameIns2bpAF),
                     AddNormalPrefix(kNameIns3bpAF),
                     AddNormalPrefix(kNameInsGt3bpAF),
                     AddTumorPrefix(kNameDel1bpIns1bpAFLogRatio),
                     AddTumorPrefix(kNameDel1bp2bpAFLogRatio),
                     AddTumorPrefix(kNameDel2bp3bpAFLogRatio),
                     AddTumorPrefix(kNameDel3bpGt3bpAFLogRatio),
                     AddTumorPrefix(kNameIns1bp2bpAFLogRatio),
                     AddTumorPrefix(kNameIns2bp3bpAFLogRatio),
                     AddTumorPrefix(kNameIns3bpGt3bpAFLogRatio),
                     AddTumorPrefix(kNameSnvAF),
                     kNameDistanceMean,
                     kNameRefDistanceMean,
                     AddTumorPrefix(kNameDuplexAF),
                     AddNormalPrefix(kNameDuplexAF),
                     kNameBamTnAfRatio,
                     AddTumorPrefix(kNameAlignmentBias),
                     kNameContext,
                     kNamePonOcc,
                     kNamePonOccRatio,
                     kNamePonMeanDuplexAf,
                     kNamePonMinDuplexAf,
                     kNamePonMaxDuplexAf);
}

vec<std::string> TumorNormalWgsCategoricalNames() {
  return ToStringVec(
      kNameRu, kNameStr, kNameSubtypeIndex, kNameVariantType, kNamePre2BpContext, kNamePost2BpContext, kNameContext);
}

vec<std::string> TumorOnlyWgsBamFeatures() {
  return ToStringVec(kNameChrom,
                     kNamePos,
                     kNameRef,
                     kNameAlt,
                     kNameSupport,
                     kNameSupportReverse,
                     kNameAlignmentBias,
                     kNameRefSupport,
                     kNameMapqMean,
                     kNameMapqMin,
                     kNameMapqLT60Ratio,
                     kNameMapqLT40Ratio,
                     kNameMapqLT30Ratio,
                     kNameMapqLT20Ratio,
                     kNameRefMapqMean,
                     kNameRefMapqMin,
                     kNameRefMapqLT60Ratio,
                     kNameRefMapqLT40Ratio,
                     kNameRefMapqLT30Ratio,
                     kNameRefMapqLT20Ratio,
                     kNameBaseqMean,
                     kNameBaseqMin,
                     kNameRefBaseqMean,
                     kNameDistanceMean,
                     kNameRefDistanceMean,
                     kNameDuplexAF,
                     kNameDuplex,
                     kNameDuplexLowbq,
                     kNameSimplex,
                     kNameNonDuplex,
                     kNameStrandBias,
                     kNamePlusOnly,
                     kNameMinusOnly,
                     kNameWeightedDepth,
                     kNameWeightedScore,
                     kNameFamilysizeSum,
                     kNameFamilysizeMean,
                     kNameContext,
                     kNameSubtypeIndex);
}

vec<std::string> TumorOnlyWgsVcfFeatures() {
  return ToStringVec(kNameChrom,
                     kNamePos,
                     kNameRef,
                     kNameAlt,
                     kNameVariantType,
                     kNameTlod,
                     kNameMpos,
                     kNamePre2BpContext,
                     kNamePost2BpContext,
                     kNameHomopolymer,
                     kNameDirepeat,
                     kNameTrirepeat,
                     kNameQuadrepeat,
                     kNameUniq3mers,
                     kNameUniq4mers,
                     kNameUniq5mers,
                     kNameUniq6mers,
                     kNamePopAF,
                     kNameRefAD,
                     kNameAltAD,
                     kNameAltAD2,
                     kNameRefAdAF,
                     kNameAltAdAF,
                     kNameAltAd2AF,
                     kNameVariantDensity,
                     kNameHapcomp,
                     kNameHapdom,
                     kNameRu,
                     kNameRpaRef,
                     kNameRpaAlt,
                     kNameStr,
                     kNameGcContent,
                     kNamePonOcc,
                     kNamePonOccRatio,
                     kNamePonMeanDuplexAf,
                     kNamePonMinDuplexAf,
                     kNamePonMaxDuplexAf);
}

vec<std::string> TumorOnlyWgsScoringNames() {
  // Single-sample WGS scoring schema. Mirrors `tumor-normal-wgs` with all normal-sample features dropped
  // (tumor/normal prefixes removed, normal LOD and tumor/normal AF-ratio features removed) and the Panel
  // of Normals features appended. PON features default to 0 when no PON VCF is supplied (popAF-style).
  return ToStringVec(kNameTlod,
                     kNameMpos,
                     kNamePopAF,
                     kNameHapcomp,
                     kNameHapdom,
                     kNameStr,
                     kNameRu,
                     kNameRpaRef,
                     kNameRpaAlt,
                     kNameSubtypeIndex,
                     kNameVariantType,
                     kNameUniq3mers,
                     kNameUniq4mers,
                     kNameUniq5mers,
                     kNameUniq6mers,
                     kNamePre2BpContext,
                     kNamePost2BpContext,
                     kNameHomopolymer,
                     kNameDirepeat,
                     kNameTrirepeat,
                     kNameQuadrepeat,
                     kNameSupport,
                     kNameRefSupport,
                     kNameMapqMean,
                     kNameMapqMin,
                     kNameMapqLT60Ratio,
                     kNameMapqLT40Ratio,
                     kNameMapqLT30Ratio,
                     kNameMapqLT20Ratio,
                     kNameRefMapqMean,
                     kNameRefMapqMin,
                     kNameRefMapqLT60Ratio,
                     kNameRefMapqLT40Ratio,
                     kNameRefMapqLT30Ratio,
                     kNameRefMapqLT20Ratio,
                     kNameBaseqMean,
                     kNameBaseqMin,
                     kNameRefBaseqMean,
                     kNameDistanceMean,
                     kNameRefDistanceMean,
                     kNameDuplexAF,
                     kNameAlignmentBias,
                     kNameContext,
                     kNamePonOcc,
                     kNamePonOccRatio,
                     kNamePonMeanDuplexAf,
                     kNamePonMinDuplexAf,
                     kNamePonMaxDuplexAf);
}

vec<std::string> TumorOnlyWgsCategoricalNames() {
  return ToStringVec(
      kNameRu, kNameStr, kNameSubtypeIndex, kNameVariantType, kNamePre2BpContext, kNamePost2BpContext, kNameContext);
}

vec<std::string> PanelOfNormalsBamFeatures() {
  return ToStringVec(kNameChrom, kNamePos, kNameRef, kNameAlt, kNameDuplexAF, kNameMapqMean);
}

vec<std::string> GermlineTaggingBamFeatures() {
  return ToStringVec(kNameChrom,
                     kNamePos,
                     kNameRef,
                     kNameAlt,
                     AddNormalPrefix(kNameMapqMean),
                     AddNormalPrefix(kNameBaseqMean),
                     AddNormalPrefix(kNameDistanceMean),
                     AddNormalPrefix(kNameSupport),
                     AddNormalPrefix(kNameRefSupport),
                     AddNormalPrefix(kNameDuplexAF),
                     kNameContext);
}

vec<std::string> GermlineTaggingVcfFeatures() {
  return ToStringVec(kNameChrom,
                     kNamePos,
                     kNameRef,
                     kNameAlt,
                     kNameTumorDP,
                     kNameNormalDP,
                     kNamePopAF,
                     kNameNormalAltAD,
                     kNameNormalAF,
                     kNameSubtypeIndex,
                     kNameVariantType,
                     kNameNalod,
                     kNameNlod,
                     kNamePre2BpContext,
                     kNamePost2BpContext,
                     kNameUniq3mers,
                     kNameUniq4mers,
                     kNameUniq5mers,
                     kNameUniq6mers,
                     kNameHomopolymer,
                     kNameDirepeat,
                     kNameTrirepeat,
                     kNameQuadrepeat,
                     kNameHapcomp,
                     kNameHapdom,
                     kNameRu,
                     kNameRpaRef,
                     kNameRpaAlt,
                     kNameStr);
}

vec<std::string> GermlineTaggingScoringNames() {
  return ToStringVec(kNameRef,
                     kNameAlt,
                     AddNormalPrefix(kNameMapqMean),
                     AddNormalPrefix(kNameBaseqMean),
                     AddNormalPrefix(kNameDistanceMean),
                     AddNormalPrefix(kNameSupport),
                     AddNormalPrefix(kNameRefSupport),
                     AddNormalPrefix(kNameDuplexAF),
                     kNameContext,
                     kNamePopAF,
                     kNameNormalAltAD,
                     kNameNormalAF,
                     kNameNormalDP,
                     kNameSubtypeIndex,
                     kNameVariantType,
                     kNameNalod,
                     kNameNlod,
                     kNamePre2BpContext,
                     kNamePost2BpContext,
                     kNameUniq3mers,
                     kNameUniq4mers,
                     kNameUniq5mers,
                     kNameUniq6mers,
                     kNameHomopolymer,
                     kNameDirepeat,
                     kNameTrirepeat,
                     kNameQuadrepeat,
                     kNameHapcomp,
                     kNameHapdom,
                     kNameRu,
                     kNameRpaRef,
                     kNameRpaAlt,
                     kNameStr);
}

vec<std::string> GermlineTaggingCategoricalNames() {
  return ToStringVec(kNameRef,
                     kNameAlt,
                     kNameVariantType,
                     kNameSubtypeIndex,
                     kNamePre2BpContext,
                     kNamePost2BpContext,
                     kNameContext,
                     kNameRu,
                     kNameStr);
}

using Json = nlohmann::json;

SnvIndelPaths GetSnvIndelPaths(const std::vector<fs::path>& paths) {
  SnvIndelPaths out_paths{};
  // distinguish SNV and indel paths based on file names
  for (const auto& p : paths) {
    const auto file_name = p.filename().string();
    const bool snv_found = file_name.find("snv") != std::string::npos;
    const bool indel_found = file_name.find("indel") != std::string::npos;
    if (snv_found && !indel_found) {
      out_paths.snv_path = p;
    } else if (indel_found && !snv_found) {
      out_paths.indel_path = p;
    } else {
      // file contains both or neither "snv"/"indel" in name; skip
    }
  }
  return out_paths;
}

SVCConfig::SVCConfig(const Workflow target_workflow) : SVCConfig(ConfigProfileForWorkflow(target_workflow)) {
}

SVCConfig::SVCConfig(const ConfigProfile profile) {
  using enum ConfigProfile;
  switch (profile) {
    case kCustom: {
      ConfigureCustomWorkflow();
      break;
    }
    case kGermlineWgsDuplex: {
      ConfigureGermlineWgsDuplexProfile();
      break;
    }
    case kGermlineWgsSimplex: {
      ConfigureGermlineWgsSimplexProfile();
      break;
    }
    case kTumorOnlyTe: {
      ConfigureTumorOnlyTeWorkflow();
      break;
    }
    case kTumorOnlyWgs: {
      ConfigureTumorOnlyWgsWorkflow();
      break;
    }
    case kTumorNormalWgs: {
      ConfigureTumorNormalWgsWorkflow();
      break;
    }
    case kGermlineTagging: {
      ConfigureGermlineTaggingWorkflow();
      break;
    }
    case kGermlineMitochondria: {
      ConfigureMitochondriaWorkflow();
      break;
    }
    case kPanelOfNormals: {
      ConfigurePanelOfNormalsWorkflow();
      break;
    }
    default: {
      throw error::Error("Unknown config profile value {}", enum_util::FormatEnumName(profile));
    }
  }
  ConfigureFeatureCols();
}

void SVCConfig::ConfigureCustomWorkflow() {
  bam_feature_names = DefaultBamFeatures();
  vcf_feature_names = DefaultVcfFeatures();
  workflow = Workflow::kCustom;
  scoring_names = kDefaultScoringNames;
  categorical_names = kDefaultCategoricalNames;
  n_classes = 0;
  min_mapq = 1;
  min_bq = kMinBaseQuality;
  min_dist = 2;
  min_family_size = 1;
  filter_homopolymer = HomopolymerFilter::kAlignmentEnd;
  min_homopolymer_length = 4;
  sequencing_protocol = SequencingProtocol::kDuplex;
  duplex_lowbq = DuplexLowbqMode::kExclude;
  use_vcf_features = true;
}

void SVCConfig::ConfigureGermlineWgsDuplexProfile() {
  bam_feature_names = GermlineWgsDuplexBamFeatures();
  vcf_feature_names = GermlineWgsDuplexVcfFeatures();
  workflow = Workflow::kGermlineWgs;
  snv_scoring_names = GermlineWgsDuplexSnvScoringNames();
  indel_scoring_names = GermlineWgsDuplexIndelScoringNames();
  indel_categorical_names = GermlineWgsDuplexIndelCategoricalScoringNames();
  snv_categorical_names = GermlineWgsDuplexSnvCategoricalScoringNames();
  n_classes = 4;
  min_mapq = 1;
  min_bq = kMinBaseQuality;
  min_dist = 2;
  min_family_size = kSimplexReadFamilySize;
  filter_homopolymer = HomopolymerFilter::kAlignmentEnd;
  min_homopolymer_length = 4;
  sequencing_protocol = SequencingProtocol::kDuplex;
  duplex_lowbq = DuplexLowbqMode::kInclude;
  use_vcf_features = true;
  snv_iterations = kGermlineWgsLargeSnvIterations;
  indel_iterations = kGermlineWgsLargeIndelIterations;
  normalize_features = FeatureNormalization::kMedianDp;
  decode_yc = YcDecodeMethod::kConsensus;
  min_base_type = BaseType::kSimplex;
  snv_model_lgbm_params = kGermlineWgsLargeSnvModelLgbmParams;
  indel_model_lgbm_params = kGermlineWgsLargeIndelModelLgbmParams;
  snv_model_lgbm_prediction_params =
      "num_threads=1 force_row_wise=true deterministic=true seed=1238845 pred_early_stop=true";
  indel_model_lgbm_prediction_params =
      "num_threads=1 force_row_wise=true deterministic=true seed=1238845 pred_early_stop=true";
  PopulateGermlineModelTrainingDefaults(*this);
}

void SVCConfig::ConfigureGermlineWgsSimplexProfile() {
  bam_feature_names = GermlineWgsSimplexBamFeatures();
  vcf_feature_names = GermlineWgsSimplexVcfFeatures();
  workflow = Workflow::kGermlineWgs;
  snv_scoring_names = GermlineWgsSimplexSnvScoringNames();
  indel_scoring_names = GermlineWgsSimplexIndelScoringNames();
  indel_categorical_names = GermlineWgsSimplexIndelCategoricalScoringNames();
  snv_categorical_names = GermlineWgsSimplexSnvCategoricalScoringNames();
  n_classes = 4;
  min_mapq = 1;
  min_bq = kMinBaseQuality;
  min_dist = 2;
  min_family_size = kSimplexReadFamilySize;
  filter_homopolymer = HomopolymerFilter::kAlignmentEnd;
  min_homopolymer_length = 4;
  sequencing_protocol = SequencingProtocol::kSimplex;
  duplex_lowbq = DuplexLowbqMode::kExclude;
  use_vcf_features = true;
  snv_iterations = kGermlineWgsLargeSnvIterations;
  indel_iterations = kGermlineWgsLargeIndelIterations;
  normalize_features = FeatureNormalization::kMedianDp;
  decode_yc = YcDecodeMethod::kNone;
  min_base_type = BaseType::kSimplex;
  snv_model_lgbm_params = kGermlineWgsLargeSnvModelLgbmParams;
  indel_model_lgbm_params = kGermlineWgsLargeIndelModelLgbmParams;
  snv_model_lgbm_prediction_params =
      "num_threads=1 force_row_wise=true deterministic=true seed=1238845 pred_early_stop=true";
  indel_model_lgbm_prediction_params =
      "num_threads=1 force_row_wise=true deterministic=true seed=1238845 pred_early_stop=true";
  PopulateGermlineModelTrainingDefaults(*this);
}

void SVCConfig::ConfigureTumorOnlyTeWorkflow() {
  bam_feature_names = TumorOnlyTeBamFeatures();
  vcf_feature_names = kTumorOnlyTeVcfFeatures;
  workflow = Workflow::kTumorOnlyTe;
  scoring_names = TumorOnlyTeScoringNames();
  categorical_names = TumorOnlyTeCategoricalNames();
  n_classes = 1;
  min_mapq = kTumorOnlyTeMinMapq;
  min_bq = kSimplexMinBaseQuality;
  min_dist = 0;
  max_variants_per_read = kTumorOnlyTeMaxVariantsPerRead;
  min_family_size = 3;
  filter_homopolymer = HomopolymerFilter::kNone;
  min_homopolymer_length = kNonGermlineMinHomopolymerLength;
  sequencing_protocol = SequencingProtocol::kUmi;
  duplex_lowbq = DuplexLowbqMode::kExclude;
  min_alt_counts = 3;
  // min AF for FFPE: 0.01, for ctDNA: 0.0
  min_af = kTumorOnlyTeMinAf;
  min_phased_af = kSomaticMinPhasedAf;
  max_phased_af = kSomaticMaxPhasedAf;
  min_weighted_counts = 4;
  hotspot_min_weighted_counts = 2;
  min_ml_score = kSomaticDefaultMlScore;
  hotspot_min_ml_score = kSomaticHotspotMlScore;
  phased = false;
  use_vcf_features = false;
  iterations = kTumorOnlyTeIterations;
  snv_iterations = 0;
  indel_iterations = 0;
  snv_min_ml_score = kSomaticDefaultMlScore;
  indel_min_ml_score = kSomaticDefaultMlScore;
  min_tumor_support = 1;
  decode_yc = YcDecodeMethod::kNone;
  model_lgbm_params = "objective=binary boosting=gbdt learning_rate=0.01 seed=1238845 num_leaves=3 num_classes=1";
  model_lgbm_prediction_params =
      "num_threads=1 force_row_wise=true deterministic=true seed=1238845 pred_early_stop=true";
}

void SVCConfig::ConfigureTumorOnlyWgsWorkflow() {
  bam_feature_names = TumorOnlyWgsBamFeatures();
  vcf_feature_names = TumorOnlyWgsVcfFeatures();
  workflow = Workflow::kTumorOnlyWgs;
  scoring_names = TumorOnlyWgsScoringNames();
  categorical_names = TumorOnlyWgsCategoricalNames();
  n_classes = 1;
  min_mapq = 1;
  min_bq = kConcordantMinBaseQuality;
  min_dist = 2;
  min_family_size = kDuplexReadFamilySize;
  filter_homopolymer = HomopolymerFilter::kNone;
  min_homopolymer_length = kNonGermlineMinHomopolymerLength;
  sequencing_protocol = SequencingProtocol::kDuplex;
  duplex_lowbq = DuplexLowbqMode::kInclude;
  use_vcf_features = true;
  iterations = kTumorOnlyWgsIterations;
  normalize_features = FeatureNormalization::kNone;
  // Single-sample alt-support gate for training and filtering. Tumor-only-wgs uses the unprefixed `support`
  // feature (no tumor/normal split).
  min_support = 1;
  snv_min_ml_score = kTumorOnlyWgsSnvMinMlScore;
  indel_min_ml_score = kTumorOnlyWgsIndelMinMlScore;
  min_af = kTumorOnlyWgsMinAf;
  max_indel_size = kTumorOnlyWgsMaxIndelSize;
  decode_yc = YcDecodeMethod::kConsensus;
  min_base_type = BaseType::kConcordant;
  model_lgbm_params = "objective=binary boosting=gbdt learning_rate=0.005 seed=1238845 num_classes=1";
  model_lgbm_prediction_params = "num_threads=1 force_row_wise=true deterministic=true seed=1238845";
  // Keyed by {sample_type}-{aligner} (see
  // resources/model-tumor-only-wgs-<sample_type>-<aligner>.txt.gz).
  model_thresholds = {
      {ModelThresholdsKey(SampleType::kCellLine, AlignerType::kBwa),
       {.snv_min_ml_score = kTumorOnlyWgsCellLineSnvMlScore,
        .indel_min_ml_score = kTumorOnlyWgsCellLineIndelMlScore,
        .min_af = kTumorOnlyWgsCellLineMinAf}},
      {ModelThresholdsKey(SampleType::kFfpe, AlignerType::kBwa),
       {.snv_min_ml_score = kTumorOnlyWgsFfpeSnvMlScore,
        .indel_min_ml_score = kTumorOnlyWgsFfpeIndelMlScore,
        .min_af = kTumorOnlyWgsFfpeMinAf}},
  };
}

void SVCConfig::ConfigureTumorNormalWgsWorkflow() {
  bam_feature_names = TumorNormalWgsBamFeatures();
  vcf_feature_names = TumorNormalWgsVcfFeatures();
  workflow = Workflow::kTumorNormalWgs;
  scoring_names = TumorNormalWgsScoringNames();
  categorical_names = TumorNormalWgsCategoricalNames();
  n_classes = 1;
  min_mapq = 1;
  min_bq = kConcordantMinBaseQuality;
  min_dist = 2;
  min_family_size = kDuplexReadFamilySize;
  filter_homopolymer = HomopolymerFilter::kNone;
  min_homopolymer_length = kNonGermlineMinHomopolymerLength;
  sequencing_protocol = SequencingProtocol::kDuplex;
  duplex_lowbq = DuplexLowbqMode::kInclude;
  min_alt_counts = 3;
  min_phased_af = kSomaticMinPhasedAf;
  max_phased_af = kSomaticMaxPhasedAf;
  min_weighted_counts = 4;
  hotspot_min_weighted_counts = 2;
  min_ml_score = kSomaticDefaultMlScore;
  hotspot_min_ml_score = kSomaticHotspotMlScore;
  phased = false;
  use_vcf_features = true;
  iterations = kTumorNormalWgsIterations;
  normalize_features = FeatureNormalization::kNone;
  snv_min_ml_score = kTumorNormalWgsSnvMinMlScore;
  indel_min_ml_score = kTumorNormalWgsIndelMinMlScore;
  min_tumor_support = 1;
  max_normal_support = 2;
  min_tumor_af = kTumorNormalWgsFfpeMinTumorAf;
  max_indel_size = kTumorNormalWgsMaxIndelSize;
  min_dp_ratio = kTumorNormalWgsMinDpRatio;
  decode_yc = YcDecodeMethod::kConsensus;
  min_base_type = BaseType::kConcordant;
  model_lgbm_params = "objective=binary boosting=gbdt learning_rate=0.005 seed=1238845 num_classes=1";
  model_lgbm_prediction_params = "num_threads=1 force_row_wise=true deterministic=true seed=1238845";
  model_thresholds = {
      {ModelThresholdsKey(SampleType::kFfpe, AlignerType::kBwa),
       {.snv_min_ml_score = kTumorNormalWgsFfpeSnvMlScore,
        .indel_min_ml_score = kTumorNormalWgsFfpeIndelMlScore,
        .min_tumor_af = kTumorNormalWgsFfpeMinTumorAf}},
      {ModelThresholdsKey(SampleType::kCellLine, AlignerType::kBwa),
       {.snv_min_ml_score = kTumorNormalWgsCellLineSnvMlScore,
        .indel_min_ml_score = kTumorNormalWgsCellLineIndelMlScore,
        .min_tumor_af = kTumorNormalWgsCellLineMinTumorAf}},
  };
}

void SVCConfig::ConfigureGermlineTaggingWorkflow() {
  bam_feature_names = GermlineTaggingBamFeatures();
  vcf_feature_names = GermlineTaggingVcfFeatures();
  workflow = Workflow::kGermlineTagging;
  scoring_names = GermlineTaggingScoringNames();
  categorical_names = GermlineTaggingCategoricalNames();
  n_classes = 3;
  min_mapq = 1;
  min_bq = kConcordantMinBaseQuality;
  min_dist = 2;
  min_family_size = 2;
  filter_homopolymer = HomopolymerFilter::kNone;
  min_homopolymer_length = kNonGermlineMinHomopolymerLength;
  sequencing_protocol = SequencingProtocol::kDuplex;
  duplex_lowbq = DuplexLowbqMode::kExclude;
  use_vcf_features = true;
  iterations = kGermlineTaggingModelIterations;
  normalize_features = FeatureNormalization::kNone;
  snv_min_ml_score = 0.0F;
  indel_min_ml_score = 0.0F;
  min_tumor_support = 1;
  decode_yc = YcDecodeMethod::kConsensus;
  min_base_type = BaseType::kConcordant;
  model_lgbm_params =
      "objective=multiclass boosting=gbdt metric=multi_logloss seed=1238845 learning_rate=0.01 bagging_fraction=0.9 "
      "bagging_freq=1 feature_fraction=0.5 num_leaves=64 is_unbalance=true num_classes=3 ";
  model_lgbm_prediction_params =
      "num_threads=1 force_row_wise=true deterministic=true seed=1238845 pred_early_stop=true";
}

void SVCConfig::ConfigureMitochondriaWorkflow() {
  bam_feature_names = GermlineMitochondriaBamFeatures();
  vcf_feature_names = GermlineMitochondriaVcfFeatures();
  workflow = Workflow::kGermlineMitochondria;
  scoring_names = {};
  categorical_names = {};
  n_classes = 0;
  min_mapq = 1;
  min_bq = kMinBaseQuality;
  min_dist = kGermlineMitochondriaMinDist;
  min_family_size = 1;
  filter_homopolymer = HomopolymerFilter::kNone;
  min_homopolymer_length = 0;
  sequencing_protocol = SequencingProtocol::kDuplex;
  use_vcf_features = false;
  normalize_features = FeatureNormalization::kNone;
  decode_yc = YcDecodeMethod::kConsensus;
  min_base_type = BaseType::kConcordant;
  min_af = kGermlineMitochondriaMinAf;
  min_dp = kGermlineMitochondriaMinDp;
  min_hom_alt_af = kGermlineMitochondriaMinHomAltAf;
  duplex_lowbq = DuplexLowbqMode::kExclude;
}

void SVCConfig::ConfigurePanelOfNormalsWorkflow() {
  bam_feature_names = PanelOfNormalsBamFeatures();
  workflow = Workflow::kPanelOfNormals;
  n_classes = 0;
  min_mapq = 1;
  min_bq = kConcordantMinBaseQuality;
  min_dist = 2;
  min_family_size = kDuplexReadFamilySize;
  filter_homopolymer = HomopolymerFilter::kNone;
  min_homopolymer_length = kNonGermlineMinHomopolymerLength;
  sequencing_protocol = SequencingProtocol::kDuplex;
  duplex_lowbq = DuplexLowbqMode::kInclude;
  use_vcf_features = false;
  decode_yc = YcDecodeMethod::kConsensus;
  min_base_type = BaseType::kConcordant;
}

void SVCConfig::ApplyGermlineModelTrainingDefaults(const GermlineModelType model_type) {
  if (workflow != Workflow::kGermlineWgs) {
    throw error::Error("Germline model type '{}' can only be applied to the germline-wgs workflow",
                       GermlineModelTypeConfigKey(model_type));
  }

  const auto key = GermlineModelTypeConfigKey(model_type);
  const auto defaults_it = germline_model_type_training_defaults.find(key);
  if (defaults_it == germline_model_type_training_defaults.end()) {
    throw error::Error("Config profile does not define germline training defaults for model type '{}'", key);
  }

  const auto& defaults = defaults_it->second;
  snv_iterations = defaults.snv_iterations;
  indel_iterations = defaults.indel_iterations;
  snv_model_lgbm_params = defaults.snv_model_lgbm_params;
  indel_model_lgbm_params = defaults.indel_model_lgbm_params;
}

bool SVCConfig::HasVcfFeatureScoringCols() const {
  // determine whether any scoring columns are (derivatives of) VCF features
  return std::ranges::any_of(scoring_cols, IsVcfFeatureColumn) ||
         std::ranges::any_of(snv_scoring_cols, IsVcfFeatureColumn) ||
         std::ranges::any_of(indel_scoring_cols, IsVcfFeatureColumn);
}

void SVCConfig::ConfigureBamFeatureCols() {
  const bool allow_tn_prefix = FeatureNamesHaveSampleContext(workflow);
  // BAM feature names can be additionally attached with a sample context prefix, e.g. `tumor_` or `normal_`
  // So, we need to check for unsupported feature names with and without the prefixes.
  const auto& unsupported_fields = FindUnsupportedBamFeatureNames(bam_feature_names, allow_tn_prefix);
  if (!unsupported_fields.empty()) {
    throw error::Error("Unsupported BAM feature names: {}", string::Join(unsupported_fields, ", "));
  }
  feature_cols.reserve(bam_feature_names.size());
  for (const auto& col_name : bam_feature_names) {
    feature_cols.push_back(GetFeatureColumn(col_name, allow_tn_prefix));
  }
}

void SVCConfig::ConfigureVcfFeatureCols() {
  const bool allow_tn_prefix = FeatureNamesHaveSampleContext(workflow);
  // VCF feature names already have hard-coded sample context prefixes.
  // So, we check for unsupported feature names without inferring sample context, and throw an error if any unsupported
  // feature names are found.
  const auto unsupported_fields = FindUnsupportedVcfFeatureNames(vcf_feature_names);
  if (!unsupported_fields.empty()) {
    throw error::Error("Unsupported VCF feature names: {}", string::Join(unsupported_fields, ", "));
  }
  vcf_feature_cols.reserve(vcf_feature_names.size());
  for (const auto& col_name : vcf_feature_names) {
    // we still want to infer sample context if the workflow allows it
    vcf_feature_cols.push_back(GetFeatureColumn(col_name, allow_tn_prefix));
  }
}

void SVCConfig::ConfigureScoringCols() {
  // Extract the scoring column enums from the names provided in the config.
  // The enums are more compact and more efficient to work with than the strings.
  const bool allow_tn_prefix = FeatureNamesHaveSampleContext(workflow);
  if (!scoring_names.empty()) {
    const auto unsupported_fields = FindUnsupportedScoringFeatureNames(scoring_names, allow_tn_prefix);
    if (!unsupported_fields.empty()) {
      throw error::Error("Unsupported scoring feature names: {}", string::Join(unsupported_fields, ", "));
    }
    scoring_cols.reserve(scoring_names.size());
    for (const auto& col_name : scoring_names) {
      scoring_cols.push_back(GetFeatureColumn(col_name, allow_tn_prefix));
    }
  }
  // germline workflow's SNV model scoring columns
  if (!snv_scoring_names.empty()) {
    const auto unsupported_fields = FindUnsupportedScoringFeatureNames(snv_scoring_names, allow_tn_prefix);
    if (!unsupported_fields.empty()) {
      throw error::Error("Unsupported SNV scoring feature names: {}", string::Join(unsupported_fields, ", "));
    }
    snv_scoring_cols.reserve(snv_scoring_names.size());
    for (const auto& col_name : snv_scoring_names) {
      snv_scoring_cols.push_back(GetFeatureColumn(col_name, allow_tn_prefix));
    }
  }
  // germline workflow's indel model scoring columns
  if (!indel_scoring_names.empty()) {
    const auto unsupported_fields = FindUnsupportedScoringFeatureNames(indel_scoring_names, allow_tn_prefix);
    if (!unsupported_fields.empty()) {
      throw error::Error("Unsupported indel scoring feature names: {}", string::Join(unsupported_fields, ", "));
    }
    indel_scoring_cols.reserve(indel_scoring_names.size());
    for (const auto& col_name : indel_scoring_names) {
      indel_scoring_cols.push_back(GetFeatureColumn(col_name, allow_tn_prefix));
    }
  }
}

void SVCConfig::ConfigureFeatureCols() {
  ConfigureBamFeatureCols();
  ConfigureVcfFeatureCols();
  ConfigureScoringCols();
}

SVCConfigCollection JsonToConfigCollection(const fs::path& config_json) {
  if (!exists(config_json)) {
    throw error::Error("JSON file not found: {}", config_json);
  }
  SVCConfigCollection config_collection;
  std::ifstream fh(config_json);
  from_json(Json::parse(fh), config_collection);
  for (auto& [_, config] : config_collection.config_profiles) {
    NormalizeGermlineModelTrainingDefaults(config);
  }
  return config_collection;
}

namespace {
/**
 * @brief Extract and initialize a config profile from a JSON file.
 * @param config_json Path of JSON file. Path must exist.
 * @param profile Concrete config profile to load.
 * @return Config for the requested profile.
 * @throws xoos::error::Error if the file, profile, or declared workflow is invalid.
 */
SVCConfig JsonToConfigHelper(const fs::path& config_json, const ConfigProfile profile) {
  if (!exists(config_json)) {
    // Must check here whether file exists instead of defining the `--config` option with `check(CLI::ExistingFile)`.
    // By design, `force_callback` would still be called even when the option was not used.
    // Therefore, `check(CLI::ExistingFile)` could lead to an error because of an empty input path.
    throw error::Error("JSON file not found: {}", config_json);
  }
  std::ifstream config_stream(config_json);
  const auto config_json_data = Json::parse(config_stream);
  SVCConfigCollection config_collection;
  from_json(config_json_data, config_collection);
  auto& configs = config_collection.config_profiles;
  if (configs.empty()) {
    throw error::Error("No workflows found in config JSON file: {}", config_json);
  }
  const std::string profile_name{enum_util::FormatEnumName(profile)};
  if (!configs.contains(profile_name)) {
    throw error::Error("Specified config profile '{}' not found in config JSON file: {}", profile_name, config_json);
  }
  const auto expected_workflow = WorkflowForProfile(profile);
  const auto& profile_json = config_json_data.at("config_profiles").at(profile_name);
  // Inspect the raw JSON value before enum deserialization because nlohmann's enum conversion falls back to the first
  // enumerator for an unknown string, which would otherwise hide misspelled workflow names.
  if (!profile_json.contains("workflow") || !profile_json.at("workflow").is_string()) {
    throw error::Error("Config profile '{}' must declare a string workflow", profile_name);
  }
  const auto declared_workflow = profile_json.at("workflow").get<std::string>();
  const auto expected_workflow_name = Json(expected_workflow).get<std::string>();
  if (declared_workflow != expected_workflow_name) {
    throw error::Error("Config profile '{}' declares workflow '{}'; expected '{}'",
                       profile_name,
                       declared_workflow,
                       expected_workflow_name);
  }
  auto config = configs.at(profile_name);
  NormalizeGermlineModelTrainingDefaults(config);
  config.ConfigureFeatureCols();
  return config;
}
}  // namespace

SVCConfig JsonToConfig(const fs::path& config_json, const ConfigProfile profile) {
  if (config_json.empty()) {
    return SVCConfig(profile);
  }
  return JsonToConfigHelper(config_json, profile);
}

SVCConfig JsonToConfig(const fs::path& config_json, const Workflow workflow) {
  return JsonToConfig(config_json, ConfigProfileForWorkflow(workflow));
}

}  // namespace xoos::svc
