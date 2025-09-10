#pragma once

#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include <xoos/enum/enum-util.h>
#include <xoos/error/error.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>

#include "core/column-names.h"
#include "core/variant-info.h"
#include "core/workflow.h"

using Json = nlohmann::json;

namespace xoos::svc {

/**
 * Synopsis:
 * This file contains default configurations for various workflows in the SVC submodule.
 * If no config JSON file is provided (via `--config`) to an SVC submodule, then the default configurations
 * (e.g. feature column names, CLI option thresholds, etc.) for the specified workflow (via `--workflow`) are taken
 * from this file.
 * All workflow-specific default configurations in this file must match those in `resources/profiles_config.json`.
 * A unit-test was designed to ensure no mismatch between the two files.
 */

// TODO: replace alternate names with their primary names and update JSON config files and model resources

// Family size for duplex reads. By definition, a duplex read must consist of 2 constituent reads.
constexpr u32 kDuplexReadFamilySize{2};

// Family size for simplex reads. By definition, a simplex read is a single read, hence the family size is 1.
constexpr u32 kSimplexReadFamilySize{1};

// Default BAM features to be computed; this is identical to the `somatic` workflow's BAM features
// TODO: derive a more generalized list that is not specific to the `somatic` workflow
static const std::vector<std::string> kDefaultFeatures{
    kNameChrom,         kNamePos,          kNameRef,          kNameAlt,           kNameMapqSum,
    kNameBaseqSum,      kNameNonDuplex,    kNameDuplex,       kNameFamilysizeSum, kNameDistanceSum,
    kNamePlusOnly,      kNameMinusOnly,    kNameRefSupport,   kNameSupport,       kNameContext,
    kNameWeightedScore, kNameBaseqMean,    kNameDistanceMean, kNameMapqMean,      kNameFamilysizeMean,
    kNameStrandBias,    kNameSubtypeIndex, kNameContextIndex};

// Default model scoring BAM features
// Also used as the model scoring BAM features for the `somatic` workflow
// TODO: create a separate variable for the `somatic` workflow
static const std::vector<std::string> kDefaultScoringCols{kNameWeightedScore,
                                                          kNameBaseqMean,
                                                          kNameDistanceMean,
                                                          kNameMapqMean,
                                                          kNameDuplex,
                                                          kNameSubtypeIndex,
                                                          kNameContextIndex,
                                                          kNameStrandBias,
                                                          kNamePlusOnly,
                                                          kNameMinusOnly,
                                                          kNameFamilysizeSum};

// Default model scoring categorical features
// Also used as the model scoring categorical features for the `somatic` workflow
// TODO: create a separate variable for the `somatic` workflow
static const std::vector<std::string> kDefaultCategoricalCols{kNameSubtypeIndex, kNameContextIndex};

// Default VCF features to be computed
static const std::vector<std::string> kDefaultVcfFeatures{kNameChrom,
                                                          kNamePos,
                                                          kNameRef,
                                                          kNameAlt,
                                                          kNameNalod,
                                                          kNameNlod,
                                                          kNameTlod,
                                                          kNameMpos,
                                                          kNameMmqRef,
                                                          kNameMmqAlt,
                                                          kNameMbqRef,
                                                          kNameMbqAlt,
                                                          kNameSubtypeIndex,
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
                                                          kNameTumorAltAD,
                                                          kNameNormalAltAD,
                                                          kNameTumorAF,
                                                          kNameNormalAF,
                                                          kNameTumorNormalAfRatio,
                                                          kNameTumorDP,
                                                          kNameNormalDP,
                                                          kNamePopAF,
                                                          kNameHapcomp,
                                                          kNameHapdom,
                                                          kNameRu,
                                                          kNameRpaRef,
                                                          kNameRpaAlt,
                                                          kNameStr,
                                                          kNameAtInterest};

// BAM features for the `unified` workflow
static const std::vector<std::string> kUnifiedFeatures{kNameChrom,
                                                       kNamePos,
                                                       kNameRef,
                                                       kNameAlt,
                                                       kNameWeightedDepth,
                                                       kNameRefWeightedDepth,
                                                       kNameRefNonHpWeightedDepth,
                                                       kNameSupport,
                                                       kNameRefNonHpSupport,
                                                       kNameRefSupport,
                                                       kNameMapqLT60Ratio,
                                                       kNameMapqLT40Ratio,
                                                       kNameMapqLT30Ratio,
                                                       kNameMapqLT20Ratio,
                                                       kNameRefNonHpMapqLT60Ratio,
                                                       kNameRefNonHpMapqLT40Ratio,
                                                       kNameRefNonHpMapqLT30Ratio,
                                                       kNameRefNonHpMapqLT20Ratio,
                                                       kNameRefMapqLT60Ratio,
                                                       kNameRefMapqLT40Ratio,
                                                       kNameRefMapqLT30Ratio,
                                                       kNameRefMapqLT20Ratio,
                                                       kNameMapqMin,
                                                       kNameMapqMax,
                                                       kNameMapqSum,
                                                       kNameMapqMean,
                                                       kNameBaseqMin,
                                                       kNameBaseqMax,
                                                       kNameBaseqSum,
                                                       kNameBaseqMean,
                                                       kNameDistanceMin,
                                                       kNameDistanceMax,
                                                       kNameDistanceSum,
                                                       kNameDistanceMean,
                                                       kNameRefNonHpBaseqMin,
                                                       kNameRefNonHpBaseqMax,
                                                       kNameRefNonHpBaseqSum,
                                                       kNameRefNonHpBaseqMean,
                                                       kNameRefNonHpMapqMin,
                                                       kNameRefNonHpMapqMax,
                                                       kNameRefNonHpMapqSum,
                                                       kNameRefNonHpMapqMean,
                                                       kNameRefDistanceMin,
                                                       kNameRefDistanceMax,
                                                       kNameRefDistanceSum,
                                                       kNameRefDistanceMean,
                                                       kNameRefMapqMin,
                                                       kNameRefMapqMax,
                                                       kNameRefMapqSum,
                                                       kNameRefMapqMean,
                                                       kNameMapqAF,
                                                       kNameBaseqAF,
                                                       kNameRefMapqAF,
                                                       kNameRefBaseqAF,
                                                       kNameRefBaseqMean,
                                                       kNameFamilysizeMean,
                                                       kNameFamilysizeLT3Ratio,
                                                       kNameFamilysizeLT5Ratio,
                                                       kNameRefFamilysizeMean,
                                                       kNameRefFamilysizeLT3Ratio,
                                                       kNameRefFamilysizeLT5Ratio,
                                                       kNameBaseqLT20Ratio,
                                                       kNameRefBaseqLT20Ratio,
                                                       kNameNonDuplex,
                                                       kNameDuplex,
                                                       kNameFamilysizeSum,
                                                       kNamePlusOnly,
                                                       kNameMinusOnly,
                                                       kNameContext,
                                                       kNameTumorSupport,
                                                       kNameTumorMapqMean,
                                                       kNameTumorBaseqMean,
                                                       kNameTumorDistanceMean,
                                                       kNameTumorRefSupport,
                                                       kNameNormalSupport,
                                                       kNameNormalMapqMean,
                                                       kNameNormalBaseqMean,
                                                       kNameNormalDistanceMean,
                                                       kNameNormalRefSupport,
                                                       kNameSupportReverse,
                                                       kNameTumorSupportReverse,
                                                       kNameNormalSupportReverse,
                                                       kNameBamTumorAF,
                                                       kNameBamNormalAF,
                                                       kNameBamRAT,
                                                       kNameAlignmentBias,
                                                       kNameTumorAlignmentBias,
                                                       kNameNormalAlignmentBias,
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
                                                       kNameADT,
                                                       kNameADTL,
                                                       kNameIndelAF,
                                                       kNameSimplex,
                                                       kNameMapqSumSimplex,
                                                       kNameMapqMeanSimplex,
                                                       kNameDistanceSumSimplex,
                                                       kNameDistanceMeanSimplex,
                                                       kNameRefSimplex,
                                                       kNameRefMapqSumSimplex,
                                                       kNameRefMapqMeanSimplex,
                                                       kNameRefDistanceSumSimplex,
                                                       kNameRefDistanceMeanSimplex};

// scoring feature names for the `germline` workflow's indel model
static const std::vector<std::string> kIndelScoringCols{kNameRef,
                                                        kNameAlt,
                                                        kNameDuplex,
                                                        kNameRefNonHpSupport,
                                                        kNameAltLen,
                                                        kNameRefSupport,
                                                        kNameDistanceMin,
                                                        kNameDistanceMax,
                                                        kNameDistanceMean,
                                                        kNameDistanceSum,
                                                        kNameBaseqAF,
                                                        kNameMapqAF,
                                                        kNameRefBaseqAF,
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
                                                        kNameSimplex,
                                                        kNameMapqSumSimplex,
                                                        kNameMapqMeanSimplex,
                                                        kNameDistanceSumSimplex,
                                                        kNameDistanceMeanSimplex,
                                                        kNameRefSimplex,
                                                        kNameRefMapqSumSimplex,
                                                        kNameRefMapqMeanSimplex,
                                                        kNameRefDistanceSumSimplex,
                                                        kNameRefDistanceMeanSimplex,
                                                        kNameAtInterest};

// categorical scoring feature names for the `germline` workflow's SNV model
static const std::vector<std::string> kSnvCategoricalScoringCols{kNameRef, kNameAlt, kNameGenotype, kNameAtInterest};

// categorical scoring feature names for the `germline` workflow's indel model
static const std::vector<std::string> kIndelCategoricalScoringCols{
    kNameRef, kNameAlt, kNamePre2BpContext, kNamePost2BpContext, kNameGenotype, kNameRu, kNameStr, kNameAtInterest};

// BAM features for the `germline` workflow
static const std::vector<std::string> kGermlineBAMFeatures{kNameChrom,
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
                                                           kNameBaseqAF,
                                                           kNameMapqAF,
                                                           kNameRefBaseqAF,
                                                           kNameRefMapqAF,
                                                           kNameSupport,
                                                           kNameWeightedDepth,
                                                           kNameBaseqMean,
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
                                                           kNameIndelAF,
                                                           kNameSimplex,
                                                           kNameMapqSumSimplex,
                                                           kNameMapqMeanSimplex,
                                                           kNameDistanceSumSimplex,
                                                           kNameDistanceMeanSimplex,
                                                           kNameRefSimplex,
                                                           kNameRefMapqSumSimplex,
                                                           kNameRefMapqMeanSimplex,
                                                           kNameRefDistanceSumSimplex,
                                                           kNameRefDistanceMeanSimplex};

// no VCF features for the `somatic` workflow
static const std::vector<std::string> kSomaticVCFFeatures{};

// VCF features for the `germline` workflow
static const std::vector<std::string> kGermlineVCFFeatures{kNameChrom,
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
                                                           kNameAtInterest};

// scoring feature names for the `germline` workflow's SNV model
static const std::vector<std::string> kSnvScoringCols{kNameRef,
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
                                                      kNameSimplex,
                                                      kNameMapqSumSimplex,
                                                      kNameMapqMeanSimplex,
                                                      kNameDistanceSumSimplex,
                                                      kNameDistanceMeanSimplex,
                                                      kNameRefSimplex,
                                                      kNameRefMapqSumSimplex,
                                                      kNameRefMapqMeanSimplex,
                                                      kNameRefDistanceSumSimplex,
                                                      kNameRefDistanceMeanSimplex,
                                                      kNameAtInterest};

// BAM features for the `somatic` workflow
static const std::vector<std::string> kSomaticFeatures{
    kNameChrom,         kNamePos,          kNameRef,          kNameAlt,           kNameMapqSum,
    kNameBaseqSum,      kNameNonDuplex,    kNameDuplex,       kNameFamilysizeSum, kNameDistanceSum,
    kNamePlusOnly,      kNameMinusOnly,    kNameRefSupport,   kNameSupport,       kNameContext,
    kNameWeightedScore, kNameBaseqMean,    kNameDistanceMean, kNameMapqMean,      kNameFamilysizeMean,
    kNameStrandBias,    kNameSubtypeIndex, kNameContextIndex};

// BAM features for the `somatic-tumor-normal` workflow
static const vec<std::string> kSomaticTumorNormalFeatures{kNameChrom,
                                                          kNamePos,
                                                          kNameRef,
                                                          kNameAlt,
                                                          kNameTumorSupport,
                                                          kNameTumorRefSupport,
                                                          kNameTumorDistanceMean,
                                                          kNameTumorMapqMean,
                                                          kNameTumorBaseqMean,
                                                          kNameNormalSupport,
                                                          kNameNormalRefSupport,
                                                          kNameNormalDistanceMean,
                                                          kNameNormalMapqMean,
                                                          kNameNormalBaseqMean,
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
                                                          kNameBamTumorAF,
                                                          kNameBamNormalAF,
                                                          kNameBamRAT,
                                                          kNameTumorSupportReverse,
                                                          kNameTumorAlignmentBias,
                                                          kNameContext};

// VCF features for the `somatic-tumor-normal` workflow
static const vec<std::string> kSomaticTumorNormalVCFFeatures{kNameChrom,         kNamePos,
                                                             kNameRef,           kNameAlt,
                                                             kNameNalod,         kNameNlod,
                                                             kNameTlod,          kNameMpos,
                                                             kNamePopAF,         kNameHapcomp,
                                                             kNameHapdom,        kNameRu,
                                                             kNameRpaRef,        kNameRpaAlt,
                                                             kNameSubtypeIndex,  kNameVariantType,
                                                             kNameUniq3mers,     kNameUniq4mers,
                                                             kNameUniq5mers,     kNameUniq6mers,
                                                             kNamePre2BpContext, kNamePost2BpContext,
                                                             kNameHomopolymer,   kNameDirepeat,
                                                             kNameTrirepeat,     kNameQuadrepeat};

struct SVCConfig {
  auto operator<=>(const SVCConfig&) const = default;

  // These functions take use the lists of feature names read in from the config or set by the constructor and generate
  // the list of matching enums to be used internally by SVC for feature computation, model training and filtering
  void GetFeatureCols();
  void GetVcfFeatureCols();
  void GetScoringCols();

  // These two functions validate the paths passed in from the command line for germline models.
  void GetGermlineModelPaths(const std::vector<fs::path>& paths);

#ifdef SOMATIC_ENABLE
  void GetSomaticTNModelPaths(const std::vector<fs::path>& paths);
#endif  // SOMATIC_ENABLE

  // This function is a wrapper to GetFeatureCols, GetVcfFeatureCols and GetScoringCols.
  void SetUpWorkflow();

  SVCConfig() = default;

  bool HasVcfFeatureScoringCols() const;

  /**
   * Create a SVCConfig from a preset workflow
   * @param target_workflow target workflow enum
   */
  explicit SVCConfig(const Workflow target_workflow) {
    switch (target_workflow) {
#ifdef SOMATIC_ENABLE
      case Workflow::kUnified: {
        feature_names = kUnifiedFeatures;
        vcf_feature_names = kDefaultVcfFeatures;
        workflow = Workflow::kUnified;
        scoring_names = kDefaultScoringCols;
        categorical_names = kDefaultCategoricalCols;
        n_classes = 4;
        min_mapq = 0;
        min_bq = 0;
        min_allowed_distance_from_end = 0;
        min_family_size = 0;
        filter_homopolymer = false;
        min_homopolymer_length = 0;
        duplex = true;
        min_alt_counts = 0;
        min_ctdna_allele_freq_threshold = 0;
        min_ffpe_allele_freq_threshold = 0;
        min_phased_allele_freq = 0;
        max_phased_allele_freq = 0;
        weighted_counts_threshold = 0;
        hotspot_weighted_counts_threshold = 0;
        ml_threshold = 0;
        hotspot_ml_threshold = 0;
        dynamic_thresholds = false;
        phased = false;
        use_vcf_features = true;
        iterations = 0;
        snv_iterations = 0;
        indel_iterations = 0;
        normalize_features = false;
        somatic_tn_snv_ml_threshold = 0.3;
        somatic_tn_indel_ml_threshold = 0.3;
        tumor_support_threshold = 1;
        decode_yc = false;
        model_lgbm_params =
            "objective=multiclass boosting=gbdt metric=multi_logloss seed=1238845 learning_rate=0.01 "
            "bagging_fraction=0.9 bagging_freq=1 feature_fraction=0.5 num_leaves=64 num_classes=4";
        break;
      }
#endif  // SOMATIC_ENABLE
      case Workflow::kGermline: {
        feature_names = kGermlineBAMFeatures;
        vcf_feature_names = kGermlineVCFFeatures;
        workflow = Workflow::kGermline;
        snv_scoring_names = kSnvScoringCols;
        indel_scoring_names = kIndelScoringCols;
        indel_categorical_names = kIndelCategoricalScoringCols;
        snv_categorical_names = kSnvCategoricalScoringCols;
        snv_model_file = "snv_model.txt";
        indel_model_file = "indel_model.txt";
        n_classes = 4;
        min_mapq = 1;
        min_bq = 6;  // filters disconcordant base in both 0,18,93 or 5,22,39
        min_allowed_distance_from_end = 2;
        min_family_size = kSimplexReadFamilySize;
        filter_homopolymer = true;
        min_homopolymer_length = 4;
        duplex = true;
        min_alt_counts = 0;
        min_ctdna_allele_freq_threshold = 0.0;
        min_ffpe_allele_freq_threshold = 0.0;
        min_phased_allele_freq = 0.0;
        max_phased_allele_freq = 0.0;
        weighted_counts_threshold = 0;
        hotspot_weighted_counts_threshold = 0;
        ml_threshold = 0.0;
        hotspot_ml_threshold = 0.0;
        dynamic_thresholds = false;
        phased = false;
        use_vcf_features = true;
        iterations = 0;
        snv_iterations = 1500;
        indel_iterations = 1500;
        normalize_features = true;
        somatic_tn_snv_ml_threshold = 0;
        somatic_tn_indel_ml_threshold = 0;
        tumor_support_threshold = 0;
        decode_yc = false;
        snv_model_lgbm_params =
            "objective=multiclass boosting=gbdt metric=multi_logloss seed=1238845 learning_rate=0.01 "
            "bagging_fraction=0.9 bagging_freq=1 feature_fraction=0.5 num_leaves=64 num_classes=4 ";
        indel_model_lgbm_params =
            "objective=multiclass boosting=gbdt metric=multi_logloss seed=1238845 learning_rate=0.02 "
            "bagging_fraction=0.9 bagging_freq=1 num_leaves=64 feature_fraction=0.5 num_classes=4 ";
        break;
      }
      case Workflow::kGermlineMultiSample: {
        feature_names = kGermlineBAMFeatures;
        vcf_feature_names = kGermlineVCFFeatures;
        workflow = Workflow::kGermlineMultiSample;
        snv_scoring_names = kSnvScoringCols;
        indel_scoring_names = kIndelScoringCols;
        indel_categorical_names = kIndelCategoricalScoringCols;
        snv_categorical_names = kSnvCategoricalScoringCols;
        snv_model_file = "snv_model.txt";
        indel_model_file = "indel_model.txt";
        n_classes = 4;
        min_mapq = 1;
        min_bq = 6;  // filters disconcordant base in both 0,18,93 or 5,22,39
        min_allowed_distance_from_end = 2;
        min_family_size = kSimplexReadFamilySize;
        filter_homopolymer = true;
        min_homopolymer_length = 4;
        duplex = true;
        min_alt_counts = 0;
        min_ctdna_allele_freq_threshold = 0.0;
        min_ffpe_allele_freq_threshold = 0.0;
        min_phased_allele_freq = 0.0;
        max_phased_allele_freq = 0.0;
        weighted_counts_threshold = 0;
        hotspot_weighted_counts_threshold = 0;
        ml_threshold = 0.0;
        hotspot_ml_threshold = 0.0;
        dynamic_thresholds = false;
        phased = false;
        use_vcf_features = true;
        iterations = 0;
        snv_iterations = 1000;
        indel_iterations = 8000;
        normalize_features = true;
        somatic_tn_snv_ml_threshold = 0;
        somatic_tn_indel_ml_threshold = 0;
        tumor_support_threshold = 0;
        decode_yc = false;
        snv_model_lgbm_params =
            "objective=multiclass boosting=gbdt metric=multi_logloss seed=1238845 learning_rate=0.0286803679209628 "
            "bagging_fraction=0.683997122877158 bagging_freq=1 feature_fraction=0.753244391002605 num_leaves=58 "
            "min_data_in_leaf=435 num_classes=4 ";
        indel_model_lgbm_params =
            "objective=multiclass boosting=gbdt metric=multi_logloss seed=1238845 learning_rate=0.0420473270669401 "
            "bagging_fraction=0.88062343020432 bagging_freq=1 num_leaves=59 feature_fraction=0.902813236067906 "
            "min_data_in_leaf=1325 num_classes=4 ";
        break;
      }
#ifdef SOMATIC_ENABLE
      case Workflow::kSomatic: {
        feature_names = kSomaticFeatures;
        vcf_feature_names = kSomaticVCFFeatures;
        workflow = Workflow::kSomatic;
        scoring_names = kDefaultScoringCols;
        categorical_names = kDefaultCategoricalCols;
        n_classes = 1;
        min_mapq = 9;
        min_bq = 18;
        min_allowed_distance_from_end = 0;
        min_family_size = 3;
        filter_homopolymer = false;
        min_homopolymer_length = 7;
        duplex = false;
        min_alt_counts = 3;
        min_ctdna_allele_freq_threshold = 0.0;
        min_ffpe_allele_freq_threshold = 0.01;
        min_phased_allele_freq = 0.001;
        max_phased_allele_freq = 0.5;
        weighted_counts_threshold = 4;
        hotspot_weighted_counts_threshold = 2;
        ml_threshold = 0.3;
        hotspot_ml_threshold = 0.01;
        dynamic_thresholds = false;
        phased = false;
        use_vcf_features = false;
        iterations = 1000;
        snv_iterations = 0;
        indel_iterations = 0;
        somatic_tn_snv_ml_threshold = 0.3;
        somatic_tn_indel_ml_threshold = 0.3;
        tumor_support_threshold = 1;
        decode_yc = false;
        model_lgbm_params = "objective=binary boosting=gbdt learning_rate=0.01 seed=1238845 num_leaves=3 num_classes=1";
        break;
      }
      case Workflow::kSomaticTumorNormal: {
        // TODO : add scoring names, categorical names
        // TODO : add value for `--max-variants-per-read`
        feature_names = kSomaticTumorNormalFeatures;
        vcf_feature_names = kSomaticTumorNormalVCFFeatures;
        workflow = Workflow::kSomaticTumorNormal;
        n_classes = 1;
        min_mapq = 1;
        min_bq = 19;
        min_allowed_distance_from_end = 2;
        min_family_size = kDuplexReadFamilySize;
        filter_homopolymer = false;
        min_homopolymer_length = 7;
        duplex = true;
        min_alt_counts = 3;
        min_ctdna_allele_freq_threshold = 0.0;
        min_ffpe_allele_freq_threshold = 0.01;
        min_phased_allele_freq = 0.001;
        max_phased_allele_freq = 0.5;
        weighted_counts_threshold = 4;
        hotspot_weighted_counts_threshold = 2;
        ml_threshold = 0.3;
        hotspot_ml_threshold = 0.01;
        dynamic_thresholds = false;
        phased = false;
        use_vcf_features = true;
        iterations = 1000;
        snv_iterations = 0;
        indel_iterations = 0;
        normalize_features = false;
        somatic_tn_snv_ml_threshold = 0.3;    // TODO : Update value with appropriate default
        somatic_tn_indel_ml_threshold = 0.3;  // TODO : Update value with appropriate default
        tumor_support_threshold = 1;
        decode_yc = false;
        break;
      }
#endif  // SOMATIC_ENABLE
      default: {
        throw error::Error("Unknown Workflow value {}", enum_util::FormatEnumName(target_workflow));
      }
    }
    GetFeatureCols();
    GetScoringCols();
    GetVcfFeatureCols();
  }

  std::vector<std::string> feature_names;
  std::vector<UnifiedFeatureCols> feature_cols;
  std::vector<std::string> scoring_names;
  std::vector<UnifiedFeatureCols> scoring_cols;
  std::vector<std::string> snv_scoring_names;
  std::vector<UnifiedFeatureCols> snv_scoring_cols;
  std::vector<std::string> indel_scoring_names;
  std::vector<std::string> indel_categorical_names;
  std::vector<std::string> snv_categorical_names;
  std::vector<UnifiedFeatureCols> indel_scoring_cols;
  std::vector<std::string> vcf_feature_names;
  std::vector<UnifiedFeatureCols> vcf_feature_cols;
  std::vector<std::string> categorical_names;
  std::vector<std::string> somatic_tn_scoring_names;
  std::vector<UnifiedFeatureCols> somatic_tn_scoring_cols;
  std::vector<std::string> somatic_tn_categorical_names;
  std::vector<std::string> germline_fail_scoring_names;
  std::vector<std::string> germline_fail_categorical_names;
  std::vector<UnifiedFeatureCols> germline_fail_scoring_cols;
  Workflow workflow = Workflow::kGermline;
  size_t n_classes = 1;

  fs::path snv_model_file;
  fs::path indel_model_file;
  fs::path germline_fail_model_file;
  fs::path model_file;

  std::string snv_model_lgbm_params;
  std::string indel_model_lgbm_params;
  std::string model_lgbm_params;
  std::string germline_fail_lgbm_params;

  u8 min_mapq{};
  u8 min_bq{};
  float min_allowed_distance_from_end{};
  u32 min_family_size{};
  bool filter_homopolymer{};
  u32 min_homopolymer_length{};
  bool duplex{};
  u32 min_alt_counts{};
  float min_ctdna_allele_freq_threshold{};
  float min_ffpe_allele_freq_threshold{};
  float min_phased_allele_freq{};
  float max_phased_allele_freq{};
  float weighted_counts_threshold{};
  float hotspot_weighted_counts_threshold{};
  float ml_threshold{};
  float hotspot_ml_threshold{};
  float somatic_tn_snv_ml_threshold{};
  float somatic_tn_indel_ml_threshold{};
  u32 tumor_support_threshold{};
  bool dynamic_thresholds{};
  bool phased{};
  bool use_vcf_features{};
  u32 iterations{};
  u32 snv_iterations{};
  u32 indel_iterations{};
  bool normalize_features{};
  bool decode_yc{};
};

// A collection of SVCConfig structs indexed by their workflow names.
struct SVCConfigCollection {
  StrMap<SVCConfig> config_profiles;
};

SVCConfigCollection JsonToConfigCollection(const fs::path& config_json);
SVCConfig JsonToConfig(const fs::path& config_json, const std::string& workflow);
SVCConfig JsonToConfig(const fs::path& config_json, Workflow workflow);

// We use the macro below to serialize the workflow enum allowing us to parse the enum in from a JSON file
NLOHMANN_JSON_SERIALIZE_ENUM(Workflow,
                             {{Workflow::kGermline, "germline"},
                              {Workflow::kGermlineMultiSample, "germline-multi-sample"}})

// The following macro allows us to read in an SVCConfig directly from JSON.
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SVCConfig,
                                   feature_names,
                                   snv_scoring_names,
                                   indel_scoring_names,
                                   indel_categorical_names,
                                   snv_categorical_names,
                                   vcf_feature_names,
                                   n_classes,
                                   workflow,
                                   min_mapq,
                                   min_bq,
                                   min_allowed_distance_from_end,
                                   min_family_size,
                                   filter_homopolymer,
                                   min_homopolymer_length,
                                   duplex,
                                   use_vcf_features,
                                   snv_iterations,
                                   indel_iterations,
                                   snv_model_file,
                                   indel_model_file,
                                   normalize_features,
                                   decode_yc,
                                   snv_model_lgbm_params,
                                   indel_model_lgbm_params)

// The following macro allows us to parse a collection of SVCConfig(s) from one JSON file, and select the workflow to
// be executed. Each SVCConfig struct corresponds to the parameter settings for one specific workflow.
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SVCConfigCollection, config_profiles)

}  // namespace xoos::svc
