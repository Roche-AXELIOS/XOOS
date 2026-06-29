#pragma once

#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include <xoos/enum/enum-util.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>
#include <xoos/yc-decode/yc-decoder.h>

#include "core/column-names.h"
#include "core/feature-normalization.h"
#include "core/variant-info.h"
#include "core/workflow.h"
#include "homopolymer-filter.h"
#include "sequencing-protocol.h"
#include "xoos/types/float.h"
#include "yc-decode-method.h"

using Json = nlohmann::json;
using BaseType = xoos::yc_decode::BaseType;

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

// Base quality values are either 5,22,39 (or 0,18,93 for legacy scale)
// Threshold of 6 filters discordant bases in both scales
constexpr u8 kSimplexMinBaseQuality{6};
// Threshold of 23 filters discordant and simplex bases in both scales
constexpr u8 kConcordantMinBaseQuality{23};

// Family size for duplex reads. By definition, a duplex read must consist of 2 constituent reads.
constexpr u32 kDuplexReadFamilySize{2};

// Family size for simplex reads. By definition, a simplex read is a single read, hence the family size is 1.
constexpr u32 kSimplexReadFamilySize{1};

// Add tumor feature name prefix to a feature name
static std::string AddTumorPrefix(const std::string& feature_name) {
  return kTumorPrefix + feature_name;
}

// Add normal feature name prefix to a feature name
static std::string AddNormalPrefix(const std::string& feature_name) {
  return kNormalPrefix + feature_name;
}

// Default BAM features to be computed
static const std::vector<std::string> kDefaultBamFeatures{
    kNameChrom, kNamePos, kNameRef, kNameAlt, kNameDuplex, kNameDuplexLowbq, kNameDuplexDP, kNameDuplexAF};

// Default VCF features to be computed
static const std::vector<std::string> kDefaultVcfFeatures{kNameChrom,
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
                                                          kNamePopAF};

static constexpr std::vector<std::string> kDefaultScoringNames{};

static constexpr std::vector<std::string> kDefaultCategoricalNames{};

// BAM features for the `tumor-only-te` workflow
static const std::vector<std::string> kTumorOnlyTeBamFeatures{
    kNameChrom,         kNamePos,          kNameRef,          kNameAlt,           kNameMapqSum,
    kNameBaseqSum,      kNameNonDuplex,    kNameDuplex,       kNameFamilysizeSum, kNameDistanceSum,
    kNamePlusOnly,      kNameMinusOnly,    kNameRefSupport,   kNameSupport,       kNameContext,
    kNameWeightedScore, kNameBaseqMean,    kNameDistanceMean, kNameMapqMean,      kNameFamilysizeMean,
    kNameStrandBias,    kNameSubtypeIndex, kNameContextIndex};

// VCF features for the `tumor-only-te` workflow
static constexpr std::vector<std::string> kTumorOnlyTeVcfFeatures{};

// Model scoring BAM features for the `tumor-only-te` workflow
static const std::vector<std::string> kTumorOnlyTeScoringNames{kNameWeightedScore,
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

// Model scoring categorical features for the `tumor-only-te` workflow
static const std::vector<std::string> kTumorOnlyTeCategoricalNames{kNameSubtypeIndex, kNameContextIndex};

// scoring feature names for the `germline` workflow's indel model
static const std::vector<std::string> kGermlineIndelScoringNames{kNameRef,
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
static const std::vector<std::string> kGermlineSnvCategoricalScoringNames{
    kNameRef, kNameAlt, kNameGenotype, kNameAtInterest};

// categorical scoring feature names for the `germline` workflow's indel model
static const std::vector<std::string> kGermlineIndelCategoricalScoringNames{
    kNameRef, kNameAlt, kNamePre2BpContext, kNamePost2BpContext, kNameGenotype, kNameRu, kNameStr, kNameAtInterest};

// BAM features for the `germline` workflow
static const std::vector<std::string> kGermlineBamFeatures{kNameChrom,
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

// VCF features for the `germline` workflow
static const std::vector<std::string> kGermlineVcfFeatures{kNameChrom,
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
static const std::vector<std::string> kGermlineSnvScoringNames{kNameRef,
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

// BAM features for the `tumor-normal-wgs` workflow
static const vec<std::string> kTumorNormalWgsBamFeatures{kNameChrom,
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
                                                         kNameDistanceMean,
                                                         kNameRefDistanceMean,
                                                         AddTumorPrefix(kNameDuplexAF),
                                                         AddNormalPrefix(kNameDuplexAF),
                                                         kNameBamTnAfRatio,
                                                         AddTumorPrefix(kNameSupportReverse),
                                                         AddTumorPrefix(kNameAlignmentBias),
                                                         kNameContext};

// VCF features for the `tumor-normal-wgs` workflow
static const vec<std::string> kTumorNormalWgsVcfFeatures{
    kNameChrom,       kNamePos,       kNameRef,       kNameAlt,       kNameTumorDP,       kNameNormalDP,
    kNameNalod,       kNameNlod,      kNameTlod,      kNameMpos,      kNamePopAF,         kNameHapcomp,
    kNameHapdom,      kNameRu,        kNameRpaRef,    kNameRpaAlt,    kNameSubtypeIndex,  kNameVariantType,
    kNameUniq3mers,   kNameUniq4mers, kNameUniq5mers, kNameUniq6mers, kNamePre2BpContext, kNamePost2BpContext,
    kNameHomopolymer, kNameDirepeat,  kNameTrirepeat, kNameQuadrepeat};

// scoring feature names for the `tumor-normal-wgs` workflow
static const vec<std::string> kTumorNormalWgsScoringNames{kNameNalod,
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
                                                          kNameDistanceMean,
                                                          kNameRefDistanceMean,
                                                          AddTumorPrefix(kNameDuplexAF),
                                                          AddNormalPrefix(kNameDuplexAF),
                                                          kNameBamTnAfRatio,
                                                          AddTumorPrefix(kNameAlignmentBias),
                                                          kNameContext};

// categorical scoring feature names for the `tumor-normal-wgs` workflow
static const vec<std::string> kTumorNormalWgsCategoricalNames{
    kNameRu, kNameStr, kNameSubtypeIndex, kNameVariantType, kNamePre2BpContext, kNamePost2BpContext, kNameContext};

static const vec<std::string> kGermlineTaggingBamFeatures{kNameChrom,
                                                          kNamePos,
                                                          kNameRef,
                                                          kNameAlt,
                                                          AddNormalPrefix(kNameMapqMean),
                                                          AddNormalPrefix(kNameBaseqMean),
                                                          AddNormalPrefix(kNameDistanceMean),
                                                          AddNormalPrefix(kNameSupport),
                                                          AddNormalPrefix(kNameRefSupport),
                                                          AddNormalPrefix(kNameDuplexAF),
                                                          kNameContext};

static const vec<std::string> kGermlineTaggingVcfFeatures{kNameChrom,
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
                                                          kNameStr};

static const vec<std::string> kGermlineTaggingScoringNames{kNameRef,
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
                                                           kNameStr};

static const vec<std::string> kGermlineTaggingCategoricalNames{kNameRef,
                                                               kNameAlt,
                                                               kNameVariantType,
                                                               kNameSubtypeIndex,
                                                               kNamePre2BpContext,
                                                               kNamePost2BpContext,
                                                               kNameContext,
                                                               kNameRu,
                                                               kNameStr};

/**
 * @brief Struct to hold file paths containing "snv" and "indel" in their file names.
 */
struct SnvIndelPaths {
  fs::path snv_path;
  fs::path indel_path;
};

/**
 * @brief Extracts and returns the SNV and Indel file paths from the provided list of paths.
 * @param paths A vector of file paths to search through.
 * @return A SnvIndelPaths struct containing the SNV and Indel file paths
 */
SnvIndelPaths GetSnvIndelPaths(const std::vector<fs::path>& paths);

/**
 * @brief Configuration structure for the SVC submodule.
 * Contains all parameters needed to set up and run the SVC submodule.
 */
struct SVCConfig {
  /**
   * @brief Default three-way comparison operator
   */
  auto operator<=>(const SVCConfig&) const = default;

  /**
   * @brief Default constructor
   */
  SVCConfig() = default;

  /**
   * @brief Create a SVCConfig from a preset workflow.
   * @param target_workflow target workflow enum
   */
  explicit SVCConfig(Workflow target_workflow);

  /**
   * Populates a vector of UnifiedFeatureCol enums that match the requested column names for BAM features
   */
  void ConfigureBamFeatureCols();

  /**
   * Populates a vector of UnifiedFeatureCol enums that match the requested column names for VCF features
   */
  void ConfigureVcfFeatureCols();

  /**
   * @brief Populates a vector of UnifiedFeatureCol enums that match the requested column names for model scoring
   * columns
   */
  void ConfigureScoringCols();

  /**
   * @brief Wrapper helper to populate the BAM and VCF feature enum sets and model training and scoring feature sets
   * with enum values once a SVCConfig has been read in from file.
   */
  void ConfigureFeatureCols();

  /**
   * @brief Checks whether any of the scoring columns are (derivatives of) VCF features
   * @return true if any scoring column is a VCF feature, false otherwise
   */
  bool HasVcfFeatureScoringCols() const;

  /**
   * @brief Configure parameters for the "custom" workflow.
   */
  void ConfigureCustomWorkflow();

  /**
   * @brief Configure parameters for the "germline" workflow.
   */
  void ConfigureGermlineWorkflow();

  /**
   * @brief Configure parameters for the "germline-multi-sample" workflow.
   */
  void ConfigureGermlineMultiSampleWorkflow();

  /**
   * @brief Configure parameters for the "tumor-only-te" workflow.
   */
  void ConfigureTumorOnlyTeWorkflow();

  /**
   * @brief Configure parameters for the "tumor-normal-wgs" workflow.
   */
  void ConfigureTumorNormalWgsWorkflow();

  /**
   * @brief Configure parameters for the "germline-tagging" workflow.
   */
  void ConfigureGermlineTaggingWorkflow();

  std::vector<std::string> bam_feature_names;
  std::vector<FeatureColumn> feature_cols;
  std::vector<std::string> scoring_names;
  std::vector<FeatureColumn> scoring_cols;
  std::vector<std::string> snv_scoring_names;
  std::vector<FeatureColumn> snv_scoring_cols;
  std::vector<std::string> indel_scoring_names;
  std::vector<std::string> indel_categorical_names;
  std::vector<std::string> snv_categorical_names;
  std::vector<FeatureColumn> indel_scoring_cols;
  std::vector<std::string> vcf_feature_names;
  std::vector<FeatureColumn> vcf_feature_cols;
  std::vector<std::string> categorical_names;
  std::vector<std::string> somatic_tn_scoring_names;
  std::vector<FeatureColumn> somatic_tn_scoring_cols;
  std::vector<std::string> somatic_tn_categorical_names;
  std::vector<std::string> germline_fail_scoring_names;
  std::vector<std::string> germline_fail_categorical_names;
  std::vector<FeatureColumn> germline_fail_scoring_cols;
  Workflow workflow = Workflow::kGermline;
  size_t n_classes = 0;

  std::string snv_model_lgbm_params;
  std::string indel_model_lgbm_params;
  std::string model_lgbm_params;
  std::string germline_fail_lgbm_params;

  std::string snv_model_lgbm_prediction_params;
  std::string indel_model_lgbm_prediction_params;
  std::string model_lgbm_prediction_params;

  u8 min_mapq{};
  u8 min_bq{};
  f32 min_dist{};
  u32 min_family_size{};
  HomopolymerFilter filter_homopolymer{HomopolymerFilter::kNone};
  u32 min_homopolymer_length{};
  SequencingProtocol sequencing_protocol{SequencingProtocol::kDuplex};
  u32 min_alt_counts{};
  f32 min_af{};
  f32 min_phased_af{};
  f32 max_phased_af{};
  f32 min_weighted_counts{};
  f32 hotspot_min_weighted_counts{};
  f32 min_ml_score{};
  f32 hotspot_min_ml_score{};
  f32 snv_min_ml_score{};
  f32 indel_min_ml_score{};
  u32 min_tumor_support{};
  bool phased{};
  bool use_vcf_features{};
  u32 iterations{};
  u32 snv_iterations{};
  u32 indel_iterations{};
  FeatureNormalization normalize_features{FeatureNormalization::kNone};
  YcDecodeMethod decode_yc = YcDecodeMethod::kNone;
  yc_decode::BaseType min_base_type = BaseType::kDiscordant;
  u32 max_variants_per_read{};
  f32 max_variants_per_read_normalized{};
};

/**
 * @brief Struct to hold a collection of SVCConfig(s) indexed by their profile names.
 * This allows us to read in multiple profiles from one JSON file, and select the desired profile at runtime.
 */
struct SVCConfigCollection {
  StrMap<SVCConfig> config_profiles;
};

/**
 * @brief Extract a collection of configs from a JSON file.
 * @param config_json Path of JSON file
 * @return Collection of configs
 */
SVCConfigCollection JsonToConfigCollection(const fs::path& config_json);

/**
 * @brief Extract and set up the config for a given workflow from a JSON file. If an empty JSON path is provided,
 * then the default config values are used for the specified workflow.
 * @param config_json Path of JSON file
 * @param workflow Workflow name
 * @return Config for the workflow
 */
SVCConfig JsonToConfig(const fs::path& config_json, const std::string& workflow);

/**
 * @brief Extract and set up the config for a given workflow from a JSON file. If an empty JSON path is provided,
 * then the default config values are used for the specified workflow.
 * @param config_json Path of JSON file
 * @param workflow Workflow enum
 * @return Config for the workflow
 */
SVCConfig JsonToConfig(const fs::path& config_json, Workflow workflow);

// macro to (de)serialize the Workflow enum to/from JSON
NLOHMANN_JSON_SERIALIZE_ENUM(Workflow,
                             {{Workflow::kGermline, "germline"},
                              {Workflow::kGermlineMultiSample, "germline-multi-sample"},
                              {Workflow::kTumorOnlyTe, "tumor-only-te"},
                              {Workflow::kTumorNormalWgs, "tumor-normal-wgs"},
                              {Workflow::kGermlineTagging, "germline-tagging"},
                              {Workflow::kCustom, "custom"}})

// macro to (de)serialize the SequencingProtocol enum to/from JSON
NLOHMANN_JSON_SERIALIZE_ENUM(SequencingProtocol,
                             {{SequencingProtocol::kDuplex, "duplex"},
                              {SequencingProtocol::kDuplexSimplex, "duplex-simplex"},
                              {SequencingProtocol::kUmi, "umi"}})

// macro to (de)serialize the YcDecodeMethod enum to/from JSON
NLOHMANN_JSON_SERIALIZE_ENUM(YcDecodeMethod,
                             {{YcDecodeMethod::kNone, "none"},
                              {YcDecodeMethod::kConsensus, "consensus"},
                              {YcDecodeMethod::kSplit, "split"}})

// macro to (de)serialize the FeatureNormalizationProtocol enum to/from JSON
NLOHMANN_JSON_SERIALIZE_ENUM(FeatureNormalization,
                             {{FeatureNormalization::kNone, "none"}, {FeatureNormalization::kMedianDp, "median-dp"}})

// macro to (de)serialize the HomopolymerFilter enum to/from JSON
NLOHMANN_JSON_SERIALIZE_ENUM(HomopolymerFilter,
                             {{HomopolymerFilter::kNone, "none"}, {HomopolymerFilter::kAlignmentEnd, "alignment-end"}})

// The following macro allows us to read in an SVCConfig directly from JSON.
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SVCConfig,
                                   bam_feature_names,
                                   snv_scoring_names,
                                   indel_scoring_names,
                                   indel_categorical_names,
                                   snv_categorical_names,
                                   vcf_feature_names,
                                   n_classes,
                                   workflow,
                                   min_mapq,
                                   min_bq,
                                   min_dist,
                                   max_variants_per_read,
                                   max_variants_per_read_normalized,
                                   min_family_size,
                                   filter_homopolymer,
                                   min_homopolymer_length,
                                   sequencing_protocol,
                                   use_vcf_features,
                                   snv_iterations,
                                   indel_iterations,
                                   normalize_features,
                                   decode_yc,
                                   min_base_type,
                                   snv_model_lgbm_params,
                                   indel_model_lgbm_params,
                                   scoring_names,
                                   categorical_names,
                                   model_lgbm_params,
                                   snv_model_lgbm_prediction_params,
                                   indel_model_lgbm_prediction_params,
                                   model_lgbm_prediction_params,
                                   iterations,
                                   min_alt_counts,
                                   min_af,
                                   min_phased_af,
                                   max_phased_af,
                                   min_weighted_counts,
                                   hotspot_min_weighted_counts,
                                   min_ml_score,
                                   hotspot_min_ml_score,
                                   phased,
                                   snv_min_ml_score,
                                   indel_min_ml_score,
                                   min_tumor_support)

// The following macro allows us to parse a collection of SVCConfig(s) from one JSON file, and select the workflow to
// be executed. Each SVCConfig struct corresponds to the parameter settings for one specific workflow.
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SVCConfigCollection, config_profiles)

}  // namespace xoos::svc

namespace xoos::yc_decode {
// Macro to serialize/deserialize BaseType enum to/from JSON.
// This must be declared in the same namespace as the enum.
// Otherwise, nlohmann::json will not be able to locate it, and it will default to integer serialization.
NLOHMANN_JSON_SERIALIZE_ENUM(BaseType,
                             {{BaseType::kDiscordant, "discordant"},
                              {BaseType::kSimplex, "simplex"},
                              {BaseType::kConcordant, "concordant"}})
}  // namespace xoos::yc_decode
