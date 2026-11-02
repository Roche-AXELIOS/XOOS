#pragma once

#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include <xoos/enum/enum-util.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>
#include <xoos/yc-decode/yc-decoder.h>

#include "config-profile.h"
#include "duplex-lowbq-mode.h"
#include "feature-normalization.h"
#include "germline-model-type.h"
#include "homopolymer-filter.h"
#include "sequencing-protocol.h"
#include "variant-info.h"
#include "workflow.h"
#include "yc-decode-method.h"

namespace xoos::svc {

using BaseType = xoos::yc_decode::BaseType;

/**
 * Synopsis:
 * This file contains default configurations for various workflows in the SVC submodule.
 * If no config JSON file is provided (via `--config`) to an SVC submodule, then the default configurations
 * (e.g. feature column names, CLI option thresholds, etc.) for the specified workflow (via `--workflow`) are taken
 * from this file.
 * All workflow-specific default configurations in this file must match those in `resources/profiles_config.json`.
 * A unit-test was designed to ensure no mismatch between the two files.
 */

// Base quality values are either 5,22,39 (or 0,18,93 for legacy scale)
// Threshold of 6 filters discordant bases in both scales
constexpr u8 kSimplexMinBaseQuality{6};
// Threshold of 0 includes all bases regardless of base quality
constexpr u8 kMinBaseQuality{0};
// Threshold of 23 filters non-concordant bases; used by tumor-normal-wgs, tumor-only-wgs, and germline-tagging
constexpr u8 kConcordantMinBaseQuality{23};

// Family size for duplex reads. By definition, a duplex read must consist of 2 constituent reads.
constexpr u32 kDuplexReadFamilySize{2};

// Family size for simplex reads. By definition, a simplex read is a single read, hence the family size is 1.
constexpr u32 kSimplexReadFamilySize{1};

/**
 * @brief Prepend the tumor sample prefix to a feature name.
 * @param feature_name Feature name to prefix
 * @return The feature name with the tumor prefix prepended
 */
std::string AddTumorPrefix(std::string_view feature_name);

/**
 * @brief Prepend the normal sample prefix to a feature name.
 * @param feature_name Feature name to prefix
 * @return The feature name with the normal prefix prepended
 */
std::string AddNormalPrefix(std::string_view feature_name);

/**
 * @brief Build a vector of strings from a list of string-like arguments.
 * @tparam Args Argument types convertible to `std::string` (e.g. `std::string_view`, `std::string`, `const char*`)
 * @param args Values to copy into the vector
 * @return Vector of strings
 * @note Used to initialize the runtime-mutable feature-name config aggregates from `constexpr std::string_view`
 * constants, since `std::string` is not implicitly constructible from `std::string_view` in a braced initializer
 * list. Accepting a variadic pack also allows mixing `std::string_view` constants with `std::string` values
 * (e.g. prefixed names from `AddTumorPrefix`).
 */
template <typename... Args>
vec<std::string> ToStringVec(Args&&... args) {
  return {std::string(std::forward<Args>(args))...};
}

// Default BAM features to be computed
vec<std::string> DefaultBamFeatures();

// Default VCF features to be computed
vec<std::string> DefaultVcfFeatures();

static constexpr std::vector<std::string> kDefaultScoringNames{};

static constexpr std::vector<std::string> kDefaultCategoricalNames{};

// BAM features for the `tumor-only-te` workflow
vec<std::string> TumorOnlyTeBamFeatures();

// VCF features for the `tumor-only-te` workflow
static constexpr std::vector<std::string> kTumorOnlyTeVcfFeatures{};

// Model scoring BAM features for the `tumor-only-te` workflow
vec<std::string> TumorOnlyTeScoringNames();

// Model scoring categorical features for the `tumor-only-te` workflow
vec<std::string> TumorOnlyTeCategoricalNames();

/** @brief Return scoring feature names for the `germline-wgs-duplex` indel model. */
vec<std::string> GermlineWgsDuplexIndelScoringNames();

/** @brief Return categorical scoring feature names for the `germline-wgs-duplex` SNV model. */
vec<std::string> GermlineWgsDuplexSnvCategoricalScoringNames();

/** @brief Return categorical scoring feature names for the `germline-wgs-duplex` indel model. */
vec<std::string> GermlineWgsDuplexIndelCategoricalScoringNames();

/** @brief Return BAM feature names for the `germline-wgs-duplex` profile. */
vec<std::string> GermlineWgsDuplexBamFeatures();

/** @brief Return VCF feature names for the `germline-wgs-duplex` profile. */
vec<std::string> GermlineWgsDuplexVcfFeatures();

/** @brief Return scoring feature names for the `germline-wgs-duplex` SNV model. */
vec<std::string> GermlineWgsDuplexSnvScoringNames();

/** @brief Return BAM feature names for the `germline-wgs-simplex` profile. */
vec<std::string> GermlineWgsSimplexBamFeatures();

/** @brief Return VCF feature names for the `germline-wgs-simplex` profile. */
vec<std::string> GermlineWgsSimplexVcfFeatures();

/** @brief Return scoring feature names for the `germline-wgs-simplex` SNV model. */
vec<std::string> GermlineWgsSimplexSnvScoringNames();

/** @brief Return scoring feature names for the `germline-wgs-simplex` indel model. */
vec<std::string> GermlineWgsSimplexIndelScoringNames();

/** @brief Return categorical scoring feature names for the `germline-wgs-simplex` SNV model. */
vec<std::string> GermlineWgsSimplexSnvCategoricalScoringNames();

/** @brief Return categorical scoring feature names for the `germline-wgs-simplex` indel model. */
vec<std::string> GermlineWgsSimplexIndelCategoricalScoringNames();

// BAM features for the `tumor-normal-wgs` workflow
vec<std::string> TumorNormalWgsBamFeatures();

// VCF features for the `tumor-normal-wgs` workflow
vec<std::string> TumorNormalWgsVcfFeatures();

// scoring feature names for the `tumor-normal-wgs` workflow
vec<std::string> TumorNormalWgsScoringNames();

// categorical scoring feature names for the `tumor-normal-wgs` workflow
vec<std::string> TumorNormalWgsCategoricalNames();

// BAM features for the `tumor-only-wgs` workflow
vec<std::string> TumorOnlyWgsBamFeatures();

// VCF features for the `tumor-only-wgs` workflow
vec<std::string> TumorOnlyWgsVcfFeatures();

// scoring feature names for the `tumor-only-wgs` workflow
vec<std::string> TumorOnlyWgsScoringNames();

// categorical scoring feature names for the `tumor-only-wgs` workflow
vec<std::string> TumorOnlyWgsCategoricalNames();

// BAM features for the `pon` workflow (minimal set for Panel of Normals generation)
vec<std::string> PanelOfNormalsBamFeatures();

vec<std::string> GermlineTaggingBamFeatures();

vec<std::string> GermlineTaggingVcfFeatures();

vec<std::string> GermlineTaggingScoringNames();

vec<std::string> GermlineTaggingCategoricalNames();

vec<std::string> GermlineMitochondriaBamFeatures();

vec<std::string> GermlineMitochondriaVcfFeatures();

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
 * @brief Optional filtering thresholds associated with one pre-trained model selection.
 * @details Entries are keyed by "{sample_type}-{aligner}" in profile JSON. Every field is optional so partial and
 * legacy configurations inherit omitted values from the profile's top-level thresholds.
 */
struct ModelThresholds {
  /** @brief Compare threshold collections field by field. */
  friend auto operator<=>(const ModelThresholds&, const ModelThresholds&) = default;

  // Model threshold entries are partial overrides. Missing fields in legacy or custom profile JSON must fall back to
  // the top-level tumor-normal-wgs config values instead of silently becoming 0.0.
  std::optional<f32> snv_min_ml_score;
  std::optional<f32> indel_min_ml_score;
  std::optional<f32> min_af;
  std::optional<f32> min_tumor_af;
};

/**
 * @brief Training defaults for one bundled germline WGS model family.
 * @details These values control training only. Feature extraction, feature schema selection, and filtering model
 * resolution are still selected by config profile, library type, and aligner.
 */
struct GermlineModelTrainingDefaults {
  /** @brief Compare training defaults field by field. */
  friend auto operator<=>(const GermlineModelTrainingDefaults&, const GermlineModelTrainingDefaults&) = default;

  u32 snv_iterations{};
  u32 indel_iterations{};
  std::string snv_model_lgbm_params;
  std::string indel_model_lgbm_params;
};

/**
 * @brief Serialize the threshold values that are explicitly present.
 * @details Explicit conversion preserves partial overrides. The lowercase name is required for nlohmann::json ADL.
 * @param json_obj Destination JSON object.
 * @param thresholds Threshold values to serialize.
 */
// NOLINTNEXTLINE(readability-identifier-naming)
inline void to_json(nlohmann::json& json_obj, const ModelThresholds& thresholds) {
  json_obj = nlohmann::json::object();
  if (thresholds.snv_min_ml_score.has_value()) {
    json_obj["snv_min_ml_score"] = *thresholds.snv_min_ml_score;
  }
  if (thresholds.indel_min_ml_score.has_value()) {
    json_obj["indel_min_ml_score"] = *thresholds.indel_min_ml_score;
  }
  if (thresholds.min_af.has_value()) {
    json_obj["min_af"] = *thresholds.min_af;
  }
  if (thresholds.min_tumor_af.has_value()) {
    json_obj["min_tumor_af"] = *thresholds.min_tumor_af;
  }
}

/**
 * @brief Deserialize a partial set of per-model thresholds.
 * @details Explicit conversion is used because the default nlohmann macros either reject missing fields or turn them
 * into numeric zero, preventing inheritance from top-level profile defaults.
 * @param json_obj Source JSON object.
 * @param thresholds Threshold collection to replace.
 */
// NOLINTNEXTLINE(readability-identifier-naming)
inline void from_json(const nlohmann::json& json_obj, ModelThresholds& thresholds) {
  thresholds = {};
  if (json_obj.contains("snv_min_ml_score") && !json_obj.at("snv_min_ml_score").is_null()) {
    thresholds.snv_min_ml_score = json_obj.at("snv_min_ml_score").get<f32>();
  }
  if (json_obj.contains("indel_min_ml_score") && !json_obj.at("indel_min_ml_score").is_null()) {
    thresholds.indel_min_ml_score = json_obj.at("indel_min_ml_score").get<f32>();
  }
  if (json_obj.contains("min_af") && !json_obj.at("min_af").is_null()) {
    thresholds.min_af = json_obj.at("min_af").get<f32>();
  }
  if (json_obj.contains("min_tumor_af") && !json_obj.at("min_tumor_af").is_null()) {
    thresholds.min_tumor_af = json_obj.at("min_tumor_af").get<f32>();
  }
}

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
   * @brief Create a SVCConfig from a concrete profile.
   * @param profile config profile enum
   * @throws xoos::error::Error if @p profile is not recognized.
   */
  explicit SVCConfig(ConfigProfile profile);

  /**
   * @brief Create a configuration for a workflow whose name uniquely identifies a profile.
   * @param target_workflow Workflow to configure.
   * @throws xoos::error::Error for germline WGS because its duplex and simplex profiles are ambiguous.
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
   * @brief Apply training defaults for one bundled germline WGS model family.
   * @details `train_model germline-wgs` uses this after selecting the concrete library-type profile and before CLI
   * iteration overrides are applied. The method intentionally mutates the top-level training fields so the rest of the
   * training code continues to consume one active configuration.
   * @param model_type Germline WGS model family whose training defaults should become active.
   * @throws xoos::error::Error if the selected model family is not configured.
   */
  void ApplyGermlineModelTrainingDefaults(GermlineModelType model_type);

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
   * @brief Configure parameters for the duplex "germline-wgs" profile.
   */
  void ConfigureGermlineWgsDuplexProfile();

  /**
   * @brief Configure parameters for the simplex "germline-wgs" profile.
   */
  void ConfigureGermlineWgsSimplexProfile();

  /**
   * @brief Configure parameters for the "tumor-only-te" workflow.
   */
  void ConfigureTumorOnlyTeWorkflow();

  /**
   * @brief Configure parameters for the "tumor-only-wgs" workflow.
   */
  void ConfigureTumorOnlyWgsWorkflow();

  /**
   * @brief Configure parameters for the "tumor-normal-wgs" workflow.
   */
  void ConfigureTumorNormalWgsWorkflow();

  /**
   * @brief Configure parameters for the "germline-tagging" workflow.
   */
  void ConfigureGermlineTaggingWorkflow();

  /**
   * @brief Configure parameters for the "germline-mitochondria" workflow.
   */
  void ConfigureMitochondriaWorkflow();

  /**
   * @brief Configure parameters for the "pon" workflow (Panel of Normals feature extraction).
   */
  void ConfigurePanelOfNormalsWorkflow();

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
  Workflow workflow = Workflow::kGermlineWgs;
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
  u32 min_dp{};
  f32 min_hom_alt_af{};
  f32 min_phased_af{};
  f32 max_phased_af{};
  f32 min_weighted_counts{};
  f32 hotspot_min_weighted_counts{};
  f32 min_ml_score{};
  f32 hotspot_min_ml_score{};
  f32 snv_min_ml_score{};
  f32 indel_min_ml_score{};
  u32 min_tumor_support{};
  u32 max_normal_support{};
  // Minimum alt support for a single-sample variant to be eligible for training/filtering (tumor-only-wgs).
  u32 min_support{};
  f32 min_tumor_af{};
  u32 max_indel_size{};
  f32 min_dp_ratio{};
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
  DuplexLowbqMode duplex_lowbq{DuplexLowbqMode::kInclude};

  // Per-model defaults keyed by "{sample_type}-{aligner}" (e.g. "ffpe-bwa").
  // Used by tumor-normal-wgs and tumor-only-wgs to tie thresholds to the auto-resolved model.
  // When the key is absent or the map is empty, the top-level config values are used.
  StrMap<ModelThresholds> model_thresholds;

  // Germline WGS training defaults keyed by model type ("large" or "small"). The top-level training fields remain the
  // active values so code that does not expose model-type selection keeps existing behavior.
  StrMap<GermlineModelTrainingDefaults> germline_model_type_training_defaults;
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
 * @brief Extract and set up a concrete profile from a JSON file. If an empty JSON path is provided,
 * then the built-in defaults are used for the specified profile.
 * @param config_json Path of JSON file
 * @param profile Config profile
 * @return Config for the profile
 */
SVCConfig JsonToConfig(const fs::path& config_json, ConfigProfile profile);

/**
 * @brief Extract and set up the config for a workflow that maps to one concrete profile.
 * @details If an empty JSON path is provided, then built-in defaults are used for the selected workflow.
 * Germline WGS is intentionally rejected because its workflow maps to multiple profiles and requires a resolved
 * library type.
 * @param config_json Path of JSON file.
 * @param workflow Workflow enum.
 * @return Config for the workflow.
 * @throws xoos::error::Error if @p workflow does not uniquely identify one configuration profile.
 */
SVCConfig JsonToConfig(const fs::path& config_json, Workflow workflow);

// macro to (de)serialize the Workflow enum to/from JSON
NLOHMANN_JSON_SERIALIZE_ENUM(Workflow,
                             {{Workflow::kGermlineWgs, "germline-wgs"},
                              {Workflow::kTumorOnlyTe, "tumor-only-te"},
                              {Workflow::kTumorOnlyWgs, "tumor-only-wgs"},
                              {Workflow::kTumorNormalWgs, "tumor-normal-wgs"},
                              {Workflow::kGermlineTagging, "germline-tagging"},
                              {Workflow::kGermlineMitochondria, "germline-mitochondria"},
                              {Workflow::kPanelOfNormals, "pon"},
                              {Workflow::kCustom, "custom"}})

// macro to (de)serialize the SequencingProtocol enum to/from JSON
NLOHMANN_JSON_SERIALIZE_ENUM(SequencingProtocol,
                             {{SequencingProtocol::kDuplex, "duplex"},
                              {SequencingProtocol::kDuplexSimplex, "duplex-simplex"},
                              {SequencingProtocol::kSimplex, "simplex"},
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

// JSON serialization for DuplexLowbqMode
NLOHMANN_JSON_SERIALIZE_ENUM(DuplexLowbqMode,
                             {{DuplexLowbqMode::kInclude, "include"}, {DuplexLowbqMode::kExclude, "exclude"}})

// JSON serialization for germline WGS model-family training defaults.
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(
    GermlineModelTrainingDefaults, snv_iterations, indel_iterations, snv_model_lgbm_params, indel_model_lgbm_params)

// The following macro allows us to read in an SVCConfig directly from JSON.
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(SVCConfig,
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
                                                duplex_lowbq,
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
                                                min_dp,
                                                min_hom_alt_af,
                                                min_phased_af,
                                                max_phased_af,
                                                min_weighted_counts,
                                                hotspot_min_weighted_counts,
                                                min_ml_score,
                                                hotspot_min_ml_score,
                                                phased,
                                                snv_min_ml_score,
                                                indel_min_ml_score,
                                                min_tumor_support,
                                                max_normal_support,
                                                min_support,
                                                min_tumor_af,
                                                max_indel_size,
                                                min_dp_ratio,
                                                model_thresholds,
                                                germline_model_type_training_defaults)

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
