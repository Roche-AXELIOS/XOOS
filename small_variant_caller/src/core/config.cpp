#include "config.h"

#include <fstream>

#include <xoos/error/error.h>
#include <xoos/util/string-functions.h>

#include "variant-info.h"

namespace xoos::svc {

SnvIndelPaths GetSnvIndelPaths(const std::vector<fs::path>& paths) {
  SnvIndelPaths out_paths{};
  // distinguish SNV and indel paths based on file names
  for (const auto& p : paths) {
    auto file_name = p.filename().string();
    const bool snv_found = file_name.find("snv") != std::string::npos;
    const bool indel_found = file_name.find("indel") != std::string::npos;
    if (snv_found && !indel_found) {
      out_paths.snv_path = p;
    } else if (indel_found && !snv_found) {
      out_paths.indel_path = p;
    }
  }
  return out_paths;
}

SVCConfig::SVCConfig(const Workflow target_workflow) {
  using enum Workflow;
  switch (target_workflow) {
    case kCustom: {
      ConfigureCustomWorkflow();
      break;
    }
    case kGermline: {
      ConfigureGermlineWorkflow();
      break;
    }
    case kGermlineMultiSample: {
      ConfigureGermlineMultiSampleWorkflow();
      break;
    }
    case kTumorOnlyTe: {
      ConfigureTumorOnlyTeWorkflow();
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
    default: {
      throw error::Error("Unknown Workflow value {}", enum_util::FormatEnumName(target_workflow));
    }
  }
  ConfigureFeatureCols();
}

void SVCConfig::ConfigureCustomWorkflow() {
  bam_feature_names = kDefaultBamFeatures;
  vcf_feature_names = kDefaultVcfFeatures;
  workflow = Workflow::kCustom;
  scoring_names = kDefaultScoringNames;
  categorical_names = kDefaultCategoricalNames;
  n_classes = 0;
  min_mapq = 1;
  min_bq = kSimplexMinBaseQuality;
  min_dist = 2;
  min_family_size = 1;
  filter_homopolymer = HomopolymerFilter::kAlignmentEnd;
  min_homopolymer_length = 4;
  sequencing_protocol = SequencingProtocol::kDuplex;
  use_vcf_features = true;
}

void SVCConfig::ConfigureGermlineWorkflow() {
  bam_feature_names = kGermlineBamFeatures;
  vcf_feature_names = kGermlineVcfFeatures;
  workflow = Workflow::kGermline;
  snv_scoring_names = kGermlineSnvScoringNames;
  indel_scoring_names = kGermlineIndelScoringNames;
  indel_categorical_names = kGermlineIndelCategoricalScoringNames;
  snv_categorical_names = kGermlineSnvCategoricalScoringNames;
  n_classes = 4;
  min_mapq = 1;
  min_bq = kSimplexMinBaseQuality;
  min_dist = 2;
  min_family_size = kSimplexReadFamilySize;
  filter_homopolymer = HomopolymerFilter::kAlignmentEnd;
  min_homopolymer_length = 4;
  sequencing_protocol = SequencingProtocol::kDuplex;
  use_vcf_features = true;
  snv_iterations = 1500;
  indel_iterations = 1500;
  normalize_features = FeatureNormalization::kMedianDp;
  decode_yc = YcDecodeMethod::kConsensus;
  min_base_type = BaseType::kSimplex;
  snv_model_lgbm_params =
      "objective=multiclass boosting=gbdt metric=multi_logloss seed=1238845 learning_rate=0.01 "
      "bagging_fraction=0.9 bagging_freq=1 feature_fraction=0.5 num_leaves=64 num_classes=4 ";
  indel_model_lgbm_params =
      "objective=multiclass boosting=gbdt metric=multi_logloss seed=1238845 learning_rate=0.02 "
      "bagging_fraction=0.9 bagging_freq=1 num_leaves=64 feature_fraction=0.5 num_classes=4 ";
  snv_model_lgbm_prediction_params =
      "num_threads=1 force_row_wise=true deterministic=true seed=1238845 pred_early_stop=true";
  indel_model_lgbm_prediction_params =
      "num_threads=1 force_row_wise=true deterministic=true seed=1238845 pred_early_stop=true";
}

void SVCConfig::ConfigureGermlineMultiSampleWorkflow() {
  ConfigureGermlineWorkflow();
  workflow = Workflow::kGermlineMultiSample;
  snv_iterations = 1000;
  indel_iterations = 8000;
  snv_model_lgbm_params =
      "objective=multiclass boosting=gbdt metric=multi_logloss seed=1238845 learning_rate=0.0286803679209628 "
      "bagging_fraction=0.683997122877158 bagging_freq=1 feature_fraction=0.753244391002605 num_leaves=58 "
      "min_data_in_leaf=435 num_classes=4 ";
  indel_model_lgbm_params =
      "objective=multiclass boosting=gbdt metric=multi_logloss seed=1238845 learning_rate=0.0420473270669401 "
      "bagging_fraction=0.88062343020432 bagging_freq=1 num_leaves=59 feature_fraction=0.902813236067906 "
      "min_data_in_leaf=1325 num_classes=4 ";
  snv_model_lgbm_prediction_params =
      "num_threads=1 force_row_wise=true deterministic=true seed=1238845 pred_early_stop=true";
  indel_model_lgbm_prediction_params =
      "num_threads=1 force_row_wise=true deterministic=true seed=1238845 pred_early_stop=true";
}

void SVCConfig::ConfigureTumorOnlyTeWorkflow() {
  bam_feature_names = kTumorOnlyTeBamFeatures;
  vcf_feature_names = kTumorOnlyTeVcfFeatures;
  workflow = Workflow::kTumorOnlyTe;
  scoring_names = kTumorOnlyTeScoringNames;
  categorical_names = kTumorOnlyTeCategoricalNames;
  n_classes = 1;
  min_mapq = 9;
  min_bq = kSimplexMinBaseQuality;
  min_dist = 0;
  max_variants_per_read = 10;
  min_family_size = 3;
  filter_homopolymer = HomopolymerFilter::kNone;
  min_homopolymer_length = 7;
  sequencing_protocol = SequencingProtocol::kUmi;
  min_alt_counts = 3;
  // min AF for FFPE: 0.01, for ctDNA: 0.0
  min_af = 0.01F;
  min_phased_af = 0.001F;
  max_phased_af = 0.5F;
  min_weighted_counts = 4;
  hotspot_min_weighted_counts = 2;
  min_ml_score = 0.3F;
  hotspot_min_ml_score = 0.01F;
  phased = false;
  use_vcf_features = false;
  iterations = 1000;
  snv_iterations = 0;
  indel_iterations = 0;
  snv_min_ml_score = 0.3F;
  indel_min_ml_score = 0.3F;
  min_tumor_support = 1;
  decode_yc = YcDecodeMethod::kNone;
  model_lgbm_params = "objective=binary boosting=gbdt learning_rate=0.01 seed=1238845 num_leaves=3 num_classes=1";
  model_lgbm_prediction_params =
      "num_threads=1 force_row_wise=true deterministic=true seed=1238845 pred_early_stop=true";
}

void SVCConfig::ConfigureTumorNormalWgsWorkflow() {
  // TODO : add scoring names, categorical names
  // TODO : add value for `--max-variants-per-read`
  bam_feature_names = kTumorNormalWgsBamFeatures;
  vcf_feature_names = kTumorNormalWgsVcfFeatures;
  workflow = Workflow::kTumorNormalWgs;
  scoring_names = kTumorNormalWgsScoringNames;
  categorical_names = kTumorNormalWgsCategoricalNames;
  n_classes = 1;
  min_mapq = 1;
  min_bq = kConcordantMinBaseQuality;
  min_dist = 2;
  min_family_size = kDuplexReadFamilySize;
  filter_homopolymer = HomopolymerFilter::kNone;
  min_homopolymer_length = 7;
  sequencing_protocol = SequencingProtocol::kDuplex;
  min_alt_counts = 3;
  min_phased_af = 0.001F;
  max_phased_af = 0.5F;
  min_weighted_counts = 4;
  hotspot_min_weighted_counts = 2;
  min_ml_score = 0.3F;
  hotspot_min_ml_score = 0.01F;
  phased = false;
  use_vcf_features = true;
  iterations = 3000;
  normalize_features = FeatureNormalization::kNone;
  // TODO : Update thresholds with appropriate default values
  snv_min_ml_score = 0.03F;
  indel_min_ml_score = 0.008F;
  min_tumor_support = 1;
  decode_yc = YcDecodeMethod::kNone;
  min_base_type = BaseType::kConcordant;
  model_lgbm_params = "objective=binary boosting=gbdt learning_rate=0.005 seed=1238845 num_classes=1";
  model_lgbm_prediction_params = "num_threads=1 force_row_wise=true deterministic=true seed=1238845";
}

void SVCConfig::ConfigureGermlineTaggingWorkflow() {
  bam_feature_names = kGermlineTaggingBamFeatures;
  vcf_feature_names = kGermlineTaggingVcfFeatures;
  workflow = Workflow::kGermlineTagging;
  scoring_names = kGermlineTaggingScoringNames;
  categorical_names = kGermlineTaggingCategoricalNames;
  n_classes = 3;
  min_mapq = 1;
  min_bq = kConcordantMinBaseQuality;
  min_dist = 2;
  min_family_size = 2;
  filter_homopolymer = HomopolymerFilter::kNone;
  min_homopolymer_length = 7;
  sequencing_protocol = SequencingProtocol::kDuplex;
  use_vcf_features = true;
  iterations = 1500;
  normalize_features = FeatureNormalization::kNone;
  snv_min_ml_score = 0.0F;
  indel_min_ml_score = 0.0F;
  min_tumor_support = 1;
  decode_yc = YcDecodeMethod::kNone;
  min_base_type = BaseType::kConcordant;
  model_lgbm_params =
      "objective=multiclass boosting=gbdt metric=multi_logloss seed=1238845 learning_rate=0.01 bagging_fraction=0.9 "
      "bagging_freq=1 feature_fraction=0.5 num_leaves=64 is_unbalance=true num_classes=3 ";
  model_lgbm_prediction_params =
      "num_threads=1 force_row_wise=true deterministic=true seed=1238845 pred_early_stop=true";
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
  for (auto& col_name : vcf_feature_names) {
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
    for (auto& col_name : scoring_names) {
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
    for (auto& col_name : snv_scoring_names) {
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
    for (auto& col_name : indel_scoring_names) {
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
  return config_collection;
}

/**
 * @brief Helper function to extract and set up the config for a given workflow from a JSON file.
 * @param config_json Path of JSON file. Path must exist.
 * @param workflow Workflow name
 * @return Config for the workflow
 */
static SVCConfig JsonToConfigHelper(const fs::path& config_json, const std::string& workflow) {
  if (!exists(config_json)) {
    // Must check here whether file exists instead of defining the `--config` option with `check(CLI::ExistingFile)`.
    // By design, `force_callback` would still be called even when the option was not used.
    // Therefore, `check(CLI::ExistingFile)` could lead to an error because of an empty input path.
    throw error::Error("JSON file not found: {}", config_json);
  }
  SVCConfigCollection config_collection = JsonToConfigCollection(config_json);
  auto& configs = config_collection.config_profiles;
  if (configs.empty()) {
    throw error::Error("No workflows found in config JSON file: {}", config_json);
  }
  if (!configs.contains(workflow)) {
    throw error::Error("Specified workflow '{}' not found in config JSON file: {}", workflow, config_json);
  }
  auto config = configs.at(workflow);
  config.ConfigureFeatureCols();
  return config;
}

SVCConfig JsonToConfig(const fs::path& config_json, const std::string& workflow) {
  if (config_json.empty()) {
    // Convert workflow from string to enum.
    auto workflow_enum = enum_util::ParseEnumName<Workflow>(workflow);
    if (workflow_enum.has_value()) {
      // Set up the default config for the workflow enum.
      return SVCConfig(workflow_enum.value());
    } else {
      throw error::Error("Specified workflow '{}' not supported");
    }
  }
  return JsonToConfigHelper(config_json, workflow);
}

SVCConfig JsonToConfig(const fs::path& config_json, const Workflow workflow) {
  if (config_json.empty()) {
    // Set up the default config for the workflow enum.
    return SVCConfig(workflow);
  }
  const std::string name = enum_util::FormatEnumName(workflow);
  return JsonToConfigHelper(config_json, name);
}

}  // namespace xoos::svc
