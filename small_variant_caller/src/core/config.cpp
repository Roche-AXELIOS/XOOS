#include "config.h"

#include <fstream>

#include <xoos/error/error.h>
#include <xoos/util/string-functions.h>

#include "variant-info.h"

namespace xoos::svc {

/**
 * Checks whether any of the scoring columns are (derivatives of) VCF features
 * @return true if any scoring column is a VCF feature, false otherwise
 */
bool SVCConfig::HasVcfFeatureScoringCols() const {
  // determine whether any scoring columns are (derivatives of) VCF features
  return std::ranges::any_of(scoring_cols, IsVcfFeatureCol) || std::ranges::any_of(snv_scoring_cols, IsVcfFeatureCol) ||
         std::ranges::any_of(indel_scoring_cols, IsVcfFeatureCol);
}

/**
 * Populates a vector of UnifiedFeatureCol enums that match the requested column names for BAM features
 */
void SVCConfig::GetFeatureCols() {
  auto unsupported_fields = FindUnsupportedBamFeatureNames(feature_names);
  if (!unsupported_fields.empty()) {
    throw error::Error("Unsupported BAM feature names: {}", string::Join(unsupported_fields, ", "));
  }
  feature_cols.reserve(feature_names.size());
  for (auto& col_name : feature_names) {
    feature_cols.push_back(GetCol(col_name));
  }
}

/**
 * Populates a vector of UnifiedFeatureCol enums that match the requested column names for VCF features
 */
void SVCConfig::GetVcfFeatureCols() {
  auto unsupported_fields = FindUnsupportedVcfFeatureNames(vcf_feature_names);
  if (!unsupported_fields.empty()) {
    throw error::Error("Unsupported VCF feature names: {}", string::Join(unsupported_fields, ", "));
  }
  vcf_feature_cols.reserve(vcf_feature_names.size());
  for (auto& col_name : vcf_feature_names) {
    vcf_feature_cols.push_back(GetCol(col_name));
  }
}

/**
 * Populates a vector of UnifiedFeatureCol enums that match the requested column names for model scoring columns
 */
void SVCConfig::GetScoringCols() {
  // Extract the scoring column enums from the names provided in the config.
  // The enums are more compact and more efficient to work with than the strings.
  if (!scoring_names.empty()) {
    auto unsupported_fields = FindUnsupportedScoringFeatureNames(scoring_names);
    if (!unsupported_fields.empty()) {
      throw error::Error("Unsupported scoring feature names: {}", string::Join(unsupported_fields, ", "));
    }
    scoring_cols.reserve(scoring_names.size());
    for (auto& col_name : scoring_names) {
      scoring_cols.push_back(GetCol(col_name));
    }
  }
  // germline workflow's SNV model scoring columns
  if (!snv_scoring_names.empty()) {
    auto unsupported_fields = FindUnsupportedScoringFeatureNames(snv_scoring_names);
    if (!unsupported_fields.empty()) {
      throw error::Error("Unsupported SNV scoring feature names: {}", string::Join(unsupported_fields, ", "));
    }
    snv_scoring_cols.reserve(snv_scoring_names.size());
    for (auto& col_name : snv_scoring_names) {
      snv_scoring_cols.push_back(GetCol(col_name));
    }
  }
  // germline workflow's indel model scoring columns
  if (!indel_scoring_names.empty()) {
    auto unsupported_fields = FindUnsupportedScoringFeatureNames(indel_scoring_names);
    if (!unsupported_fields.empty()) {
      throw error::Error("Unsupported indel scoring feature names: {}", string::Join(unsupported_fields, ", "));
    }
    indel_scoring_cols.reserve(indel_scoring_names.size());
    for (auto& col_name : indel_scoring_names) {
      indel_scoring_cols.push_back(GetCol(col_name));
    }
  }
}

/**
 * Checks the paths to model files passed in to ensure that they follow the specified naming conventions
 * One file must be for SNVs and one for InDels.
 * The SNV model file must contain 'snv' in the name while the InDel must contain 'indel' in its name
 */
void SVCConfig::GetGermlineModelPaths(const std::vector<fs::path>& paths) {
  // distinguish SNV and indel models based on file names
  bool snv_found = false;
  bool indel_found = false;
  for (const auto& p : paths) {
    auto file_name = p.filename().string();
    if (file_name.find("snv") != std::string::npos) {
      snv_model_file = p;
      snv_found = true;
    } else if (file_name.find("indel") != std::string::npos) {
      indel_model_file = p;
      indel_found = true;
    }
  }
  if (!snv_found && !indel_found) {
    throw error::Error("Germline SNV and indel model file names must contain 'snv' and 'indel', respectively");
  }
  if (!snv_found) {
    throw error::Error("Germline SNV model file name must contain 'snv'");
  }
  if (!indel_found) {
    throw error::Error("Germline indel model file name must contain 'indel'");
  }
}

#ifdef SOMATIC_ENABLE
/**
 * Checks the paths to model files passed in to ensure that they follow the specified naming conventions
 * One file must be for Somatic and one for Germline Fail.
 * The Somatic model file must contain 'somatic' in the name while the Germline Fail must contain 'germline' in its name
 */
void SVCConfig::GetSomaticTNModelPaths(const std::vector<fs::path>& paths) {
  // distinguish Somatic and Germline Fail models based on file names
  bool somatic_found = false;
  bool germline_found = false;
  for (auto& p : paths) {
    auto file_name = p.filename().string();
    if (file_name.find("somatic") != std::string::npos) {
      model_file = p;
      somatic_found = true;
    } else if (file_name.find("indel") != std::string::npos) {
      germline_fail_model_file = p;
      germline_found = true;
    }
  }
  if (!somatic_found && !germline_found) {
    throw error::Error(
        "Germline Fail and Somatic model file names must contain 'germline' and 'somatic', respectively");
  }
  if (!somatic_found) {
    throw error::Error("Somatic model file name must contain 'somatic'");
  }
  if (!germline_found) {
    throw error::Error("Germline Fail model file name must contain 'germline'");
  }
}
#endif  // SOMATIC_ENABLE

/**
 * Wrapper helper to populate the BAM and VCF feature enum sets and model training and scoring feature sets with enum
 * values once a SVCConfig has been read in from file.
 */
void SVCConfig::SetUpWorkflow() {
  GetFeatureCols();
  GetVcfFeatureCols();
  GetScoringCols();
}

/**
 * @brief Extract a collection of configs from a JSON file.
 * @param config_json Path of JSON file
 * @return Collection of configs
 */
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
  config.SetUpWorkflow();
  return config;
}

/**
 * @brief Extract and set up the config for a given workflow from a JSON file. If an empty JSON path is provided,
 * then the default config values are used for the specified workflow.
 * @param config_json Path of JSON file
 * @param workflow Workflow name
 * @return Config for the workflow
 */
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

/**
 * @brief Extract and set up the config for a given workflow from a JSON file. If an empty JSON path is provided,
 * then the default config values are used for the specified workflow.
 * @param config_json Path of JSON file
 * @param workflow Workflow enum
 * @return Config for the workflow
 */
SVCConfig JsonToConfig(const fs::path& config_json, const Workflow workflow) {
  if (config_json.empty()) {
    // Set up the default config for the workflow enum.
    return SVCConfig(workflow);
  }
  const std::string name = enum_util::FormatEnumName(workflow);
  return JsonToConfigHelper(config_json, name);
}

}  // namespace xoos::svc
