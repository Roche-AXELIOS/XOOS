#include "model-metadata.h"

#include <fstream>

#include "core/variant-feature-extraction.h"
#include "util/log-util.h"
#include "xoos/error/error.h"
#include "xoos/types/int.h"

using Json = nlohmann::json;

namespace xoos::svc {

// Suffix for the model metadata file, which is saved as `<model_path>.metadata.json`.
const std::string kMetadataFileNameSuffix = ".metadata.json";

/**
 * @brief Serialize model metadata, the specified code version and model architecture version, to a JSON file.
 * The metadata file is saved as `<model_path>.metadata.json`.
 * @param model_path Path to the model file (without `.dvc` extension).
 * @param code_version Version string of the codebase.
 * @param model_architecture_version Architecture version of the model, which is used to verify compatibility with the
 * current code version.
 */
void SerializeModelMetadata(const fs::path& model_path,
                            const std::string& code_version,
                            const u64 model_architecture_version) {
  const ModelMetadata metadata{.model_architecture_version = model_architecture_version, .code_version = code_version};
  Json json;
  to_json(json, metadata);
  std::ofstream fh(model_path.c_str() + kMetadataFileNameSuffix);
  fh << json.dump(2);
}

/**
 * @brief Serialize model metadata, the internal model architecture version and the specified code version, to a JSON
 * file.
 * @param model_path Path to the model file (without `.dvc` extension).
 * @param code_version Version string of the codebase.
 */
void SerializeModelMetadata(const fs::path& model_path, const std::string& code_version) {
  // Use the hardcoded model architecture version for compatibility checks
  SerializeModelMetadata(model_path, code_version, kModelArchitectureVersion);
}

/**
 * @brief Verify the model metadata against the current code version.
 * This function reads the model metadata file and checks if the build ID matches the current internal build ID.
 * It also checks if the code version from the model metadata matches the current code version.
 * If the build ID is older than the current internal build ID, an error is thrown.
 * If the code version does not match, a warning is logged.
 * @param model_path Path to the model file (without `.dvc` extension).
 * @param code_version Current code version.
 * @throws std::runtime_error if the model build ID is older than the current internal build ID.
 * @note The model metadata file is expected to be named `<model_path>.metadata.json`.
 *       It should contain the keys `build_id` and `code_version`.
 */
void VerifyModelMetadata(const fs::path& model_path, const std::string& code_version) {
  fs::path metadata_path = model_path.c_str() + kMetadataFileNameSuffix;
  if (!fs::exists(metadata_path)) {
    throw error::Error("Model metadata file '{}' does not exist. Please ensure the model is properly trained.",
                       metadata_path);
  }

  // Read the metadata file into JSON format
  std::ifstream fh(metadata_path);
  ModelMetadata metadata;
  from_json(Json::parse(fh), metadata);

  if (metadata.code_version != code_version) {
    // If the code version does not match, log a warning because different code versions may still be compatible.
    WarnAsErrorIfSet("The model code version '{}' and the current code version '{}' do not match.",
                     metadata.code_version,
                     code_version);
  }
  if (metadata.model_architecture_version != kModelArchitectureVersion) {
    throw error::Error(
        "Incompatible model architecture version '{}'. Current model architecture version is '{}'."
        " Please retrain your model.",
        metadata.model_architecture_version,
        kModelArchitectureVersion);
  }
}

/**
 * @brief Verify feature names in score calculator against the specified scoring columns.
 * This function retrieves the model feature names from the score calculator and compares them with the expected
 * scoring columns. If there is a mismatch, an error is thrown indicating the index and name of the mismatched features.
 * @param cal Score calculator
 * @param scoring_cols Vector of scoring column enums
 * @throws std::runtime_error if the feature names in the score calculator do not match the expected scoring columns.
 */
void VerifyModelFeatureNames(const ScoreCalculator& cal, const vec<UnifiedFeatureCols>& scoring_cols) {
  const auto num_cols = scoring_cols.size();
  auto model_feat_names = cal.GetModelFeatureNames(static_cast<int>(num_cols));
  for (size_t i = 0; i < num_cols; ++i) {
    if (GetCol(model_feat_names.at(i)) != scoring_cols.at(i)) {
      throw error::Error("The scoring feature name at index '{}' in the config '{}' and the model '{}' do not match",
                         i,
                         GetFeatureName(scoring_cols.at(i)),
                         model_feat_names.at(i));
    }
  }
}

/**
 * @brief Verify model compatibility by checking metadata and feature names.
 * This function checks if the model metadata matches the current code version and verifies that the feature names
 * in the model match the expected scoring columns.
 * @param model_path Path to the model file (without `.dvc` extension).
 * @param scoring_cols Vector of scoring column enums.
 * @throws std::runtime_error if model file not found, metadata is incompatible or feature names do not match.
 */
void VerifyModelCompatibility(const fs::path& model_path, const vec<UnifiedFeatureCols>& scoring_cols) {
  if (!fs::exists(model_path)) {
    throw error::Error("Model file '{}' does not exist.", model_path);
  }
  const ScoreCalculator cal(model_path, GetFeatureVecLength(scoring_cols));
  VerifyModelFeatureNames(cal, scoring_cols);
}

}  // namespace xoos::svc
