#pragma once

#include <nlohmann/json.hpp>

#include "score-calculator.h"
#include "variant-info.h"
#include "xoos/types/vec.h"

namespace xoos::svc {

// Model metadata to be stored in JSON format
struct ModelMetadata {
  u64 model_architecture_version;
  std::string code_version;
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ModelMetadata, model_architecture_version, code_version)

// The model architecture version is used to ensure that the model file is compatible with the current code version.
// It is in YYYYMMDDHH format, which was chosen to make it easy to identify when the build ID was last updated.
// It is hardcoded, and it should only be updated to deprecate older pre-train models (due to major changes to feature
// computation, model structure, LightGBM parameters, etc.).
constexpr u64 kModelArchitectureVersion = 2025073100;

void SerializeModelMetadata(const fs::path& model_path,
                            const std::string& code_version,
                            u64 model_architecture_version);
void SerializeModelMetadata(const fs::path& model_path, const std::string& code_version);
void VerifyModelMetadata(const fs::path& model_path, const std::string& code_version);
void VerifyModelFeatureNames(const ScoreCalculator& cal, const vec<UnifiedFeatureCols>& scoring_cols);
void VerifyModelCompatibility(const fs::path& model_path, const vec<UnifiedFeatureCols>& scoring_cols);

}  // namespace xoos::svc
