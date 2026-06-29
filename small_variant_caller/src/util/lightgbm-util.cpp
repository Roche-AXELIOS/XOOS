#include "lightgbm-util.h"

#include <utility>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>

#include "xoos/types/int.h"

namespace xoos::svc::lightgbm {

void BoosterHandleDeleter::operator()(BoosterHandle handle) const {
  if (LGBM_BoosterFree(handle) == -1) {
    Logging::WarnNoExcept("Failed to free Booster handle: {}", LGBM_GetLastError());
  }
}

void FastConfigHandleDeleter::operator()(FastConfigHandle fast_config) const {
  if (LGBM_FastConfigFree(fast_config) == -1) {
    Logging::WarnNoExcept("Failed to free FastConfig handle: {}", LGBM_GetLastError());
  }
}

void DatasetHandleDeleter::operator()(DatasetHandle handle) const {
  if (LGBM_DatasetFree(handle) == -1) {
    Logging::WarnNoExcept("Failed to free Dataset handle: {}", LGBM_GetLastError());
  }
}

void BoosterLoadModelFromString(const std::string& model_str, s32* out_num_iterations, BoosterHandle* out) {
  if (LGBM_BoosterLoadModelFromString(model_str.c_str(), out_num_iterations, out) == -1) {
    throw error::Error("Failed to load the LightGBM model from string: {}", LGBM_GetLastError());
  }
}

void BoosterCreateFromModelFile(const fs::path& filename, s32* out_num_iterations, BoosterHandle* out) {
  if (LGBM_BoosterCreateFromModelfile(filename.c_str(), out_num_iterations, out) == -1) {
    throw error::Error("Failed to create Booster from model file: {}", LGBM_GetLastError());
  }
}

s32 BoosterGetNumClasses(BoosterHandle handle) {
  s32 num_classes = 0;
  if (LGBM_BoosterGetNumClasses(handle, &num_classes) == -1) {
    throw error::Error("Failed to get number of classes from Booster: {}", LGBM_GetLastError());
  }
  return num_classes;
}

void BoosterPredictForMatSingleRowFastInit(BoosterHandle handle,
                                           const s32 predict_type,
                                           const s32 start_iteration,
                                           const s32 num_iteration,
                                           const s32 data_type,
                                           const s32 ncol,
                                           const std::string& parameter,
                                           FastConfigHandle* out_fast_config) {
  const auto result = LGBM_BoosterPredictForMatSingleRowFastInit(
      handle, predict_type, start_iteration, num_iteration, data_type, ncol, parameter.c_str(), out_fast_config);
  if (result == -1) {
    throw error::Error("Failed to initialize FastConfig for Booster prediction: {}", LGBM_GetLastError());
  }
}

void BoosterPredictForMatSingleRowFast(BoosterHandle handle, const f64* data, s64* out_len, f64* out_result) {
  if (LGBM_BoosterPredictForMatSingleRowFast(handle, data, out_len, out_result) == -1) {
    throw error::Error("Failed to predict for single row with FastConfig: {}", LGBM_GetLastError());
  }
}

void BoosterGetFeatureNames(BoosterHandle handle,
                            const s32 len,
                            s32* out_len,
                            const size_t buffer_len,
                            size_t* out_buffer_len,
                            char** out_strs) {
  const auto result = LGBM_BoosterGetFeatureNames(handle, len, out_len, buffer_len, out_buffer_len, out_strs);
  if (result == -1) {
    throw error::Error("Failed to get feature names from Booster: {}", LGBM_GetLastError());
  }
}

s32 BoosterGetNumFeature(BoosterHandle handle) {
  s32 num_features = 0;
  if (LGBM_BoosterGetNumFeature(handle, &num_features) == -1) {
    throw error::Error("Failed to get number of features from Booster: {}", LGBM_GetLastError());
  }
  return num_features;
}

vec<std::string> BoosterGetFeatureNames(BoosterHandle handle, size_t max_feature_name_len) {
  const auto num_features = lightgbm::BoosterGetNumFeature(handle);

  vec<vec<char>> string_buffers(num_features, vec<char>(max_feature_name_len));
  vec<char*> c_strings(num_features);
  for (s32 i = 0; i < num_features; ++i) {
    c_strings[i] = string_buffers[i].data();
  }

  s32 out_len;
  size_t out_buffer_len;
  BoosterGetFeatureNames(handle, num_features, &out_len, max_feature_name_len, &out_buffer_len, c_strings.data());
  assert(out_len == num_features);

  vec<std::string> result;
  result.reserve(out_len);
  for (s32 i = 0; i < out_len; ++i) {
    result.emplace_back(c_strings[i]);
  }
  return result;
}

vec<std::string> BoosterGetFeatureNames(BoosterHandle handle) {
  constexpr s32 kMaxFeatureNameLen = 256;
  return BoosterGetFeatureNames(handle, kMaxFeatureNameLen);
}

void DatasetCreateFromMat(const f64* data,
                          s32 data_type,
                          s32 nrow,
                          s32 ncol,
                          s32 is_row_major,
                          const std::string& parameters,
                          DatasetHandle reference,
                          DatasetHandle* out) {
  if (LGBM_DatasetCreateFromMat(data, data_type, nrow, ncol, is_row_major, parameters.c_str(), reference, out) == -1) {
    throw error::Error("Failed to create dataset from matrix: {}", LGBM_GetLastError());
  }
}

void DatasetSetFeatureNames(DatasetHandle handle, const char** feature_names, s32 num_feature_names) {
  if (LGBM_DatasetSetFeatureNames(handle, feature_names, num_feature_names) == -1) {
    throw error::Error("Failed to set feature names for Dataset: {}", LGBM_GetLastError());
  }
}

void DatasetSetField(
    DatasetHandle handle, const std::string& field_name, const void* field_data, s32 num_element, s32 type) {
  if (LGBM_DatasetSetField(handle, field_name.c_str(), field_data, num_element, type) == -1) {
    throw error::Error("Failed to set field '{}' in dataset: {}", field_name, LGBM_GetLastError());
  }
}

void BoosterCreate(DatasetHandle train_data, const std::string& parameters, BoosterHandle* out) {
  if (LGBM_BoosterCreate(train_data, parameters.c_str(), out) == -1) {
    throw error::Error("Failed to create Booster from dataset: {}", LGBM_GetLastError());
  }
}

void BoosterUpdateOneIter(BoosterHandle handle, s32* is_finished) {
  if (LGBM_BoosterUpdateOneIter(handle, is_finished) == -1) {
    throw error::Error("Failed to update Booster: {}", LGBM_GetLastError());
  }
}

bool BoosterUpdateOneIter(BoosterHandle handle) {
  s32 is_finished = 0;
  BoosterUpdateOneIter(handle, &is_finished);
  return is_finished == 1;
}

s64 BoosterSaveModelToString(
    BoosterHandle booster, s32 start_iteration, s32 num_iteration, s32 importance_type, vec<char>& out_result) {
  s64 out_len = 0;
  const auto result = LGBM_BoosterSaveModelToString(booster,
                                                    start_iteration,
                                                    num_iteration,
                                                    importance_type,
                                                    static_cast<s64>(out_result.size()),
                                                    &out_len,
                                                    out_result.data());
  if (result == -1) {
    throw error::Error("Failed to save model to string: {}", LGBM_GetLastError());
  }
  return out_len;
}

std::string BoosterSaveModelToString(BoosterHandle booster,
                                     s32 start_iteration,
                                     s32 num_iteration,
                                     s32 importance_type) {
  vec<char> model_buffer(1024 * 1024);
  const auto out_len = BoosterSaveModelToString(booster, start_iteration, num_iteration, importance_type, model_buffer);
  if (std::cmp_less(model_buffer.size(), out_len)) {
    model_buffer.resize(out_len);
    BoosterSaveModelToString(booster, start_iteration, num_iteration, importance_type, model_buffer);
  }
  return {model_buffer.data(), static_cast<size_t>(out_len)};
}

void BoosterSaveModel(
    BoosterHandle booster, s32 start_iteration, s32 num_iteration, s32 importance_type, const fs::path& output_file) {
  const auto ret = LGBM_BoosterSaveModel(booster, start_iteration, num_iteration, importance_type, output_file.c_str());
  if (ret == -1) {
    throw error::Error("Failed to save LightGBM model to {}: {}", output_file, LGBM_GetLastError());
  }
}

}  // namespace xoos::svc::lightgbm
