#include "lightgbm-util.h"

#include <cstddef>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>

#include "xoos/types/int.h"

namespace xoos::svc::lightgbm {

void BoosterHandleDeleter::operator()(
    BoosterHandle handle) const {  // NOSONAR: const on void* typedef yields void *const
  if (LGBM_BoosterFree(handle) == -1) {
    Logging::WarnNoExcept("Failed to free Booster handle: {}", LGBM_GetLastError());
  }
}

void FastConfigHandleDeleter::operator()(
    FastConfigHandle fast_config) const {  // NOSONAR: const on void* typedef yields void *const
  if (LGBM_FastConfigFree(fast_config) == -1) {
    Logging::WarnNoExcept("Failed to free FastConfig handle: {}", LGBM_GetLastError());
  }
}

void DatasetHandleDeleter::operator()(
    DatasetHandle handle) const {  // NOSONAR: const on void* typedef yields void *const
  if (LGBM_DatasetFree(handle) == -1) {
    Logging::WarnNoExcept("Failed to free Dataset handle: {}", LGBM_GetLastError());
  }
}

void BoosterLoadModelFromString(const std::string& model_str, s32* const out_num_iterations, BoosterHandle* const out) {
  if (LGBM_BoosterLoadModelFromString(model_str.c_str(), out_num_iterations, out) == -1) {
    throw error::Error("Failed to load the LightGBM model from string: {}", LGBM_GetLastError());
  }
}

void BoosterCreateFromModelFile(const fs::path& filename, s32* const out_num_iterations, BoosterHandle* const out) {
  if (LGBM_BoosterCreateFromModelfile(filename.c_str(), out_num_iterations, out) == -1) {
    throw error::Error("Failed to create Booster from model file: {}", LGBM_GetLastError());
  }
}

s32 BoosterGetNumClasses(BoosterHandle handle) {  // NOSONAR: const on void* typedef yields void *const, not const void*
  s32 num_classes = 0;
  if (LGBM_BoosterGetNumClasses(handle, &num_classes) == -1) {
    throw error::Error("Failed to get number of classes from Booster: {}", LGBM_GetLastError());
  }
  return num_classes;
}

void BoosterPredictForMatSingleRowFastInit(BoosterHandle handle,  // NOSONAR: const on void* typedef yields void *const
                                           const s32 predict_type,
                                           const s32 start_iteration,
                                           const s32 num_iteration,
                                           const s32 data_type,
                                           const s32 ncol,
                                           const std::string& parameter,
                                           FastConfigHandle* const out_fast_config) {
  const auto result = LGBM_BoosterPredictForMatSingleRowFastInit(
      handle, predict_type, start_iteration, num_iteration, data_type, ncol, parameter.c_str(), out_fast_config);
  if (result == -1) {
    throw error::Error("Failed to initialize FastConfig for Booster prediction: {}", LGBM_GetLastError());
  }
}

void BoosterPredictForMatSingleRowFast(BoosterHandle handle,  // NOSONAR: const on void* typedef yields void *const
                                       const f64* const data,
                                       s64* const out_len,
                                       f64* const out_result) {
  if (LGBM_BoosterPredictForMatSingleRowFast(handle, data, out_len, out_result) ==
      -1) {  // NOSONAR(cpp:S5356) — void* cast required by the LightGBM C API signature
    throw error::Error("Failed to predict for single row with FastConfig: {}", LGBM_GetLastError());
  }
}

void BoosterGetFeatureNames(BoosterHandle handle,  // NOSONAR: const on void* typedef yields void *const
                            const s32 len,
                            s32* const out_len,
                            const size_t buffer_len,
                            size_t* const out_buffer_len,
                            char** const out_strs) {
  const auto result = LGBM_BoosterGetFeatureNames(handle, len, out_len, buffer_len, out_buffer_len, out_strs);
  if (result == -1) {
    throw error::Error("Failed to get feature names from Booster: {}", LGBM_GetLastError());
  }
}

s32 BoosterGetNumFeature(BoosterHandle handle) {  // NOSONAR: const on void* typedef yields void *const
  s32 num_features = 0;
  if (LGBM_BoosterGetNumFeature(handle, &num_features) == -1) {
    throw error::Error("Failed to get number of features from Booster: {}", LGBM_GetLastError());
  }
  return num_features;
}

vec<std::string> BoosterGetFeatureNames(
    BoosterHandle handle, const size_t max_feature_name_len) {  // NOSONAR: const on void* typedef yields void *const
  const auto num_features = lightgbm::BoosterGetNumFeature(handle);

  vec<vec<char>> string_buffers(ToUnsigned(num_features), vec<char>(max_feature_name_len));
  vec<char*> c_strings(ToUnsigned(num_features));
  for (s32 i = 0; i < num_features; ++i) {
    c_strings[ToUnsigned(i)] = string_buffers[ToUnsigned(i)].data();
  }

  s32 out_len = 0;
  size_t out_buffer_len = 0;
  BoosterGetFeatureNames(handle, num_features, &out_len, max_feature_name_len, &out_buffer_len, c_strings.data());
  if (out_len != num_features) {
    throw error::Error("Feature name count mismatch: expected {}, got {}", num_features, out_len);
  }

  vec<std::string> result;
  result.reserve(ToUnsigned(out_len));
  for (s32 i = 0; i < out_len; ++i) {
    result.emplace_back(c_strings[ToUnsigned(i)]);
  }
  return result;
}

vec<std::string> BoosterGetFeatureNames(BoosterHandle handle) {  // NOSONAR: const on void* typedef yields void *const
  constexpr s32 kMaxFeatureNameLen = 256;
  return BoosterGetFeatureNames(handle, kMaxFeatureNameLen);
}

void DatasetCreateFromMat(const f64* const data,  // NOSONAR: signature mirrors LGBM_DatasetCreateFromMat C API
                          const s32 data_type,
                          const s32 nrow,
                          const s32 ncol,
                          const s32 is_row_major,
                          const std::string& parameters,
                          DatasetHandle reference,  // NOSONAR: const on void* typedef yields void *const
                          DatasetHandle* const out) {
  if (LGBM_DatasetCreateFromMat(data, data_type, nrow, ncol, is_row_major, parameters.c_str(), reference, out) ==
      -1) {  // NOSONAR(cpp:S5356) — void* cast required by the LightGBM C API signature
    throw error::Error("Failed to create dataset from matrix: {}", LGBM_GetLastError());
  }
}

void DatasetSetFeatureNames(DatasetHandle handle,  // NOSONAR: const on void* typedef yields void *const
                            const char** const feature_names,
                            const s32 num_feature_names) {
  if (LGBM_DatasetSetFeatureNames(handle, feature_names, num_feature_names) == -1) {
    throw error::Error("Failed to set feature names for Dataset: {}", LGBM_GetLastError());
  }
}

void DatasetSetField(DatasetHandle handle,  // NOSONAR: const on void* typedef yields void *const
                     const std::string& field_name,
                     const void* const field_data,  // NOSONAR: void* required by LGBM_DatasetSetField C API
                     const s32 num_element,
                     const s32 type) {
  if (LGBM_DatasetSetField(handle, field_name.c_str(), field_data, num_element, type) == -1) {
    throw error::Error("Failed to set field '{}' in dataset: {}", field_name, LGBM_GetLastError());
  }
}

void BoosterCreate(DatasetHandle train_data,
                   const std::string& parameters,
                   BoosterHandle* const out) {  // NOSONAR: const on void* typedef yields void *const
  if (LGBM_BoosterCreate(train_data, parameters.c_str(), out) == -1) {
    throw error::Error("Failed to create Booster from dataset: {}", LGBM_GetLastError());
  }
}

void BoosterUpdateOneIter(BoosterHandle handle,  // NOSONAR: const on void* typedef yields void *const
                          s32* const is_finished) {
  if (LGBM_BoosterUpdateOneIter(handle, is_finished) == -1) {
    throw error::Error("Failed to update Booster: {}", LGBM_GetLastError());
  }
}

bool BoosterUpdateOneIter(BoosterHandle handle) {  // NOSONAR: const on void* typedef yields void *const
  s32 is_finished = 0;
  BoosterUpdateOneIter(handle, &is_finished);
  return is_finished == 1;
}

s64 BoosterSaveModelToString(BoosterHandle booster,  // NOSONAR: const on void* typedef yields void *const
                             const s32 start_iteration,
                             const s32 num_iteration,
                             const s32 importance_type,
                             vec<char>& out_result) {
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

std::string BoosterSaveModelToString(BoosterHandle booster,  // NOSONAR: const on void* typedef yields void *const
                                     const s32 start_iteration,
                                     const s32 num_iteration,
                                     const s32 importance_type) {
  constexpr vec<char>::size_type kInitialModelBufferSize{vec<char>::size_type{1024} * 1024};
  vec<char> model_buffer(kInitialModelBufferSize);
  const auto out_len = BoosterSaveModelToString(booster, start_iteration, num_iteration, importance_type, model_buffer);
  if (std::cmp_less(model_buffer.size(), out_len)) {
    model_buffer.resize(ToUnsigned(out_len));
    BoosterSaveModelToString(booster, start_iteration, num_iteration, importance_type, model_buffer);
  }
  return {model_buffer.data(), static_cast<size_t>(out_len)};
}

void BoosterSaveModel(BoosterHandle booster,  // NOSONAR: const on void* typedef yields void *const
                      const s32 start_iteration,
                      const s32 num_iteration,
                      const s32 importance_type,
                      const fs::path& output_file) {
  const auto ret = LGBM_BoosterSaveModel(booster, start_iteration, num_iteration, importance_type, output_file.c_str());
  if (ret == -1) {
    throw error::Error("Failed to save LightGBM model to {}: {}", output_file, LGBM_GetLastError());
  }
}

}  // namespace xoos::svc::lightgbm
