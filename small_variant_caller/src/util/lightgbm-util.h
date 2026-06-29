#pragma once

#include <LightGBM/c_api.h>

#include <xoos/types/fs.h>
#include <xoos/types/vec.h>

#include "xoos/types/float.h"
#include "xoos/types/int.h"

namespace xoos::svc::lightgbm {

struct BoosterHandleDeleter {
  void operator()(BoosterHandle handle) const;
};

using BoosterPtr = std::unique_ptr<std::remove_pointer_t<BoosterHandle>, BoosterHandleDeleter>;

struct FastConfigHandleDeleter {
  void operator()(FastConfigHandle fast_config) const;
};

using FastConfigPtr = std::unique_ptr<std::remove_pointer_t<FastConfigHandle>, FastConfigHandleDeleter>;

struct DatasetHandleDeleter {
  void operator()(DatasetHandle handle) const;
};

using DatasetPtr = std::unique_ptr<std::remove_pointer_t<DatasetHandle>, DatasetHandleDeleter>;

/// @copydoc LGBM_BoosterLoadModelFromString
void BoosterLoadModelFromString(const std::string& model_str, s32* out_num_iterations, BoosterHandle* out);

/// @copydoc LGBM_BoosterCreateFromModelFile
void BoosterCreateFromModelFile(const fs::path& filename, s32* out_num_iterations, BoosterHandle* out);

/// @copydoc LGBM_BoosterGetNumClasses
s32 BoosterGetNumClasses(BoosterHandle handle);

/// @copydoc LGBM_BoosterPredictForMatSingleRowFastInit
void BoosterPredictForMatSingleRowFastInit(BoosterHandle handle,
                                           s32 predict_type,
                                           s32 start_iteration,
                                           s32 num_iteration,
                                           s32 data_type,
                                           s32 ncol,
                                           const std::string& parameter,
                                           FastConfigHandle* out_fast_config);

/// @copydoc LGBM_BoosterPredictForMatSingleRowFast
void BoosterPredictForMatSingleRowFast(BoosterHandle handle, const f64* data, s64* out_len, f64* out_result);

/// @copydoc LGBM_BoosterGetNumFeature
s32 BoosterGetNumFeature(BoosterHandle handle);

/// @copydoc LGBM_BoosterGetFeatureNames
void BoosterGetFeatureNames(
    BoosterHandle handle, s32 len, s32* out_len, size_t buffer_len, size_t* out_buffer_len, char** out_strs);

/// @copydoc LGBM_BoosterGetFeatureNames
vec<std::string> BoosterGetFeatureNames(BoosterHandle handle, size_t max_feature_name_len);

/// @copydoc LGBM_BoosterGetFeatureNames
vec<std::string> BoosterGetFeatureNames(BoosterHandle handle);

/// @copydoc LGBM_DatasetCreateFromMat
void DatasetCreateFromMat(const f64* data,
                          s32 data_type,
                          s32 nrow,
                          s32 ncol,
                          s32 is_row_major,
                          const std::string& parameters,
                          DatasetHandle reference,
                          DatasetHandle* out);

/// @copydoc LGBM_DatasetSetFeatureNames
void DatasetSetFeatureNames(DatasetHandle handle, const char** feature_names, s32 num_feature_names);

/// @copydoc LGBM_DatasetSetField
void DatasetSetField(
    DatasetHandle handle, const std::string& field_name, const void* field_data, s32 num_element, s32 type);

/// @copydoc LGBM_BoosterCreate
void BoosterCreate(DatasetHandle train_data, const std::string& parameters, BoosterHandle* out);

/// @copydoc LGBM_BoosterUpdateOneIter
void BoosterUpdateOneIter(BoosterHandle handle, s32* is_finished);

/// @copydoc LGBM_BoosterUpdateOneIter
bool BoosterUpdateOneIter(BoosterHandle handle);

s64 BoosterSaveModelToString(
    BoosterHandle booster, s32 start_iteration, s32 num_iteration, s32 importance_type, vec<char>& out_result);

/// @copydoc LGBM_BoosterSaveModelToString
std::string BoosterSaveModelToString(BoosterHandle booster,
                                     s32 start_iteration,
                                     s32 num_iteration,
                                     s32 importance_type);

/// @copydoc LGBM_BoosterSaveModel
void BoosterSaveModel(
    BoosterHandle booster, s32 start_iteration, s32 num_iteration, s32 importance_type, const fs::path& output_file);

}  // namespace xoos::svc::lightgbm
