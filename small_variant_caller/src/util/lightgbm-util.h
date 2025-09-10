#pragma once

#include <LightGBM/c_api.h>

#include <xoos/types/fs.h>
#include <xoos/types/vec.h>

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
void BoosterLoadModelFromString(const std::string& model_str, int* out_num_iterators, BoosterHandle* out);

/// @copydoc LGBM_BoosterCreateFromModelFile
void BoosterCreateFromModelFile(const fs::path& filename, int* out_num_iterations, BoosterHandle* out);

/// @copydoc LGBM_BoosterGetNumClasses
int BoosterGetNumClasses(BoosterHandle handle);

/// @copydoc LGBM_BoosterPredictForMatSingleRowFastInit
void BoosterPredictForMatSingleRowFastInit(BoosterHandle handle,
                                           int predict_type,
                                           int start_iteration,
                                           int num_iteration,
                                           int data_type,
                                           int32_t ncol,
                                           const std::string& parameter,
                                           FastConfigHandle* out_fast_config);

/// @copydoc LGBM_BoosterPredictForMatSingleRowFast
void BoosterPredictForMatSingleRowFast(BoosterHandle handle, const double* data, int64_t* out_len, double* out_result);

/// @copydoc LGBM_BoosterGetNumFeature
int BoosterGetNumFeature(BoosterHandle handle);

/// @copydoc LGBM_BoosterGetFeatureNames
void BoosterGetFeatureNames(
    BoosterHandle handle, int len, int* out_len, size_t buffer_len, size_t* out_buffer_len, char** out_strs);

/// @copydoc LGBM_BoosterGetFeatureNames
vec<std::string> BoosterGetFeatureNames(BoosterHandle handle, size_t max_feature_name_len);

/// @copydoc LGBM_BoosterGetFeatureNames
vec<std::string> BoosterGetFeatureNames(BoosterHandle handle);

/// @copydoc LGBM_DatasetCreateFromMat
void DatasetCreateFromMat(const double* data,
                          int data_type,
                          int nrow,
                          int ncol,
                          int is_row_major,
                          const std::string& parameters,
                          DatasetHandle reference,
                          DatasetHandle* out);

/// @copydoc LGBM_DatasetSetFeatureNames
void DatasetSetFeatureNames(DatasetHandle handle, const char** feature_names, int num_feature_names);

/// @copydoc LGBM_DatasetSetField
void DatasetSetField(
    DatasetHandle handle, const std::string& field_name, const void* field_data, int num_element, int type);

/// @copydoc LGBM_BoosterCreate
void BoosterCreate(DatasetHandle train_data, const std::string& parameters, BoosterHandle* out);

/// @copydoc LGBM_BoosterUpdateOneIter
void BoosterUpdateOneIter(BoosterHandle handle, int* is_finished);

/// @copydoc LGBM_BoosterUpdateOneIter
bool BoosterUpdateOneIter(BoosterHandle handle);

int64_t BoosterSaveModelToString(
    BoosterHandle booster, int start_iteration, int num_iteration, int importance_type, vec<char>& out_result);

/// @copydoc LGBM_BoosterSaveModelToString
std::string BoosterSaveModelToString(BoosterHandle booster,
                                     int start_iteration,
                                     int num_iteration,
                                     int importance_type);

/// @copydoc LGBM_BoosterSaveModel
void BoosterSaveModel(
    BoosterHandle booster, int start_iteration, int num_iteration, int importance_type, const fs::path& output_file);

}  // namespace xoos::svc::lightgbm
