#pragma once

#include <map>
#include <string>

#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "core/variant-id.h"
#include "util/locked-tsv-writer.h"

namespace xoos::svc {

/** @brief Default positive label used when truth labels are not provided for positive samples. */
constexpr u32 kDefaultPositiveLabel = 1;

/** @brief Default negative label used for all negative samples. */
constexpr u32 kDefaultNegativeLabel = 0;

/**
 * @brief TrainingDataSet2 is a data structure for storing training data for machine learning models. It includes the
 * feature vectors, labels, and metadata needed for training. The class provides methods for inserting training data,
 * writing the data to a TSV file, and managing downsampling of over-represented labels in the training data.
 */
struct TrainingDataSet {
 public:
  // length of feature vector for one variant
  size_t feature_vec_size{};
  // training labels as floats to be compatible with LightGBM
  vec<f32> labels{};
  // training data matrix in numeric format, linearized (row-major) to be compatible with LightGBM
  vec<f64> data_matrix{};
  // labels that need downsampling; stored as integer to avoid floating point precision issues when comparing labels
  std::unordered_set<u32> downsampling_labels{};
  // map from VariantId to one or more indexes in the training labels
  std::unordered_map<VariantId, vec<size_t>> vid_to_label_idx{};

  /**
   * @brief Get the number of data points (i.e., variants) in the training dataset. This is equal to the size of the
   * training labels vector.
   * @return The number of data points in the training dataset
   */
  size_t GetDataPointCount() const;

  /**
   * @brief Get the counts of each label in the training data. The counts are returned as a map where the key is the
   * label and the value is the count of that label in the training data.
   * @return Map of label counts in the training data
   */
  std::map<u32, u64> GetLabelCounts() const;

  /**
   * @brief Determine the labels that need to be downsampled based on the distribution of labels in the training data.
   * Labels that are over-represented in the training data (i.e., have higher-than-expected counts) will be added to the
   * set of downsampling labels. If there are no labels in the training data, no labels will be added to the set of
   * downsampling labels.
   */
  void DetermineDownsamplingLabels();

  /**
   * @brief Decide whether to keep the existing sample or replace it with a new sample for a given variant and label
   * index. This is used to achieve downsampling of over-represented labels in the training data. The decision is based
   * on the position of the variant, which is used as a hash value to ensure that the same variants are consistently
   * kept or replaced across different runs, reducing stochasticity in multi-sample models.
   * @param vid The VariantId of the variant for which to decide whether to keep or replace the sample
   * @param feature_vec The feature vector of the new sample being inserted
   * @param label_idx The index of the label for the existing sample in the training labels vector
   * @return true if the existing sample was replaced with the new sample, false if the existing sample was kept and
   * the new sample was discarded
   */
  bool DecideAndSwapRow(const VariantId& vid, const vec<f64>& feature_vec, size_t label_idx);

  /**
   * @brief Insert a row of training data for a specific variant. The training data is a matrix where each row
   * represents a variant and each column represents a feature. The labels are the ground truth genotypes or variant
   * classifications for each variant.
   * @post If this is the first row inserted into the dataset, the feature vector size will be set based on the size of
   * the inserted feature vector.
   * @post If the label for the variant is in the downsampling set and there is already a sample for the variant with
   * the same label, one of the samples will be retained for training and the other sample will be skipped.
   * @post If the row is inserted, the feature vector and label will be added to the training data. If the row is
   * skipped due to downsampling, the training data will remain unchanged.
   * @param vid The VariantId of the variant to insert
   * @param feature_vec The feature vector for the variant to insert
   * @param label The label for the variant to insert
   * @return true if the row was inserted or replaced an existing sample, false if the row was skipped due to
   * downsampling
   * @throws error::Error if the feature vector is empty
   * @throws error::Error if the size of the feature vector does not match the expected feature vector size
   */
  bool InsertRow(const VariantId& vid, const vec<f64>& feature_vec, u32 label);

  /**
   * @brief Write the training data to a TSV file using the provided writer.
   * @post The training data is written to the specified TSV file.
   * @post The header row includes the variant ID, label, and feature names.
   * @post Each subsequent row represents a variant and includes the variant ID, label, and feature values.
   * @param writer The writer to use for writing the training data to a TSV file
   * @param feature_names The names of the features corresponding to each column in the data matrix
   * @param germline_genotype Flag indicating whether to convert numeric labels to genotype strings for germline
   * workflows.
   * @throws error::Error if the size of the feature names does not match the expected feature vector size
   */
  void WriteData(LockedTsvWriter& writer, const vec<std::string>& feature_names, bool germline_genotype) const;

  /**
   * @brief Clear all training data from the dataset. This is used to reset the dataset before re-using it for training
   * a different model type.
   */
  void ClearData();
};

}  // namespace xoos::svc
