#include "training-data-set.h"

#include <cmath>
#include <ranges>

#include "core/genotype.h"
#include "xoos/error/error.h"

namespace xoos::svc {

size_t TrainingDataSet::GetDataPointCount() const {
  return labels.size();
}

std::map<u32, u64> TrainingDataSet::GetLabelCounts() const {
  std::map<u32, u64> label_counts;
  for (const auto label : labels) {
    ++label_counts[static_cast<u32>(label)];
  }
  return label_counts;
}

void TrainingDataSet::DetermineDownsamplingLabels() {
  downsampling_labels.clear();
  if (labels.empty()) {
    return;
  }
  // count labels in the training data
  u64 total_count = 0;
  const auto& label_counts = GetLabelCounts();
  for (const u64 count : label_counts | std::views::values) {
    total_count += count;
  }
  // expected count for each label if the training data is perfectly balanced
  const u64 expected_count = total_count / label_counts.size();
  // find labels with higher-than-expected counts in the training data and add them to the set of downsampling labels
  for (const auto& [label, count] : label_counts) {
    if (count > expected_count) {
      downsampling_labels.emplace(label);
    }
  }
}

/**
 * @brief Helper function to compare a float label with an integer label for equality, considering floating point
 * precision issues.
 * @param label_f The float label to compare
 * @param label_u The integer label to compare
 * @return true if the labels are considered equal, false otherwise
 */
static bool Equal(const f32 label_f, const u32 label_u) {
  static constexpr f32 kEpsilon = 0.001f;
  // Check if the absolute difference is less than or equal to epsilon
  return std::abs(label_f - static_cast<f32>(label_u)) <= kEpsilon;
}

bool TrainingDataSet::DecideAndSwapRow(const VariantId& vid, const vec<f64>& feature_vec, const size_t label_idx) {
  // Use the module operator to decide whether to keep the existing sample or replace with the new sample.
  // Variant position is usually much larger than the number of samples, and it is sufficiently unique for all
  // variants in the training data. So, we can use it as a hash value in this context. Contrary to random
  // sampling, this helps reduce stochasticity in multi-sample models.
  if (vid.pos % 2 == 0) {
    // Overwrite the feature values of the existing sample in the data matrix with those of the new sample.
    const auto data_matrix_idx = label_idx * feature_vec_size;
    for (size_t i = 0; i < feature_vec_size; ++i) {
      data_matrix.at(data_matrix_idx + i) = feature_vec.at(i);
    }
    return true;
  }
  // Keep the existing sample and skip the new sample.
  return false;
}

bool TrainingDataSet::InsertRow(const VariantId& vid, const vec<f64>& feature_vec, const u32 label) {
  if (feature_vec.empty()) {
    throw error::Error(
        "Scoring feature vector is empty. Please check the scoring feature names in your configuration.");
  }

  if (feature_vec_size == 0) {
    // this is the first row to be inserted
    // set the feature vector size based on this row's feature vector size
    feature_vec_size = feature_vec.size();
  }

  if (feature_vec.size() != feature_vec_size) {
    throw error::Error(
        "Scoring feature vector size ({}) does not match the expected size ({}). "
        "Please check the scoring feature names in your configuration.",
        feature_vec.size(),
        feature_vec_size);
  }

  if (downsampling_labels.contains(label) && vid_to_label_idx.contains(vid)) {
    // this is a label that needs downsampling, and we already have a sample for this variant
    // check if the existing sample has the same label
    for (const auto label_idx : vid_to_label_idx.at(vid)) {
      if (Equal(labels.at(label_idx), label)) {
        // the variant already has this label
        // Keep either the existing sample or the new sample, to achieve downsampling of this label in the training data
        return DecideAndSwapRow(vid, feature_vec, label_idx);
      }
    }
  }
  // insert as a new data point for the variant
  vid_to_label_idx[vid].push_back(labels.size());
  data_matrix.insert(data_matrix.end(), feature_vec.begin(), feature_vec.end());
  labels.push_back(static_cast<f32>(label));
  return true;
}

void TrainingDataSet::WriteData(LockedTsvWriter& writer,
                                const vec<std::string>& feature_names,
                                const bool germline_genotype) const {
  if (feature_names.size() != feature_vec_size) {
    throw error::Error(
        "Scoring feature names size ({}) does not match the expected feature vector size ({}). "
        "Please check the scoring feature names in your configuration.",
        feature_names.size(),
        feature_vec_size);
  }

  // write header row to file
  // VARIANT, LABEL, FEATURE_1, FEATURE_2, ...
  vec<std::string> header{"VARIANT", "LABEL"};
  header.insert(header.end(), feature_names.begin(), feature_names.end());
  writer.AppendRow(header);

  // iterate through variants and write data rows to file
  for (const auto& [vid, label_idxs] : vid_to_label_idx) {
    for (const auto idx : label_idxs) {
      const auto label = labels[idx];
      const auto feature_vec_start = idx * feature_vec_size;
      const auto feature_vec_end = feature_vec_start + feature_vec_size;
      // variant ID is written as CHROM:POS:REF>ALT for better readability in the output file
      vec<std::string> data_row = {vid.ToString()};
      if (germline_genotype) {
        // convert numeric label to genotype string for germline workflow
        const auto genotype = IntToGenotype(static_cast<u32>(label));
        data_row.emplace_back(GenotypeToString(genotype));
      } else {
        data_row.emplace_back(std::to_string(label));
      }
      // append feature values to the data row
      for (size_t i = feature_vec_start; i < feature_vec_end; ++i) {
        data_row.emplace_back(std::to_string(data_matrix.at(i)));
      }
      writer.AppendRow(data_row);
    }
  }
}

void TrainingDataSet::ClearData() {
  data_matrix.clear();
  labels.clear();
  vid_to_label_idx.clear();
  downsampling_labels.clear();
  feature_vec_size = 0;
}

}  // namespace xoos::svc
