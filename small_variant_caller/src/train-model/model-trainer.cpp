#include "model-trainer.h"

#include <zlib.h>
#include <zstd.h>

#include <LightGBM/c_api.h>
#include <fmt/format.h>

#include <csv.hpp>

#include <xoos/compress/compress.h>
#include <xoos/error/error.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/log/logging.h>
#include <xoos/util/string-functions.h>

#include "util/lightgbm-util.h"
#include "util/log-util.h"

namespace xoos::svc {

// Same variant of a genotype repeated in multiple samples may lead to over-representation in the training data.
// This may lead to the model converging too early, which yields a lower than expected F1 score. To avoid this, we
// reduce redundancy by keeping only one sample for each variant with the same genotype. Sample count is used to
// determine if we need to for germline models.
ModelTrainer::ModelTrainer(TrainingDataSet data, SVCConfig config, Workflow workflow)
    : _dataset(std::move(data)),
      _config(std::move(config)),
      _workflow(workflow),
      _is_germline(IsGermlineWorkflow(_workflow)) {
}

/**
 * @brief Get the workflow-specific LightGBM configuration string from the model training configuration.
 * @param workflow Workflow type
 * @param cat_names List of categorical feature names
 * @param config Model training configuration and parameters
 * @param var_type Variant type
 * @param num_threads Number of threads for model training
 * @return LightGBM configuration string
 */
static std::string GetLightGBMConfig(const Workflow workflow,
                                     const std::vector<std::string>& cat_names,
                                     const SVCConfig& config,
                                     const VariantGroup var_type,
                                     const size_t num_threads) {
  using enum Workflow;
  std::string params{fmt::format("num_threads={} force_row_wise=true deterministic=true ", num_threads)};
  switch (workflow) {
    case kGermlineMultiSample:
    case kGermline: {
      // Use Different settings for indel vs SNV models
      params += var_type == VariantGroup::kIndelOnly ? config.indel_model_lgbm_params : config.snv_model_lgbm_params;
      break;
    }
    case kTumorOnlyTe:
    case kTumorNormalWgs:
    case kGermlineTagging: {
      params += config.model_lgbm_params;
      break;
    }
    default:
      throw error::Error("Unsupported workflow");
  }
  if (!cat_names.empty()) {
    params += " categorical_feature=name:" + string::Join(cat_names, ",");
  }
  return params;
}

/**
 * @brief Helper function to serialize the LightGBM booster to a string.
 * @param booster LightGBM booster
 * @return Model content as a string
 */
static std::string SaveModelToString(BoosterHandle booster) {
  return lightgbm::BoosterSaveModelToString(booster, 0, -1, C_API_FEATURE_IMPORTANCE_SPLIT);
}

void ModelTrainer::SetupScoringNamesAndCols(const VariantGroup var_group) {
  switch (var_group) {
    case VariantGroup::kSnvOnly:
      _scoring_names = _config.snv_scoring_names;
      _scoring_cols = _config.snv_scoring_cols;
      _cat_names = _config.snv_categorical_names;
      break;
    case VariantGroup::kIndelOnly:
      _scoring_names = _config.indel_scoring_names;
      _scoring_cols = _config.indel_scoring_cols;
      _cat_names = _config.indel_categorical_names;
      break;
    default:
      _scoring_names = _config.scoring_names;
      _scoring_cols = _config.scoring_cols;
      _cat_names = _config.categorical_names;
      break;
  }
}

void ModelTrainer::CleanupTrainer() {
  _scoring_names.clear();
  _scoring_cols.clear();
  _cat_names.clear();
}

void ModelTrainer::Train(const fs::path& output_file,
                         const size_t num_threads,
                         const u32 num_rounds,
                         const VariantGroup var_group,
                         const std::optional<fs::path>& output_training_data) {
  // 1. Set up LightGBM scoring column names and categorical names based on the variant type.
  // 2. Write training data to a TSV file if the output_training_data path is provided in the parameters.
  // 3. Print the label distribution of the training data for debugging and verification purposes.
  // 4. Create LightGBM training dataset.
  // 5. Create LightGBM booster.
  // 6. Train the model for the specified number of rounds.
  // 7. Save the model to the output file.

  // set up scoring names, scoring columns, and categorical names
  SetupScoringNamesAndCols(var_group);

  // set up lightGBM column names
  vec<const char*> columns;
  columns.reserve(_scoring_names.size());
  for (const auto& col : _scoring_names) {
    columns.push_back(col.c_str());
  }

  if (output_training_data.has_value()) {
    auto writer = LockedTsvWriter(output_training_data.value());
    Logging::Info("Writing training data to {}", output_training_data.value());
    _dataset.WriteData(writer, _scoring_names, _is_germline);
  }

  {
    // Extract and print the label counts of training data
    auto label_counts = _dataset.GetLabelCounts();
    vec<std::string> items{};
    if (_is_germline) {
      for (const auto& [label, count] : label_counts) {
        const auto genotype = IntToGenotype(label);
        items.emplace_back(fmt::format("{}:{}", GenotypeToString(genotype), count));
      }
    } else {
      for (const auto& [label, count] : label_counts) {
        items.emplace_back(fmt::format("{}:{}", label, count));
      }
    }
    Logging::Info("Label counts: {}", string::Join(items, ", "));
    if (label_counts.size() != _config.n_classes) {
      WarnAsErrorIfSet(
          "Number of unique labels ({}) and classes ({}) are different!", label_counts.size(), _config.n_classes);
    }
  }
  auto lightgbm_config_str = GetLightGBMConfig(_workflow, _cat_names, _config, var_group, num_threads);
  Logging::Info("Training with config: {}", lightgbm_config_str);

  const auto num_rows = static_cast<s32>(_dataset.labels.size());
  const auto num_cols = static_cast<s32>(_dataset.feature_vec_size);
  DatasetHandle train_data = nullptr;
  lightgbm::DatasetCreateFromMat(
      _dataset.data_matrix.data(), C_API_DTYPE_FLOAT64, num_rows, num_cols, 1, "", nullptr, &train_data);
  _train_data.reset(train_data);

  lightgbm::DatasetSetFeatureNames(_train_data.get(), columns.data(), num_cols);
  lightgbm::DatasetSetField(_train_data.get(), "label", _dataset.labels.data(), num_rows, C_API_DTYPE_FLOAT32);

  BoosterHandle booster = nullptr;
  lightgbm::BoosterCreate(_train_data.get(), lightgbm_config_str, &booster);
  _booster.reset(booster);

  for (u32 train_round = 0; train_round < num_rounds; ++train_round) {
    if (0 == (train_round % 100)) {
      Logging::Info("Training round {}", train_round);
    }
    if (lightgbm::BoosterUpdateOneIter(_booster.get())) {
      Logging::Info("Training finished at round {}", train_round);
      break;
    }
  }
  Logging::Info("Saving model to {}", output_file);
  if (compress::IsCompressed(output_file)) {
    compress::Compress(SaveModelToString(_booster.get()), output_file);
  } else {
    lightgbm::BoosterSaveModel(_booster.get(), 0, -1, C_API_FEATURE_IMPORTANCE_SPLIT, output_file);
  }

  CleanupTrainer();
}

}  // namespace xoos::svc
