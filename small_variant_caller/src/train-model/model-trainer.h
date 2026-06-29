#pragma once

#include <optional>
#include <string>

#include <util/lightgbm-util.h>

#include <xoos/types/int.h>

#include "core/config.h"
#include "core/variant-info.h"
#include "core/workflow.h"
#include "training-data-set.h"

namespace xoos::svc {

class ModelTrainer {
 public:
  ModelTrainer(TrainingDataSet data, SVCConfig config, Workflow workflow);

  /**
   * @brief Train a LightGBM model based on the model training configurations. A subset of the pre-computed features are
   * numericalized and passed to LightGBM for training. The model is then trained on the data matrix using the specified
   * labels and LightGBM parameters, determined by the user selected workflow.
   * @param output_file The model file to save
   * @param num_threads Number of threads to use
   * @param num_rounds The maximum number of rounds of training.
   * @param var_group Type of variant to train
   * @param output_training_data Write training data to file
   */
  void Train(const fs::path& output_file,
             size_t num_threads,
             u32 num_rounds,
             VariantGroup var_group,
             const std::optional<fs::path>& output_training_data);

  /**
   * @brief Set up LightGBM scoring column names and categorical names based on the variant type, provided config file
   * and specified workflow
   * @param var_group Type of variant to train
   */
  void SetupScoringNamesAndCols(VariantGroup var_group);

  /**
   * @brief Clean up the trainer by clearing all internal data structures used during training. As we re-use the object
   * for training different model types internal data needs to be cleared before re-use.
   */
  void CleanupTrainer();

 private:
  const TrainingDataSet _dataset{};
  const SVCConfig _config{};
  Workflow _workflow{Workflow::kGermline};
  const bool _is_germline{};

  lightgbm::DatasetPtr _train_data{};
  lightgbm::BoosterPtr _booster{};
  vec<std::string> _scoring_names{};
  vec<FeatureColumn> _scoring_cols{};
  vec<std::string> _cat_names{};

  static constexpr u32 kNumTrainingRounds = 10000;
  static constexpr u32 kNumWriteOutput = 10000;
};

}  // namespace xoos::svc
