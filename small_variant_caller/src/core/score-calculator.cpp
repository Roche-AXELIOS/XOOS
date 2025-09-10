#include "score-calculator.h"

#include <cassert>

#include <xoos/error/error.h>

#include "core/genotype.h"
#include "util/compress/compress.h"
#include "util/lightgbm-util.h"

namespace xoos::svc {

ScoreCalculator::ScoreCalculator(const fs::path& model_file, const size_t ncol) {
  int iterations = 0;
  BoosterHandle booster;
  if (IsCompressed(model_file)) {
    lightgbm::BoosterLoadModelFromString(Decompress(model_file), &iterations, &booster);
  } else {
    lightgbm::BoosterCreateFromModelFile(model_file, &iterations, &booster);
  }
  _booster.reset(booster);
  _num_classes = lightgbm::BoosterGetNumClasses(_booster.get());
  FastConfigHandle fast_config;
  lightgbm::BoosterPredictForMatSingleRowFastInit(
      _booster.get(),
      C_API_PREDICT_NORMAL,
      0,
      0,
      C_API_DTYPE_FLOAT64,
      static_cast<int32_t>(ncol),
      "num_threads=1 force_row_wise=true deterministic=true seed=1238845 pred_early_stop=true",
      &fast_config);
  _fast_config.reset(fast_config);
}

/**
 * @brief Extract feature names from the LGBM model.
 * @param ncol Expected number of feature names
 * @return Vector of feature names extracted
 */
vec<std::string> ScoreCalculator::GetModelFeatureNames(const size_t ncol) const {
  return lightgbm::BoosterGetFeatureNames(_booster.get());
}

double ScoreCalculator::CalculateScore(const vec<double>& features) const {
  double score;
  int64_t out_len = 0;
  lightgbm::BoosterPredictForMatSingleRowFast(_fast_config.get(), features.data(), &out_len, &score);
  assert(out_len == 1);
  return score;
}

Genotype ScoreCalculator::CalculateScoreGermline(const vec<double>& features) const {
  vec<double> scores(_num_classes, 0);
  int64_t out_len = 0;
  lightgbm::BoosterPredictForMatSingleRowFast(_fast_config.get(), features.data(), &out_len, scores.data());
  assert(std::cmp_equal(out_len, _num_classes));
  int best_class = 0;
  double best_score = scores[0];
  for (int i = 1; i < _num_classes; ++i) {
    if (scores[i] > best_score) {
      best_score = scores[i];
      best_class = i;
    }
  }
  return IntToGenotype(static_cast<u64>(best_class));
}

}  // namespace xoos::svc
