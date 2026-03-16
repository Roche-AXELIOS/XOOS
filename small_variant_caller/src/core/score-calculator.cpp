#include "score-calculator.h"

#include <cassert>
#include <numeric>

#include <xoos/error/error.h>

#include "core/genotype.h"
#include "util/compress/compress.h"
#include "util/lightgbm-util.h"

namespace xoos::svc {

s32 CalculatePhredQualityScore(const double pred_prob) {
  // Phred quality scores in VCF are capped at 255
  static constexpr s32 kMaxPhredQuality = 255;
  if (pred_prob >= 1.0) {
    return kMaxPhredQuality;
  }
  if (pred_prob <= 0.0) {
    return 0;
  }
  // Phred quality score is defined as -10 * log10(1 - P(pred)), where P(pred) is the prediction probability
  // Note that 1 - P(pred) is the probability of error
  return std::min(static_cast<s32>(std::round(-10.0 * std::log10(1.0 - pred_prob))), kMaxPhredQuality);
}

s32 CalculateVariantQuality(const vec<double>& scores) {
  // The probability of the reference genotype is at index 0.
  // We need at least two scores to calculate the variant quality score, otherwise we cannot calculate the probability
  // of the variant (which is 1 - P(ref)).
  if (scores.size() < 2) {
    return 0;
  }
  // Sum the probabilities of all non-reference genotypes starting from index 1 to get P(variant).
  return CalculatePhredQualityScore(std::accumulate(scores.begin() + 1, scores.end(), 0.0));
}

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

GenotypeScore ScoreCalculator::CalculateScoreGermline(const vec<double>& features) const {
  vec<double> scores(_num_classes, 0);
  int64_t out_len = 0;
  lightgbm::BoosterPredictForMatSingleRowFast(_fast_config.get(), features.data(), &out_len, scores.data());
  assert(std::cmp_equal(out_len, _num_classes));
  // default best class is 0 (GT=0/0) if no class has a score above `min_score`
  size_t best_class = 0;
  double best_score = scores[0];
  for (size_t i = 1; std::cmp_less(i, _num_classes); ++i) {
    if (scores[i] > best_score) {
      best_score = scores[i];
      best_class = i;
    }
  }
  const auto gq = CalculatePhredQualityScore(scores[best_class]);
  const auto qual = CalculateVariantQuality(scores);
  return {IntToGenotype(best_class), scores[best_class], gq, qual};
}

}  // namespace xoos::svc
