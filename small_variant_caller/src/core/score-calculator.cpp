#include "score-calculator.h"

#include <cassert>

#include <xoos/compress/compress.h>
#include <xoos/error/error.h>

#include "core/genotype.h"
#include "util/lightgbm-util.h"

namespace xoos::svc {

ScoreCalculator::ScoreCalculator(const fs::path& model_file, const size_t ncol, const std::string& prediction_params) {
  s32 iterations = 0;
  BoosterHandle booster;
  if (compress::IsCompressed(model_file)) {
    lightgbm::BoosterLoadModelFromString(compress::Decompress(model_file), &iterations, &booster);
  } else {
    lightgbm::BoosterCreateFromModelFile(model_file, &iterations, &booster);
  }
  _booster.reset(booster);
  _num_classes = lightgbm::BoosterGetNumClasses(_booster.get());
  FastConfigHandle fast_config;
  lightgbm::BoosterPredictForMatSingleRowFastInit(_booster.get(),
                                                  C_API_PREDICT_NORMAL,
                                                  0,
                                                  0,
                                                  C_API_DTYPE_FLOAT64,
                                                  static_cast<s32>(ncol),
                                                  prediction_params,
                                                  &fast_config);
  _fast_config.reset(fast_config);
}

/**
 * @brief Extract feature names from the LGBM model.
 * @return Vector of feature names extracted
 */
vec<std::string> ScoreCalculator::GetModelFeatureNames() const {
  return lightgbm::BoosterGetFeatureNames(_booster.get());
}

f64 ScoreCalculator::CalculateScore(const vec<f64>& features) const {
  f64 score;
  s64 out_len = 0;
  lightgbm::BoosterPredictForMatSingleRowFast(_fast_config.get(), features.data(), &out_len, &score);
  assert(out_len == 1);
  return score;
}

GenotypeScore ScoreCalculator::CalculateScoreGermline(const vec<f64>& features) const {
  return CalculateScoreGermline(features, 0.0);
}

f64 CalculatePhredQualityScore(const f64 prob) {
  // Phred quality scores in VCF are capped at 255
  static constexpr f64 kMaxPhredQuality = 255;
  if (prob >= 1.0) {
    return kMaxPhredQuality;
  }
  if (prob <= 0.0) {
    return 0;
  }
  // Phred quality score is defined as -10 * log10(1 - P), where P is the probability value
  // Note that 1 - P is the probability of error
  return std::min(-10.0 * std::log10(1.0 - prob), kMaxPhredQuality);
}

s32 CalculateGenotypeQuality(const f64 pred_prob) {
  return static_cast<s32>(std::round(CalculatePhredQualityScore(pred_prob)));
}

f32 CalculateVariantQuality(const vec<f64>& scores) {
  // The probability of the reference genotype is at index 0.
  // The probability of any variant genotype is at index > 0.
  // Sum the probabilities of all variant genotypes to get the overall probability of a variant existing at this
  // position.
  f64 var_score_sum = 0;
  for (size_t i = 1; i < scores.size(); ++i) {
    const f64 score = scores[i];
    if (score > 0.0) {
      var_score_sum += score;
    }
  }
  const auto qual = CalculatePhredQualityScore(var_score_sum);
  // round qual to 2 decimal places to avoid returning very long decimal values
  const auto rounded_qual = std::round(qual * 100.0) / 100.0;
  return static_cast<f32>(rounded_qual);
}

GenotypeScore ScoreCalculator::CalculateScoreGermline(const vec<f64>& features, const f64 min_score) const {
  vec<f64> scores(_num_classes, 0);
  s64 out_len = 0;
  lightgbm::BoosterPredictForMatSingleRowFast(_fast_config.get(), features.data(), &out_len, scores.data());
  assert(std::cmp_equal(out_len, _num_classes));
  // default best class is 0 (GT=0/0) if no class has a score above `min_score`
  size_t best_class = 0;
  f64 best_score = std::max(min_score, scores[0]);
  for (size_t i = 1; std::cmp_less(i, _num_classes); ++i) {
    if (scores[i] > best_score) {
      best_score = scores[i];
      best_class = i;
    }
  }
  const auto gq = CalculateGenotypeQuality(scores[best_class]);
  const auto qual = CalculateVariantQuality(scores);
  // not returning `best_score` as it is always >= min_score
  return {IntToGenotype(best_class), scores[best_class], gq, qual};
}

}  // namespace xoos::svc
