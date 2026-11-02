#include "score-calculator.h"

#include <cmath>

#include <xoos/compress/compress.h>
#include <xoos/error/error.h>

#include "core/genotype.h"
#include "util/lightgbm-util.h"

namespace xoos::svc {

ScoreCalculator::ScoreCalculator(const fs::path& model_file,
                                 const size_t ncol,
                                 const std::string& prediction_params,
                                 const bool predict_contributions) {
  s32 iterations = 0;
  BoosterHandle booster = nullptr;
  if (compress::IsCompressed(model_file)) {
    lightgbm::BoosterLoadModelFromString(compress::Decompress(model_file), &iterations, &booster);
  } else {
    lightgbm::BoosterCreateFromModelFile(model_file, &iterations, &booster);
  }
  _booster.reset(booster);
  _num_classes = lightgbm::BoosterGetNumClasses(_booster.get());

  FastConfigHandle normal_config{};
  lightgbm::BoosterPredictForMatSingleRowFastInit(_booster.get(),
                                                  C_API_PREDICT_NORMAL,
                                                  0,
                                                  0,
                                                  C_API_DTYPE_FLOAT64,
                                                  static_cast<s32>(ncol),
                                                  prediction_params,
                                                  &normal_config);
  _predict_normal_config.reset(normal_config);

  _predict_contributions = predict_contributions;
  if (predict_contributions) {
    // Set up a separate configuration for predicting feature contributions (SHAP values)
    FastConfigHandle contrib_config{};
    lightgbm::BoosterPredictForMatSingleRowFastInit(_booster.get(),
                                                    C_API_PREDICT_CONTRIB,
                                                    0,
                                                    0,
                                                    C_API_DTYPE_FLOAT64,
                                                    static_cast<s32>(ncol),
                                                    prediction_params,
                                                    &contrib_config);
    _predict_contrib_config.reset(contrib_config);
  }
}

vec<std::string> ScoreCalculator::GetModelFeatureNames() const {
  return lightgbm::BoosterGetFeatureNames(_booster.get());
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
  constexpr f64 kPhredScale{-10.0};
  return std::min(kPhredScale * std::log10(1.0 - prob), kMaxPhredQuality);
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

PredictionScore ScoreCalculator::CalculateScore(const vec<f64>& features) const {
  f64 score = 0;
  vec<f64> shap_values{};

  s64 out_len = 0;
  lightgbm::BoosterPredictForMatSingleRowFast(_predict_normal_config.get(), features.data(), &out_len, &score);
  if (out_len != 1) {
    throw error::Error("LightGBM prediction output length mismatch: expected 1, got {}", out_len);
  }

  if (_predict_contributions) {
    // The output results consist of SHAP values for all features followed by the base value.
    const auto result_len = features.size() + 1;
    shap_values = vec<f64>(result_len, 0);
    out_len = 0;
    lightgbm::BoosterPredictForMatSingleRowFast(
        _predict_contrib_config.get(), features.data(), &out_len, shap_values.data());
    if (!std::cmp_equal(out_len, result_len)) {
      throw error::Error("LightGBM SHAP value output length mismatch: expected {}, got {}", result_len, out_len);
    }
  }

  const auto qual = CalculatePhredQualityScore(score);
  // round qual to 2 decimal places to avoid returning very long decimal values
  const auto rounded_qual = static_cast<f32>(std::round(qual * 100.0) / 100.0);

  return {.genotype = Genotype::kGTNA,
          .probability = score,
          .genotype_quality = 0,
          .variant_quality = rounded_qual,
          .shap_values = std::move(shap_values)};
}

PredictionScore ScoreCalculator::CalculateGermlineScore(const vec<f64>& features, const f64 min_score) const {
  vec<f64> scores(static_cast<size_t>(_num_classes), 0);
  s64 out_len = 0;
  lightgbm::BoosterPredictForMatSingleRowFast(_predict_normal_config.get(), features.data(), &out_len, scores.data());
  if (!std::cmp_equal(out_len, _num_classes)) {
    throw error::Error(
        "LightGBM germline prediction output length mismatch: expected {}, got {}", _num_classes, out_len);
  }
  // default best class is 0 (GT=0/0) if no class has a score above `min_score`
  size_t best_class = 0;
  f64 best_score = std::max(min_score, scores[0]);
  const auto num_classes = static_cast<size_t>(_num_classes);
  for (size_t i = 1; i < num_classes; ++i) {
    if (scores[i] > best_score) {
      best_score = scores[i];
      best_class = i;
    }
  }
  best_score = scores[best_class];

  vec<f64> shap_values{};
  if (_predict_contributions) {
    // The flattened output results consist of multiple blocks of values, one block per class, from class 1 to N.
    const auto num_features = features.size();
    const auto result_len = (num_features + 1) * static_cast<size_t>(_num_classes);
    vec<f64> out_result(result_len, 0);
    out_len = 0;
    lightgbm::BoosterPredictForMatSingleRowFast(
        _predict_contrib_config.get(), features.data(), &out_len, out_result.data());
    if (!std::cmp_equal(out_len, result_len)) {
      throw error::Error(
          "LightGBM germline SHAP contribution output length mismatch: expected {}, got {}", result_len, out_len);
    }

    // Each block contains the SHAP values for all features followed by the base value.
    const auto block_size = static_cast<s64>(num_features + 1);
    // Extract SHAP values and base value for the best class
    const auto block_start = static_cast<s64>(best_class) * block_size;
    shap_values = vec<f64>(out_result.begin() + block_start, out_result.begin() + block_start + block_size);
  }

  const auto gq = CalculateGenotypeQuality(best_score);
  const auto qual = CalculateVariantQuality(scores);

  return {.genotype = IntToGenotype(best_class),
          .probability = best_score,
          .genotype_quality = gq,
          .variant_quality = qual,
          .shap_values = std::move(shap_values)};
}

PredictionScore ScoreCalculator::CalculateGermlineScore(const vec<f64>& features) const {
  return CalculateGermlineScore(features, 0.0);
}

}  // namespace xoos::svc
