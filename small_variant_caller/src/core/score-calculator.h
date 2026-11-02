#pragma once

#include <string>

#include <xoos/types/fs.h>
#include <xoos/types/vec.h>

#include "core/genotype.h"
#include "util/lightgbm-util.h"
#include "xoos/types/float.h"

namespace xoos::svc {

/**
 * @brief Struct to hold the predicted genotype and associated scores for a germline variant.
 * @see Genotype
 * @see ScoreCalculator
 */
struct PredictionScore {
  // predicted genotype based on prediction probability score
  Genotype genotype = Genotype::kGTNA;
  // prediction probability score for the predicted genotype
  f64 probability = 0;
  // Phred score for genotype quality, calculated as -10 * log10(1 - P(genotype)), where P(genotype) is the prediction
  // probability score for the predicted genotype
  s32 genotype_quality = 0;
  // Phred score for variant quality, calculated as -10 * log10(1 - P(variant)), where P(variant) is the sum of
  // probabilities of all non-reference genotypes
  f32 variant_quality = 0;
  // vector of SHAP values for each feature followed by the base margin (the last value in the vector) for the predicted
  // genotype
  vec<f64> shap_values;
  // AF-based genotype from the afq_gt BAM feature, used to rescue multiallele_partner variants
  Genotype afq_gt = Genotype::kGT00;
};

class ScoreCalculator {
 public:
  ScoreCalculator(const fs::path& model_file,
                  size_t ncol,
                  const std::string& prediction_params,
                  bool predict_contributions);

  /**
   * @brief Extract feature names from the LightGBM model.
   * @return Vector of feature names extracted
   */
  vec<std::string> GetModelFeatureNames() const;

  /**
   * @brief Calculate the prediction probability score for a single set of features.
   * @pre The LightGBM booster is configured to perform binary classification.
   * @post If `_predict_contributions` is true, also calculate SHAP values for each feature and include them in the
   * output `PredictionScore`. Otherwise, only calculate the overall prediction score without SHAP values.
   * @param features Vector of feature values for a single variant.
   * @return PredictionScore containing the predicted genotype, probability score, genotype quality, variant quality,
   * and optionally SHAP values for each feature.
   */
  PredictionScore CalculateScore(const vec<f64>& features) const;

  /**
   * @brief Calculate prediction probability score for germline WGS
   * workflows.
   * @pre The LightGBM booster is configured to perform multi-class classification with the number of classes specified
   * by `_num_classes`.
   * @post If `_predict_contributions` is true, also calculate SHAP values for each feature and include them in the
   * output `PredictionScore`. Otherwise, only calculate the overall prediction score without SHAP values.
   * @param features Vector of feature values for a single variant.
   * @return PredictionScore containing the predicted genotype, probability score, genotype quality, variant quality,
   * and optionally SHAP values for each feature.
   */
  PredictionScore CalculateGermlineScore(const vec<f64>& features) const;

  /**
   * @brief Calculate prediction probability score with a minimum score threshold for germline WGS
   * workflows.
   * @pre The LightGBM booster is configured to perform multi-class classification with the number of classes specified
   * by `_num_classes`.
   * @post If `_predict_contributions` is true, also calculate SHAP values for each feature and include them in the
   * output `PredictionScore`. Otherwise, only calculate the overall prediction score without SHAP values.
   * @param features Vector of feature values for a single variant.
   * @param min_score Minimum score threshold for returning a non-reference genotype. If the highest predicted score
   * among all classes is below this threshold, the predicted genotype will be set to reference (GT=0/0) regardless of
   * the actual scores.
   * @return PredictionScore containing the predicted genotype, probability score, genotype quality, variant quality,
   * and optionally SHAP values for each feature.
   */
  PredictionScore CalculateGermlineScore(const vec<f64>& features, f64 min_score) const;

 private:
  // LightGBM booster handle
  lightgbm::BoosterPtr _booster;
  // Fast config handle for normal prediction
  lightgbm::FastConfigPtr _predict_normal_config;
  // Fast config handle for feature contribution prediction (SHAP values)
  lightgbm::FastConfigPtr _predict_contrib_config;
  // Number of classes in classification
  s32 _num_classes{};
  // Flag whether feature contribution prediction is enabled
  bool _predict_contributions{};
};

/**
 * @brief Calculate the Phred quality score based on a probability value. Phred quality
 * score is calculated as -10 * log10(1 - P), where P is the probability value.
 * @pre Probability value is between 0 and 1 inclusive.
 * @post Phred quality score is capped at 255 if the probability value is >= 1, and is 0 if the probability is <= 0.
 * @param prob Probability value
 * @return Phred quality score
 */
f64 CalculatePhredQualityScore(f64 prob);

/**
 * @brief Calculate the genotype quality score based on the prediction probability score for the predicted genotype.
 * Genotype quality score is calculated as -10 * log10(1 - P(genotype)), where P(genotype) is the prediction probability
 * score for the predicted genotype.
 * @pre Prediction probability score for the predicted genotype is between 0 and 1 inclusive.
 * @post Genotype quality score is capped at 255 if the prediction probability for the predicted genotype is 1 or very
 * close to 1, and is 0 if the prediction probability for the predicted genotype is 0 or very close to 0.
 * @param pred_prob Prediction probability score for the predicted genotype
 * @return Genotype quality score
 */
s32 CalculateGenotypeQuality(f64 pred_prob);

/**
 * @brief Calculate the variant quality score based on the prediction probabilities of all variant genotypes. The
 * variant quality score is calculated as -10 * log10(1 - P(variant)), where P(variant) is the sum of probabilities of
 * all non-reference genotypes.
 * @post Variant quality score is capped at 255 if the sum of probabilities of all non-reference genotypes is 1 or very
 * close to 1, and is 0 if the sum of probabilities of all non-reference genotypes is 0 or very close to 0.
 * @param scores Vector of prediction probabilities for all genotypes, where the probability of the reference genotype
 * is at index 0 and the probabilities of non-reference genotypes start from index 1.
 * @return Variant quality score
 */
f32 CalculateVariantQuality(const vec<f64>& scores);

}  // namespace xoos::svc
