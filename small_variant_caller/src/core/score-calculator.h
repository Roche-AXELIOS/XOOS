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
struct GenotypeScore {
  // predicted genotype based on prediction probability score
  Genotype genotype;
  // prediction probability score for the predicted genotype
  f64 probability;
  // Phred score for genotype quality, calculated as -10 * log10(1 - P(genotype)), where P(genotype) is the prediction
  // probability score for the predicted genotype
  s32 genotype_quality;
  // Phred score for variant quality, calculated as -10 * log10(1 - P(variant)), where P(variant) is the sum of
  // probabilities of all non-reference genotypes
  f32 variant_quality;
};

class ScoreCalculator {
 public:
  ScoreCalculator(const fs::path& model_file, size_t ncol, const std::string& prediction_params);
  vec<std::string> GetModelFeatureNames() const;
  f64 CalculateScore(const vec<f64>& features) const;
  GenotypeScore CalculateScoreGermline(const vec<f64>& features) const;
  GenotypeScore CalculateScoreGermline(const vec<f64>& features, f64 min_score) const;

 private:
  lightgbm::BoosterPtr _booster{};
  lightgbm::FastConfigPtr _fast_config{};
  s32 _num_classes{};
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
