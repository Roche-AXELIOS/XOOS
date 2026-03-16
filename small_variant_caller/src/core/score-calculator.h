#pragma once

#include <string>

#include <xoos/types/fs.h>
#include <xoos/types/vec.h>

#include "core/genotype.h"
#include "util/lightgbm-util.h"

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
  double probability;
  // Phred score for genotype quality, calculated as -10 * log10(1 - P(genotype)), where P(genotype) is the prediction
  // probability score for the predicted genotype
  s32 genotype_quality;
  // Phred score for variant quality, calculated as -10 * log10(1 - P(variant)), where P(variant) is the sum of
  // probabilities of all non-reference genotypes
  s32 variant_quality;
};

/**
 * @brief Convert a prediction probability score to Phred quality score, which is calculated as
 * -10 * log10(1 - P(pred)), where P(pred) is the prediction probability score.
 * @pre Prediction probability score is between 0 and 1 inclusive.
 * @post Phred quality score is capped at 255 if the prediction probability is 1 or very close to 1, and is 0 if the
 * prediction probability is 0 or very close to 0.
 * @param pred_prob Prediction probability score
 * @return Phred quality score
 */
s32 CalculatePhredQualityScore(double pred_prob);

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
s32 CalculateVariantQuality(const vec<double>& scores);

class ScoreCalculator {
 public:
  ScoreCalculator(const fs::path& model_file, size_t ncol);
  vec<std::string> GetModelFeatureNames(size_t ncol) const;
  double CalculateScore(const vec<double>& features) const;
  GenotypeScore CalculateScoreGermline(const vec<double>& features) const;

 private:
  lightgbm::BoosterPtr _booster{};
  lightgbm::FastConfigPtr _fast_config{};
  int _num_classes{};
};

}  // namespace xoos::svc
