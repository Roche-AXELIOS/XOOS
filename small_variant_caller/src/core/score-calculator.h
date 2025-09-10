#pragma once

#include <string>

#include <xoos/types/fs.h>
#include <xoos/types/vec.h>

#include "core/genotype.h"
#include "util/lightgbm-util.h"

namespace xoos::svc {

class ScoreCalculator {
 public:
  ScoreCalculator(const fs::path& model_file, size_t ncol);
  vec<std::string> GetModelFeatureNames(size_t ncol) const;
  double CalculateScore(const vec<double>& features) const;
  Genotype CalculateScoreGermline(const vec<double>& features) const;

 private:
  lightgbm::BoosterPtr _booster{};
  lightgbm::FastConfigPtr _fast_config{};
  int _num_classes{};
};

}  // namespace xoos::svc
