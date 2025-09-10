#pragma once
#include <set>
#include <string>
#include <unordered_map>
#include <utility>

#include <util/lightgbm-util.h>

#include <xoos/io/fasta-reader.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>

#include "core/config.h"
#include "core/variant-info.h"
#include "core/workflow.h"

namespace xoos::svc {

using TruthFeatures = StrUnorderedMap<                                              // chrom
    std::unordered_map<u32,                                                         // pos
                       std::pair<Genotype,                                          // GT
                                 std::set<std::pair<std::string, std::string>>>>>;  // REF, ALT

TruthFeatures GetTruthFeatures(const fs::path& vcf);

// Specify which group of variants to train on.
enum class VariantGroup {
  kAll,
  kSnvOnly,
  kIndelOnly
};

class ModelTrainer {
 public:
  ModelTrainer(ChromToVariantInfoMapWithLabel& features,
               RefInfoMapMultiSample& ref_features,
               ChromToVcfFeaturesMapMultiSample& vcf_features,
               SVCConfig& config,
               Workflow workflow,
               vec<StrUnorderedMap<u32>>& normalize_targets);

  void Train(const fs::path& output_file,
             size_t num_threads = 1,
             uint32_t num_rounds = kNumTrainingRounds,
             VariantGroup var_group = VariantGroup::kAll,
             bool reduce_redundancy = false,
             bool write_data = false);

 private:
  const ChromToVariantInfoMapWithLabel& _features{};
  const RefInfoMapMultiSample _ref_features{};
  const ChromToVcfFeaturesMapMultiSample _vcf_features{};
  SVCConfig _config{};
  Workflow _workflow{Workflow::kGermline};
  vec<StrUnorderedMap<u32>> _normalize_targets;
  lightgbm::DatasetPtr _train_data{};
  lightgbm::BoosterPtr _booster{};

  static constexpr u32 kNumTrainingRounds = 10000;
};

}  // namespace xoos::svc
