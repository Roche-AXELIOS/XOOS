#include "train-model.h"

#include <algorithm>
#include <string>
#include <unordered_map>
#include <utility>

#include <fmt/format.h>

#include <xoos/error/error.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/log/logging.h>
#include <xoos/types/str-container.h>

#ifdef SOMATIC_ENABLE
#include "core/filtering.h"
#endif  // SOMATIC_ENABLE
#include "core/variant-info-serializer.h"
#include "model-trainer.h"
#include "util/num-utils.h"
#include "util/seq-util.h"

namespace xoos::svc {

/**
 * @brief Extract the median DP for each chromosome using VCF features.
 * @param vcf_features Map of chromosome -> position -> variant ID -> VCF feature
 * @return Map of chromosome to median DP
 */
static StrUnorderedMap<u32> GetChromosomeMedianDP(ChromToVcfFeaturesMap& vcf_features) {
  // Assumptions:
  // 1. Input VCF feature file contains the `normal_dp` column.
  // 2. VCF features were not filtered in any way; there is one VCF feature for each variant in the VCF.
  StrUnorderedMap<u32> chrom_to_dp{};
  for (auto& [chr, pos_data] : vcf_features) {
    vec<u32> vals{};
    for (const auto& vid_data : std::views::values(pos_data)) {
      // Extract DP from the first variant only because variants at the same position share the same DP value
      vals.emplace_back(vid_data.begin()->second.normal_dp);
    }
    if (!vals.empty()) {
      chrom_to_dp[chr] = Median(vals);
    }
  }
  return chrom_to_dp;
}

/**
 * @brief Extract positive training data from BAM/VCF features files and truth VCF.
 * @param param Model training parameters containing input files and settings
 * @param var_bam_features Variant BAM features to be updated
 * @param ref_bam_features Reference BAM features to be updated
 * @param var_vcf_features Variant VCF features to be updated
 * @param median_dps Vector of each sample's chromosomal median DP to be updated
 * @return Number of positive data points and number of samples
 */
static std::pair<size_t, size_t> GetPositiveData(const TrainModelParam& param,
                                                 ChromToVariantInfoMapWithLabel& var_bam_features,
                                                 RefInfoMapMultiSample& ref_bam_features,
                                                 ChromToVcfFeaturesMapMultiSample& var_vcf_features,
                                                 vec<StrUnorderedMap<u32>>& median_dps) {
  // For each sample:
  // 1. Load VCF features if available and extract the median DP for each chromosome from the VCF features.
  // 2. Load variant and reference BAM features.
  // 3. For each variant, check if it is a true positive or false positive based on the truth VCF.
  // 4. Store the variant features and labels in the output maps.

  size_t num_positives = 0;
  const size_t num_samples = param.truth_vcfs.size();
  const bool has_vcf_features = !param.positive_vcf_features.empty();
  for (u32 sid = 0; sid < num_samples; ++sid) {
    ChromToVcfFeaturesMap vcf_features;
    if (has_vcf_features) {
      const auto& vcf_feat_file = param.positive_vcf_features[sid];
      Logging::Info("Loading positive VCF features from {}", vcf_feat_file);
      vcf_features = VariantInfoSerializer::LoadVcfFeatures(vcf_feat_file);
      // Use VCF features to extract the median DP for each chromosome in this sample
      auto chrom_to_dps = GetChromosomeMedianDP(vcf_features);
      median_dps.emplace_back(chrom_to_dps);
      // Log median DP for each chromosome in this sample
      vec<std::string> items;
      items.reserve(chrom_to_dps.size());
      for (auto& [chrom, dp] : chrom_to_dps) {
        items.emplace_back(fmt::format("{}:{}", chrom, dp));
      }
      std::ranges::sort(items);
    }
    ChromToVariantInfoMap var_features;
    RefInfoMap ref_features;
    TruthFeatures truth_features;
    {
      const auto& bam_feat_file = param.positive_features[sid];
      const auto& truth_vcf = param.truth_vcfs[sid];
      Logging::Info("Loading positive BAM features from {}", bam_feat_file);
#ifdef SOMATIC_ENABLE
      if (param.workflow == Workflow::kSomatic) {
        const auto& [var_infos, ref_infos] = VariantInfoSerializer::LoadFeatures(bam_feat_file);
        ref_features = ref_infos;
        Logging::Info("Integrating truth VCF {}", truth_vcf);
        // Remove BAM features that are not in the truth VCF
        var_features = FilterVariantsByVcf(truth_vcf, var_infos);
      }
#endif  // SOMATIC_ENABLE
      if (param.workflow == Workflow::kGermline || param.workflow == Workflow::kGermlineMultiSample) {
        // use VCF feature positions to filter BAM features
        const auto& [var_infos, ref_infos] = VariantInfoSerializer::LoadFeatures(bam_feat_file, vcf_features);
        ref_features = ref_infos;
        var_features = var_infos;
        Logging::Info("Integrating truth VCF {}", truth_vcf);
        truth_features = GetTruthFeatures(truth_vcf);
      } else {
        throw error::Error("Unsupported workflow");
      }
    }
    for (auto& [chr, pos_data] : var_features) {
      for (auto& [pos, vid_data] : pos_data) {
        bool has_snv_ins = false;
        bool has_del = false;
        if (param.workflow == Workflow::kGermline || param.workflow == Workflow::kGermlineMultiSample) {
          for (auto& [vid, var_feat] : vid_data) {
            if (!vcf_features[chr][pos].contains(vid)) {
              continue;
            }
            if (vid.type == VariantType::kDeletion) {
              has_del = true;
            } else {
              has_snv_ins = true;
            }
            auto& vcf_feature = vcf_features[chr][pos][vid];
            const auto& [gt, ref_alt_set] = truth_features[chr][pos];
            if (gt != Genotype::kGTNA) {
              // Unsupported genotype should NEVER be treated as a negative or false positive in the training data!
              // They should be excluded from training data because the true genotype is uncertain.
              // Otherwise, they are introducing noise in the training data.
              const auto& [ref_trimmed, alt_trimmed] = TrimVariant(vid.ref, vid.alt);
              std::string label{"0"};
              if (gt != Genotype::kGT00 && gt != Genotype::kGT0 &&
                  ref_alt_set.contains(std::make_pair(ref_trimmed, alt_trimmed))) {
                label = GenotypeToIntString(gt);
                num_positives++;
              }
              var_bam_features[chr][pos][vid].emplace_back(var_feat, label, sid);
              var_vcf_features[chr][pos][vid].emplace_back(vcf_feature, sid);
            }
          }
        } else {
#ifdef SOMATIC_ENABLE
          for (auto& [vid, var_feat] : vid_data) {
            auto& vcf_feature = vcf_features[chr][pos][vid];
            var_bam_features[chr][pos][vid].emplace_back(var_feat, "1", sid);
            var_vcf_features[chr][pos][vid].emplace_back(vcf_feature, sid);
            num_positives++;
          }
#endif  // SOMATIC_ENABLE
        }
        if (has_snv_ins) {
          ref_bam_features[chr][sid][pos] = ref_features[chr][pos];
        }
        if (has_del) {
          // deletion reference feature position needs to be offset by 1
          ref_bam_features[chr][sid][pos + 1] = ref_features[chr][pos + 1];
        }
      }
    }
  }
  return std::make_pair(num_positives, num_samples);
}

#ifdef SOMATIC_ENABLE
/**
 * @brief Extract negative training data from BAM/VCF features files.
 * @param param Model training parameters containing input files and settings
 * @param var_bam_features Variant BAM features to be updated
 * @param ref_bam_features Reference BAM features to be updated
 * @param var_vcf_features Variant VCF features to be updated
 * @param num_samples Number of samples
 * @return Number of negative data points
 */
size_t GetNegativeData(const TrainModelParam& param,
                       ChromToVariantInfoMapWithLabel& var_bam_features,
                       RefInfoMapMultiSample& ref_bam_features,
                       ChromToVcfFeaturesMapMultiSample& var_vcf_features,
                       const u32 num_samples) {
  // This function is currently only relevant to the `somatic` workflow
  size_t num_negatives = 0;
  for (u32 i = 0; i < param.negative_features.size(); i++) {
    auto& bam_feat_file = param.negative_features[i];
    Logging::Info("Loading healthy(negative) features from {}", bam_feat_file);
    auto [var_info_map, ref_info_map] = VariantInfoSerializer::LoadFeatures(bam_feat_file);
    ChromToVcfFeaturesMap vcf_features;
    if (!param.negative_vcf_features.empty()) {
      auto& vcf_feat_file = param.negative_vcf_features[i];
      vcf_features = VariantInfoSerializer::LoadVcfFeatures(vcf_feat_file);
    }
    // the negative variants are those in the healthy samples with a weighted score < 4
    for (const auto& [chr, pos_data] : var_info_map) {
      for (auto& [pos, vid_data] : pos_data) {
        // TODO : Update how blocklist filtering is done to speed up feature parsing
        bool skip = false;
        for (auto& region : param.blocklist) {
          const u64 start = region.start >= 0 ? static_cast<u64>(region.start) : 0;
          const u64 end = region.end >= 0 ? static_cast<u64>(region.end) : 0;
          if (region.chromosome == chr && start <= pos && pos < end) {
            skip = true;
            break;
          }
        }
        if (skip) {
          continue;
        }
        auto sid = num_samples + i;
        ref_bam_features[chr][sid][pos] = ref_info_map[chr][pos];
        for (auto& [vid, var_feat] : vid_data) {
          // Score based filtering of negative features done to only select negative data points that are noise
          // True variants may also occur within the negatives that need to be excluded via this filter
          if (var_feat.weighted_score < param.max_score) {
            auto& vcf_feat = vcf_features[chr][pos][vid];
            var_bam_features[chr][pos][vid].emplace_back(var_feat, "0", sid);
            var_vcf_features[chr][pos][vid].emplace_back(vcf_feat, sid);
            num_negatives++;
          }
        }
      }
    }
  }
  return num_negatives;
}
#endif  // SOMATIC_ENABLE

/**
 * @brief Trains germline SNV and Indel models.
 * @param param CLI parameters and input/output file paths
 * @param model_config Configuration from JSON config file
 */
static void TrainGermline(const TrainModelParam& param, SVCConfig& model_config) {
  // Read in truth vcfs and positive variants and check if features files have matching feature columns
  ChromToVariantInfoMapWithLabel data;
  RefInfoMapMultiSample ref_data;
  ChromToVcfFeaturesMapMultiSample vcf_data;
  vec<StrUnorderedMap<u32>> median_dps;
  auto [num_true_pos, sample_count] = GetPositiveData(param, data, ref_data, vcf_data, median_dps);
  if (num_true_pos == 0) {
    throw error::Error("No true positives in training data");
  }
  if (param.normalize_features && median_dps.size() != sample_count) {
    throw error::Error(
        "Number of normalize target(s) {} and sample count {} are not equal", median_dps.size(), sample_count);
  }
  // Same variant of a genotype repeated in multiple samples may lead to over-represenation in the training data.
  // This may lead to the model converging too early, which yields a lower than expected F1 score. To avoid this, we
  // reduce redundancy by keeping only one sample for each variant with the same genotype.
  const bool reduce_redundancy = sample_count > 1;
  ModelTrainer trainer(data, ref_data, vcf_data, model_config, param.workflow, median_dps);
  const std::string version = param.command_line.has_value() ? param.command_line->version : "unknown";
  trainer.Train(model_config.indel_model_file,
                param.threads,
                param.indel_iterations,
                VariantGroup::kIndelOnly,
                reduce_redundancy,
                param.write_training_data_tsv);
  trainer.Train(model_config.snv_model_file,
                param.threads,
                param.snv_iterations,
                VariantGroup::kSnvOnly,
                reduce_redundancy,
                param.write_training_data_tsv);
}

/**
 * @brief Train model using the specified parameters
 * @param param Model training parameters and input/output file paths.
 */
void TrainModel(const TrainModelParam& param) {
  SVCConfig model_config = param.config;

  if (param.workflow == Workflow::kGermline || param.workflow == Workflow::kGermlineMultiSample) {
    // Perform germline-specific model training
    if (!param.output_file.empty()) {
      if (param.output_file.size() == 2) {
        model_config.GetGermlineModelPaths(param.output_file);
      } else if (model_config.indel_model_file.empty() || model_config.snv_model_file.empty()) {
        throw error::Error("Please specify two output model file paths or set SNV and Indel models in the config");
      }
    }
    TrainGermline(param, model_config);
    return;
  }
#ifdef SOMATIC_ENABLE
  if (param.output_file.size() != 1) {
    throw error::Error("Please specify one output model path");
  }
  const auto output_model_file = param.output_file[0];

  // Load BAM/VCF feature files and integrate with truth VCFs
  ChromToVariantInfoMapWithLabel var_bam_features;
  RefInfoMapMultiSample ref_bam_features;
  ChromToVcfFeaturesMapMultiSample vcf_features;
  u32 num_samples = 0;
  vec<StrUnorderedMap<u32>> median_dps;
  if (param.workflow == Workflow::kSomatic) {
    auto [num_positives, sample_count] =
        GetPositiveData(param, var_bam_features, ref_bam_features, vcf_features, median_dps);
    if (num_positives == 0) {
      throw error::Error("No positive variants found");
    }
    auto num_negatives = GetNegativeData(param, var_bam_features, ref_bam_features, vcf_features, sample_count);
    if (num_negatives == 0) {
      throw error::Error("No negative variants found");
    }
    num_samples = sample_count;
  } else {
    auto [num_data_points, sample_count] =
        GetPositiveData(param, var_bam_features, ref_bam_features, vcf_features, median_dps);
    if (num_data_points == 0) {
      throw error::Error("No training variants found");
    }
    num_samples = sample_count;
  }

  if (param.normalize_features && median_dps.size() != num_samples) {
    throw error::Error(
        "Number of normalize target(s) {} and sample count {} are not equal", median_dps.size(), num_samples);
  }
  // Train model
  ModelTrainer trainer(var_bam_features, ref_bam_features, vcf_features, model_config, param.workflow, median_dps);
  trainer.Train(output_model_file, param.threads, param.iterations);
#endif  // SOMATIC_ENABLE
}

}  // namespace xoos::svc
