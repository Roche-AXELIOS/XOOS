#pragma once

#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "core/command-line-info.h"
#include "core/config.h"
#include "core/workflow.h"
#include "training-data-set.h"
#include "util/region-util.h"

namespace xoos::svc {

/**
 * @brief Parameters for training a variant calling model.
 */
struct TrainModelParam {
  vec<fs::path> positive_features{};
  vec<fs::path> positive_vcf_features{};
  vec<fs::path> negative_features{};
  vec<fs::path> negative_vcf_features{};
  vec<fs::path> truth_vcfs{};
  fs::path output_file{};
  fs::path snv_output_file{};
  fs::path indel_output_file{};
  SVCConfig config{};
  std::optional<fs::path> config_file;
  std::optional<StrMap<vec<Interval>>> blocklist{};
  Workflow workflow{};
  u32 max_score{};
  u32 iterations{};
  u32 snv_iterations{};
  u32 indel_iterations{};
  size_t threads{};
  FeatureNormalization normalize_features{FeatureNormalization::kNone};
  std::optional<fs::path> output_training_data{};
  std::optional<fs::path> snv_output_training_data{};
  std::optional<fs::path> indel_output_training_data{};
  std::optional<CommandLineInfo> command_line;
};

void TrainModel(const TrainModelParam& param);

class TrainingController {
 public:
  TrainingController(TrainModelParam param, SVCConfig config);

  /**
   * @brief Validate that VCF features files are compatible with feature normalization settings.
   * @pre If feature normalization is set to `median-dp`, VCF features files must contain depth information
   * for both tumor and normal samples (if applicable).
   * @throws error::Error if any VCF features file is incompatible with the normalization settings.
   */
  void ValidateVcfFeaturesFilesForNormalization() const;

  /**
   * @brief Extract chromosome median depth values from VCF features for feature normalization.
   * @pre VCF features files must be compatible with the normalization settings, as validated by
   * `ValidateVcfFeaturesFilesForNormalization()`.
   * @post A ChromMedianDepth struct is returned containing the median depth values for each chromosome.
   * If a chromosome is not found in the VCF features, its depth values in the returned struct will be zero.
   * @param vcf_feats Map of VariantId to VcfFeature containing the extracted VCF features.
   * @return ChromMedianDepth struct containing median depth values for each chromosome.
   */
  ChromMedianDepth GetChromMedianDepthFromVcfFeatures(const VarIdToVcfFeatures& vcf_feats) const;

  /**
   * @brief Train SNV and Indel models for `germline` or `germline-multi-sample` workflow.
   * @pre The number of positive BAM feature files, positive VCF feature files, and truth VCF files must be the same.
   * @pre The order of positive BAM feature files, positive VCF feature files, and truth VCF files must correspond to
   * the same set of samples. For example, `positive_features[0]`, `positive_vcf_features[0]`, and `truth_vcfs[0]`
   * should all correspond to the same sample.
   * @pre There must be at least one positive sample in the training data.
   * @post Two model files are created at the specified output paths.
   */
  void TrainGermline();

  /**
   * @brief Train a model for `germline-tagging` workflow.
   * @pre The number of positive BAM feature files, positive VCF feature files, and truth VCF files must be the same.
   * @pre The order of positive BAM feature files, positive VCF feature files, and truth VCF files must correspond to
   * the same set of samples. For example, `positive_features[0]`, `positive_vcf_features[0]`, and `truth_vcfs[0]`
   * should all correspond to the same sample.
   * @pre There must be at least one positive sample in the training data.
   * @post A model file is created at the specified output path.
   */
  void TrainGermlineTagging();

  /**
   * @brief Train a model for `tumor-only-te` workflow.
   * @pre The number of positive BAM feature files and truth VCF files must be the same.
   * @pre The order of positive BAM feature files  and truth VCF files must correspond to
   * the same set of samples. For example, `positive_features[0]` and `truth_vcfs[0]`
   * should all correspond to the same sample.
   * @pre There must be at least one positive sample and one negative sample in the training data.
   * @note The number of negative samples does not need to match the number of positive samples.
   * @post A model file is created at the specified output path.
   */
  void TrainTumorOnlyTe();

  /**
   * @brief Train a model for `tumor-normal-wgs` workflow.
   * @pre The number of positive BAM feature files, positive VCF feature files, and truth VCF files must be the same.
   * @pre The order of positive BAM feature files, positive VCF feature files, and truth VCF files must correspond to
   * the same set of samples. For example, `positive_features[0]`, `positive_vcf_features[0]`, and `truth_vcfs[0]`
   * should all correspond to the same sample.
   * @pre The number of negative BAM feature files and negative VCF feature files must be the same.
   * @pre The order of negative BAM feature files and negative VCF feature files must correspond to the same set of
   * samples. For example, `negative_features[0]` and `negative_vcf_features[0]` should both correspond to the same
   * sample.
   * @pre There must be at least one positive sample and one negative sample in the training data.
   * @note The number of negative samples does not need to match the number of positive samples.
   * @post A model file is created at the specified output path.
   */
  void TrainTumorNormalWgs();

 private:
  /**
   * @brief Extract positive training data in `tumor-only-te` workflow for the specified sample index.
   * @post Training data is extracted from the BAM features.
   * @post Variants not found in the truth VCF are excluded from the training data.
   * @param sid Sample index for which to extract positive training data from the input files.
   * @param columns Feature columns to use for extracting features from the input files.
   * @param data TrainingDataSet2 object to insert the extracted training data into.
   */
  void GetPositiveDataForTumorOnlyTe(size_t sid, const vec<FeatureColumn>& columns, TrainingDataSet& data);

  /**
   * @brief Extract positive training data in `germline` or `germline-multi-sample` workflow for the specified sample
   * index.
   * @post Training data is extracted for variants having both BAM features and VCF features.
   * @post Variants with GT=0/1, GT=1/1, or GT=1/2 in the truth VCF are assigned separate positive labels.
   * @post Variants not found in the truth VCF are assigned the negative label.
   * @post Variants with unsupported genotypes in the truth VCF will be excluded from the training data.
   * @post If feature normalization is enabled in the configuration, features will be normalized using the median depth
   * values calculated from the VCF features.
   * @post SNVs and indels are extracted separately from the input files, and inserted into separate TrainingDataSet2
   * objects for training separate SNV and indel models.
   * @param sid Sample index for which to extract positive training data from the input files.
   * @param snv_columns Feature columns to use for extracting features for SNVs from the input files.
   * @param indel_columns Feature columns to use for extracting features for indels from the input files.
   * @param snv_data TrainingDataSet2 object to insert the extracted SNV training data into.
   * @param indel_data TrainingDataSet2 object to insert the extracted indel training data into.
   */
  void GetPositiveDataForGermline(size_t sid,
                                  const vec<FeatureColumn>& snv_columns,
                                  const vec<FeatureColumn>& indel_columns,
                                  TrainingDataSet& snv_data,
                                  TrainingDataSet& indel_data);

  /**
   * @brief Extract positive training data in `tumor-normal-wgs` workflow for the specified sample index.
   * @post Training data is extracted for variants having both BAM features and VCF features.
   * @post Variants found in the truth VCF are assigned the default positive label.
   * @post Variants not found in the truth VCF are excluded from the training data.
   * @post Variants with insufficient tumor support in the BAM features are excluded from the training data.
   * @post If feature normalization is enabled in the configuration, features will be normalized using the median depth
   * values calculated from the VCF features.
   * @param sid Sample index for which to extract positive training data from the input files.
   * @param columns Feature columns to use for extracting features from the input files.
   * @param data TrainingDataSet2 object to insert the extracted training data into.
   */
  void GetPositiveDataForTumorNormalWgs(size_t sid, const vec<FeatureColumn>& columns, TrainingDataSet& data);

  /**
   * @brief Extract positive training data in `germline-tagging` workflow for the specified sample index.
   * @post Training data is extracted for variants having both BAM features and VCF features.
   * @post Variants with GT=0/1 or GT=1/1 in the truth VCF are assigned separate positive labels.
   * @post Variants not found in the truth VCF are assigned the negative label.
   * @post Variants with unknown genotype or GT=1/2 in the truth VCF are excluded from the training data.
   * @param sid Sample index for which to extract positive training data from the input files.
   * @param columns Feature columns to use for extracting features from the input files.
   * @param data TrainingDataSet2 object to insert the extracted training data into.
   */
  void GetPositiveDataForGermlineTagging(size_t sid, const vec<FeatureColumn>& columns, TrainingDataSet& data);

  /**
   * @brief Extract negative training data in `tumor-only-te` workflow for the specified sample index.
   * @post Training data is extracted from BAM features.
   * @post Training data is assigned the default negative label.
   * @post Training data may be filtered using the blocklist regions specified in the configuration, if provided.
   * @post Training data with weighted score greater than the maximum score specified in the configuration will be
   * excluded.
   * @param sid Sample index for which to extract negative training data from the input files.
   * @param columns Feature columns to use for extracting features from the input files.
   * @param data TrainingDataSet2 object to insert the extracted training data into.
   * @see kDefaultNegativeLabel
   */
  void GetNegativeDataForTumorOnlyTe(size_t sid, const vec<FeatureColumn>& columns, TrainingDataSet& data);

  /**
   * @brief Extract negative training data in `tumor-normal-wgs` workflow for the specified sample index.
   * @post Training data is extracted for variants having both BAM features and VCF features.
   * @post Training data is assigned the default negative label.
   * @post If feature normalization is enabled in the configuration, features will be normalized using the median depth
   * values calculated from the VCF features.
   * @param sid Sample index for which to extract negative training data from the input files.
   * @param columns Feature columns to use for extracting features from the input files.
   * @param data TrainingDataSet2 object to insert the extracted training data into.
   * @see kDefaultNegativeLabel
   */
  void GetNegativeDataForTumorNormalWgs(size_t sid, const vec<FeatureColumn>& columns, TrainingDataSet& data);

  const TrainModelParam _param;
  SVCConfig _config;
  const Workflow _workflow;
};

}  // namespace xoos::svc
