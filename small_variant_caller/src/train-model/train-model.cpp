#include "train-model.h"

#include <string>
#include <utility>

#include <fmt/format.h>

#include <xoos/error/error.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/log/logging.h>
#include <xoos/types/str-container.h>

#include "core/variant-feature-extraction.h"
#include "core/variant-info-serializer.h"
#include "model-trainer.h"
#include "truth-label-set.h"
#include "util/file-util.h"
#include "util/seq-util.h"

namespace xoos::svc {

TrainingController::TrainingController(TrainModelParam param, SVCConfig config)
    : _param(std::move(param)), _config(std::move(config)), _workflow(_param.workflow) {
  if (_param.normalize_features == FeatureNormalization::kMedianDp) {
    // Ensure that VCF features files are compatible with median depth normalization.
    ValidateVcfFeaturesFilesForNormalization();
  }
}

/**
 * @brief Extract a set of VariantIds from a VCF file.
 * @param vcf_path Path to the VCF file.
 * @return Set of VariantIds extracted from the VCF.
 */
static std::unordered_set<VariantId> GetVariantIdSet(const fs::path& vcf_path) {
  std::unordered_set<VariantId> var_ids;
  const io::VcfReader reader(vcf_path);
  io::VcfRecordPtr record;
  while ((record = reader.GetNextRecord())) {
    const auto& ref = record->Allele(0);
    if (!ContainsOnlyACTG(ref)) {
      continue;
    }
    const auto& chrom = record->Chromosome();
    const auto pos = static_cast<u64>(record->Position());
    for (auto i = 1; i < record->NumAlleles(); ++i) {
      const auto& alt = record->Allele(i);
      if (!ContainsOnlyACTG(alt)) {
        continue;
      }
      const auto& [ref_trimmed, alt_trimmed] = TrimVariant(ref, alt);
      var_ids.emplace(chrom, pos, ref_trimmed, alt_trimmed);
    }
  }
  return var_ids;
}

void TrainingController::ValidateVcfFeaturesFilesForNormalization() const {
  const bool has_sample_context = FeatureNamesHaveSampleContext(_workflow);
  for (const auto& vcf_feat_file : _param.positive_vcf_features) {
    VariantInfoSerializer::ValidateFeatureFileHeaderForNormalization(vcf_feat_file, has_sample_context);
  }
  for (const auto& vcf_feat_file : _param.negative_vcf_features) {
    VariantInfoSerializer::ValidateFeatureFileHeaderForNormalization(vcf_feat_file, has_sample_context);
  }
}

ChromMedianDepth TrainingController::GetChromMedianDepthFromVcfFeatures(const VarIdToVcfFeatures& vcf_feats) const {
  ChromMedianDepth result{};
  if (_param.normalize_features == FeatureNormalization::kMedianDp) {
    // Extract chromosome median depth values from VCF features
    result = GetChromosomeMedianDepth(vcf_feats);
    // Check whether the calculated median depth values are valid for normalization.
    ValidateChromMedianDepthForNormalization(result, FeatureNamesHaveSampleContext(_workflow));
  }
  return result;
}

void TrainingController::GetPositiveDataForTumorOnlyTe(const size_t sid,
                                                       const vec<FeatureColumn>& columns,
                                                       TrainingDataSet& data) {
  const auto& truth_vcf = _param.truth_vcfs.at(sid);
  Logging::Info("Loading truth VCF from {}", truth_vcf);
  const auto truth_vids = GetVariantIdSet(truth_vcf);

  const auto& bam_feat_file = _param.positive_features.at(sid);
  Logging::Info("Loading positive BAM features from {}", bam_feat_file);
  const auto generator = VariantInfoSerializer::BamFeatureTupleGenerator(bam_feat_file);
  while (auto record = generator()) {
    constexpr DepthTuple kZeroDepth{};
    const auto& [vid, bam_feat] = record.value();
    if (!truth_vids.empty() && !truth_vids.contains(vid)) {
      continue;
    }
    // Use zero depth here because we cannot normalize features without VCF features
    const auto feat_vec = GetFeatureVec(columns, vid, bam_feat, kZeroVcfFeature, kZeroDepth);
    data.InsertRow(vid, feat_vec, kDefaultPositiveLabel);
  }
}

void TrainingController::GetPositiveDataForGermline(const size_t sid,
                                                    const vec<FeatureColumn>& snv_columns,
                                                    const vec<FeatureColumn>& indel_columns,
                                                    TrainingDataSet& snv_data,
                                                    TrainingDataSet& indel_data) {
  const auto& vcf_feat_file = _param.positive_vcf_features.at(sid);
  Logging::Info("Loading positive VCF features from {}", vcf_feat_file);
  const auto& vcf_feats = VariantInfoSerializer::LoadVcfFeatures(vcf_feat_file, false);

  const auto& truth_vcf = _param.truth_vcfs.at(sid);
  Logging::Info("Loading truth VCF from {}", truth_vcf);
  const auto truth_genotypes = GetGenotypeToVariantIds(truth_vcf);

  const ChromMedianDepth median_depths = GetChromMedianDepthFromVcfFeatures(vcf_feats);

  const auto& bam_feat_file = _param.positive_features.at(sid);
  Logging::Info("Loading positive BAM features from {}", bam_feat_file);
  const auto generator = VariantInfoSerializer::BamFeatureTupleGenerator(bam_feat_file);
  while (auto record = generator()) {
    const auto& [vid, bam_feat] = record.value();
    if (!vcf_feats.empty() && !vcf_feats.contains(vid)) {
      continue;
    }
    const auto genotype = GetGenotypeForVariant(vid, truth_genotypes).value_or(Genotype::kGT00);
    if (genotype == Genotype::kGTNA) {
      // GT=NA is an unsupported genotype.
      // Unsupported genotypes should NEVER be treated as a negative or false positive in the training data!
      // They should be excluded from training data because the true genotype is uncertain.
      // Otherwise, they are introducing noise in the training data.
      continue;
    }
    const auto label = GenotypeToInt(genotype);
    const auto& vcf_feat = vcf_feats.at(vid);
    const auto& depth = median_depths.GetValue(vid.chrom);
    if (CheckVariantType(VariantGroup::kSnvOnly, vid.type)) {
      const auto feat_vec = GetFeatureVec(snv_columns, vid, bam_feat, vcf_feat, depth);
      snv_data.InsertRow(vid, feat_vec, label);
    } else {
      const auto feat_vec = GetFeatureVec(indel_columns, vid, bam_feat, vcf_feat, depth);
      indel_data.InsertRow(vid, feat_vec, label);
    }
  }
}

void TrainingController::GetPositiveDataForTumorNormalWgs(const size_t sid,
                                                          const vec<FeatureColumn>& columns,
                                                          TrainingDataSet& data) {
  const auto& truth_vcf = _param.truth_vcfs.at(sid);
  Logging::Info("Loading truth VCF from {}", truth_vcf);
  const auto truth_vids = GetVariantIdSet(truth_vcf);
  if (truth_vids.empty()) {
    throw error::Error("No valid variants found in truth VCF: {}\nPlease check your input truth VCF file(s).",
                       truth_vcf);
  }

  // Load VCF features
  const auto& vcf_feat_file = _param.positive_vcf_features.at(sid);
  Logging::Info("Loading positive VCF features from {}", vcf_feat_file);
  const auto& vcf_feats = VariantInfoSerializer::LoadVcfFeatures(vcf_feat_file, true);
  if (vcf_feats.empty()) {
    throw error::Error(
        "No valid VCF features found in positive VCF features file: {}\nPlease check your input VCF features file(s).",
        vcf_feat_file);
  }

  const ChromMedianDepth median_depths = GetChromMedianDepthFromVcfFeatures(vcf_feats);

  const auto& bam_feat_file = _param.positive_features.at(sid);
  Logging::Info("Loading positive BAM features from {}", bam_feat_file);
  const auto generator = VariantInfoSerializer::TumorNormalBamFeatureTupleGenerator(bam_feat_file);
  while (auto record = generator()) {
    const auto& [vid, bam_feat] = record.value();
    // Only include variants that have truth label, VCF feature, and BAM feature available, to ensure the quality of
    // training data.
    if (!truth_vids.contains(vid) || !vcf_feats.contains(vid)) {
      continue;
    }
    // skip variants that do not have sufficient support in the tumor sample
    if (bam_feat.tumor_var_feat.support < _config.min_tumor_support) {
      continue;
    }
    const auto& vcf_feat = vcf_feats.at(vid);
    const auto& depth = median_depths.GetValue(vid.chrom);
    const auto feat_vec = GetFeatureVec(columns, vid, bam_feat, vcf_feat, depth);
    data.InsertRow(vid, feat_vec, kDefaultPositiveLabel);
  }
}

void TrainingController::GetPositiveDataForGermlineTagging(const size_t sid,
                                                           const vec<FeatureColumn>& columns,
                                                           TrainingDataSet& data) {
  const auto& vcf_feat_file = _param.positive_vcf_features.at(sid);
  Logging::Info("Loading positive VCF features from {}", vcf_feat_file);
  const auto& vcf_feats = VariantInfoSerializer::LoadVcfFeatures(vcf_feat_file, true);

  const ChromMedianDepth median_depths = GetChromMedianDepthFromVcfFeatures(vcf_feats);

  const auto& truth_vcf = _param.truth_vcfs.at(sid);
  Logging::Info("Loading truth VCF from {}", truth_vcf);
  const auto truth_genotypes = GetGenotypeToVariantIds(truth_vcf);

  const auto& bam_feat_file = _param.positive_features.at(sid);
  Logging::Info("Loading positive BAM features from {}", bam_feat_file);
  const auto generator = VariantInfoSerializer::TumorNormalBamFeatureTupleGenerator(bam_feat_file);
  while (auto record = generator()) {
    const auto& [vid, bam_feat] = record.value();
    if (!vcf_feats.empty() && !vcf_feats.contains(vid)) {
      continue;
    }
    const auto genotype = GetGenotypeForVariant(vid, truth_genotypes).value_or(Genotype::kGT00);
    if (genotype == Genotype::kGTNA || genotype == Genotype::kGT12) {
      // GT=NA and GT=1/2 (multi-allelic with 2 alt alleles) are unsupported genotypes.
      continue;
    }
    const auto label = GenotypeToInt(genotype);
    const auto& vcf_feat = vcf_feats.at(vid);
    const auto& depth = median_depths.GetValue(vid.chrom);
    const auto feat_vec = GetFeatureVec(columns, vid, bam_feat, vcf_feat, depth);
    data.InsertRow(vid, feat_vec, label);
  }
}

static bool IsVariantAtIntervals(const VariantId& vid, const ChromIntervalsMap& chrom_intervals) {
  if (!chrom_intervals.contains(vid.chrom)) {
    return false;
  }
  const auto& intervals = chrom_intervals.at(vid.chrom);
  for (const auto& [start, end] : intervals) {
    if (start <= vid.pos && vid.pos < end) {
      return true;
    }
    if (vid.pos < start) {
      // Intervals are sorted, no need to check further
      break;
    }
  }
  return false;
}

void TrainingController::GetNegativeDataForTumorOnlyTe(const size_t sid,
                                                       const vec<FeatureColumn>& columns,
                                                       TrainingDataSet& data) {
  const auto& bam_feat_file = _param.negative_features.at(sid);
  Logging::Info("Loading negative BAM features from {}", bam_feat_file);
  const auto generator = VariantInfoSerializer::BamFeatureTupleGenerator(bam_feat_file);
  while (auto record = generator()) {
    constexpr DepthTuple kZeroDepth{};
    const auto& [vid, bam_feat] = record.value();
    if (_param.blocklist.has_value() && IsVariantAtIntervals(vid, _param.blocklist.value())) {
      continue;
    }
    if (bam_feat.var_feat.weighted_score >= _param.max_score) {
      continue;
    }
    const auto feat_vec = GetFeatureVec(columns, vid, bam_feat, kZeroVcfFeature, kZeroDepth);
    data.InsertRow(vid, feat_vec, kDefaultNegativeLabel);
  }
}

void TrainingController::GetNegativeDataForTumorNormalWgs(const size_t sid,
                                                          const vec<FeatureColumn>& columns,
                                                          TrainingDataSet& data) {
  const auto& vcf_feat_file = _param.negative_vcf_features.at(sid);
  Logging::Info("Loading negative VCF features from {}", vcf_feat_file);
  const auto& vcf_feats = VariantInfoSerializer::LoadVcfFeatures(vcf_feat_file, true);

  const ChromMedianDepth median_depths = GetChromMedianDepthFromVcfFeatures(vcf_feats);

  const auto& bam_feat_file = _param.negative_features.at(sid);
  Logging::Info("Loading negative BAM features from {}", bam_feat_file);
  const auto generator = VariantInfoSerializer::TumorNormalBamFeatureTupleGenerator(bam_feat_file);
  while (auto record = generator()) {
    const auto& [vid, bam_feat] = record.value();
    if (!vcf_feats.empty() && !vcf_feats.contains(vid)) {
      continue;
    }
    const auto& vcf_feat = vcf_feats.at(vid);
    const auto& depth = median_depths.GetValue(vid.chrom);
    const auto feat_vec = GetFeatureVec(columns, vid, bam_feat, vcf_feat, depth);
    data.InsertRow(vid, feat_vec, kDefaultNegativeLabel);
  }
}

/**
 * @brief Verify that the output file path is not empty and create parent directory if it does not exist.
 * @param output_file Path of the output file
 */
static void VerifyOutputFilePath(const fs::path& output_file) {
  if (output_file.empty()) {
    throw error::Error("Please specify output model file path");
  }
  CreateParentDirectoryIfNotExists(output_file);
}

void TrainingController::TrainTumorNormalWgs() {
  VerifyOutputFilePath(_param.output_file);

  TrainingDataSet data;

  // Load training data for positive samples
  for (size_t sid = 0; sid < _param.positive_features.size(); ++sid) {
    GetPositiveDataForTumorNormalWgs(sid, _config.scoring_cols, data);
  }

  // Validate the number of positive training data points loaded
  const auto num_positive_data_points = data.GetDataPointCount();
  if (num_positive_data_points == 0) {
    throw error::Error(
        "No positive training data points were loaded. Please check your input files and configuration.");
  }

  // Load training data for negative samples
  for (size_t sid = 0; sid < _param.negative_features.size(); ++sid) {
    // If negative label needs downsampling, then we can add it to downsampling labels before loading sample 1.
    // Example:
    //   if (sid == 1) {
    //     data.downsampling_labels.insert(kDefaultNegativeLabel);
    //   }
    // No need to update downsampling labels for sample 0 because it is the first sample to be loaded.
    // No need to update downsampling labels for samples after sample 1 because there is only one negative label.

    GetNegativeDataForTumorNormalWgs(sid, _config.scoring_cols, data);
  }

  // Validate the number of negative training data points loaded
  const auto num_negative_data_points = data.GetDataPointCount() - num_positive_data_points;
  if (num_negative_data_points == 0) {
    throw error::Error(
        "No negative training data points were loaded. Please check your input files and configuration.");
  }

  ModelTrainer trainer(data, _config, _workflow);
  trainer.Train(_param.output_file, _param.threads, _param.iterations, VariantGroup::kAll, _param.output_training_data);
}

void TrainingController::TrainTumorOnlyTe() {
  VerifyOutputFilePath(_param.output_file);

  TrainingDataSet data;

  // Load training data for postive samples
  for (size_t sid = 0; sid < _param.positive_features.size(); ++sid) {
    GetPositiveDataForTumorOnlyTe(sid, _config.scoring_cols, data);
  }
  const auto num_positive_data_points = data.GetDataPointCount();
  if (num_positive_data_points == 0) {
    throw error::Error(
        "No positive training data points were loaded. Please check your input files and configuration.");
  }

  // Load training data for negative samples
  for (size_t sid = 0; sid < _param.negative_features.size(); ++sid) {
    GetNegativeDataForTumorOnlyTe(sid, _config.scoring_cols, data);
  }

  // Validate the number of negative training data points loaded
  const auto num_negative_data_points = data.GetDataPointCount() - num_positive_data_points;
  if (num_negative_data_points == 0) {
    throw error::Error(
        "No negative training data points were loaded. Please check your input files and configuration.");
  }

  ModelTrainer trainer(data, _config, _workflow);
  trainer.Train(_param.output_file, _param.threads, _param.iterations, VariantGroup::kAll, _param.output_training_data);
}

void TrainingController::TrainGermline() {
  // Check that we have valid model paths specified
  if (_param.indel_output_file.empty()) {
    throw error::Error("Please specify output indel model file path");
  }
  if (_param.snv_output_file.empty()) {
    throw error::Error("Please specify output SNV model file path");
  }
  CreateParentDirectoryIfNotExists(_param.indel_output_file);
  CreateParentDirectoryIfNotExists(_param.snv_output_file);

  TrainingDataSet snv_data;
  TrainingDataSet indel_data;

  // Load positive training data for the remaining samples
  for (size_t sid = 0; sid < _param.positive_features.size(); ++sid) {
    if (sid > 0) {
      // (Re)set the downsampling labels based on the distribution of labels currently in the training data.
      // This is repeated for every sample after the first sample is loaded, to account for potential differences in
      // the distribution of labels between samples.
      snv_data.DetermineDownsamplingLabels();
      indel_data.DetermineDownsamplingLabels();
    }
    GetPositiveDataForGermline(sid, _config.snv_scoring_cols, _config.indel_scoring_cols, snv_data, indel_data);
  }

  // Validate the number of training data points loaded for SNVs and indels
  const auto num_snv_data_points = snv_data.GetDataPointCount();
  const auto num_indel_data_points = indel_data.GetDataPointCount();
  if (num_snv_data_points == 0) {
    throw error::Error("No SNV training data points were loaded. Please check your input files and configuration.");
  }
  if (num_indel_data_points == 0) {
    throw error::Error("No indel training data points were loaded. Please check your input files and configuration.");
  }

  // Train model for SNVs
  {
    ModelTrainer snv_trainer(snv_data, _config, _workflow);
    snv_trainer.Train(_param.snv_output_file,
                      _param.threads,
                      _param.snv_iterations,
                      VariantGroup::kSnvOnly,
                      _param.snv_output_training_data);
    // Clear SNV training data from memory before training indel model, to reduce overall memory usage since SNV
    // training data is not needed for training indel model.
    snv_data.ClearData();
  }

  // Train model for indels
  ModelTrainer indel_trainer(indel_data, _config, _workflow);
  indel_trainer.Train(_param.indel_output_file,
                      _param.threads,
                      _param.indel_iterations,
                      VariantGroup::kIndelOnly,
                      _param.indel_output_training_data);
}

void TrainingController::TrainGermlineTagging() {
  VerifyOutputFilePath(_param.output_file);

  TrainingDataSet data;

  // Load positive training data for the remaining samples
  for (size_t sid = 0; sid < _param.positive_features.size(); ++sid) {
    GetPositiveDataForGermlineTagging(sid, _config.scoring_cols, data);
  }

  // Validate the number of training data points loaded
  const auto num_positive_data_points = data.GetDataPointCount();
  if (num_positive_data_points == 0) {
    throw error::Error(
        "No positive training data points were loaded. Please check your input files and configuration.");
  }

  ModelTrainer trainer(data, _config, _workflow);
  trainer.Train(_param.output_file, _param.threads, _param.iterations, VariantGroup::kAll, _param.output_training_data);
}

/**
 * @brief Train model using the specified parameters
 * @param param Model training parameters and input/output file paths.
 */
void TrainModel(const TrainModelParam& param) {
  using enum Workflow;
  TrainingController controller(param, param.config);
  switch (param.workflow) {
    case kGermlineMultiSample:
    case kGermline:
      controller.TrainGermline();
      break;
    case kGermlineTagging:
      controller.TrainGermlineTagging();
      break;
    case kTumorOnlyTe:
      controller.TrainTumorOnlyTe();
      break;
    case kTumorNormalWgs:
      controller.TrainTumorNormalWgs();
      break;
    default:
      throw error::Error("Unsupported workflow");
  }
}

}  // namespace xoos::svc
