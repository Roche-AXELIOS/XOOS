#include "model-trainer.h"

#include <zlib.h>
#include <zstd.h>

#include <optional>

#include <LightGBM/c_api.h>
#include <fmt/format.h>

#include <csv.hpp>

#include <xoos/error/error.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/log/logging.h>
#include <xoos/util/string-functions.h>

#include "core/variant-feature-extraction.h"
#include "util/compress/compress.h"
#include "util/lightgbm-util.h"
#include "util/locked-tsv-writer.h"
#include "util/log-util.h"
#include "util/seq-util.h"

namespace xoos::svc {

/**
 * @brief Extract the genotypes for ground truth variants.
 * @param vcf VCF file containing ground truth variants and their genotypes.
 * @return Ground truth variant genotypes
 */
TruthFeatures GetTruthFeatures(const fs::path& vcf) {
  TruthFeatures truth_variants;
  const io::VcfReader reader(vcf);
  while (const auto& vcf_record = reader.GetNextRecord()) {
    const auto& ref = vcf_record->Allele(0);
    if (ContainsOnlyACTG(ref)) {
      auto& [genotype, ref_alt_set] = truth_variants[vcf_record->Chromosome()][vcf_record->Position()];
      // Unsupported genotypes in the ground truth are still kept because variants with these genotypes should be
      // excluded from the training data.
      genotype = StringToGenotype(vcf_record->GetGTField());
      // Multi-allelics are rare compared to heterozygous REF-ALT and homozygous ALT.
      // To include as many true alleles as possible in the training data, store the ground truth GT label for each
      // individual ALT allele separately.
      // This is intended to handle two scenarios:
      // 1. Multi-allelics with 1 true allele and 1 false allele.
      // 2. Multi-allelics with 1 SNV and 1 indel. The germline workflow has separate training data for SNV and indel.
      const int num_alleles = vcf_record->NumAlleles();
      for (int i = 1; i < num_alleles; ++i) {
        const auto& [ref_trimmed, alt_trimmed] = TrimVariant(ref, vcf_record->Allele(i));
        if (ContainsOnlyACTG(alt_trimmed)) {
          ref_alt_set.emplace(ref_trimmed, alt_trimmed);
        }
      }
    }
  }
  return truth_variants;
}

ModelTrainer::ModelTrainer(ChromToVariantInfoMapWithLabel& features,
                           RefInfoMapMultiSample& ref_features,
                           ChromToVcfFeaturesMapMultiSample& vcf_features,
                           SVCConfig& config,
                           Workflow workflow,
                           vec<StrUnorderedMap<u32>>& normalize_targets)
    : _features(features),
      _ref_features(ref_features),
      _vcf_features(vcf_features),
      _config(config),
      _workflow(workflow),
      _normalize_targets(normalize_targets) {
}

/**
 * @brief Get the workflow-specific LightGBM configuration string from the model training configuration.
 * @param workflow Workflow type
 * @param cat_names List of categorical feature names
 * @param config Model training configuration and parameters
 * @param var_type Variant type
 * @param num_threads Number of threads for model training
 * @return LightGBM configuration string
 */
static std::string GetLightGBMConfig(const Workflow workflow,
                                     const std::vector<std::string>& cat_names,
                                     const SVCConfig& config,
                                     const VariantGroup var_type,
                                     const size_t num_threads) {
  std::string lightgbm_config_str{fmt::format("num_threads={} force_row_wise=true deterministic=true ", num_threads)};
  switch (workflow) {
    case Workflow::kGermlineMultiSample:
    case Workflow::kGermline:
      // Use Different settings for indel vs SNV models
      if (var_type == VariantGroup::kIndelOnly) {
        lightgbm_config_str += config.indel_model_lgbm_params;
      } else {
        lightgbm_config_str += config.snv_model_lgbm_params;
      }
      break;
    default:
      throw error::Error("Specified Workflow Not Implemented");
  }
  if (!cat_names.empty()) {
    lightgbm_config_str += " categorical_feature=name:" + string::Join(cat_names, ",");
  }
  return lightgbm_config_str;
}

/**
 * @brief Helper function to add rows to output data vector for writing to file
 * @param data_rows Vector or formatted rows for output
 * @param feature_vec Vector of features for the row
 * @param label Genotype label for the row
 * @param vid Variant ID
 * @param is_germline Flag whether it is germline workflow
 */
static void WriteDataRowHelper(vec<vec<std::string>>& data_rows,
                               const vec<double>& feature_vec,
                               const std::string& label,
                               const VariantId& vid,
                               const bool is_germline) {
  // write data row to file
  vec<std::string> data_row = {vid.chrom, std::to_string(vid.pos), vid.ref, vid.alt};
  if (is_germline) {
    const auto genotype = IntToGenotype(std::stoul(label));
    data_row.emplace_back(GenotypeToString(genotype));
  } else {
    data_row.emplace_back(label);
  }
  for (const auto val : feature_vec) {
    data_row.emplace_back(std::to_string(val));
  }
  data_rows.emplace_back(data_row);
}

/**
 * @brief Helper function to write rows of data to file
 * @param data_file Filename of output file to append to
 * @param data_rows Vector or formatted rows for output
 */
static void FlushDataHelper(LockedTsvWriter& writer, const vec<vec<std::string>>& data_rows) {
  writer.AppendRows(data_rows);
  writer.Flush();
}

/**
 * @brief Helper function to serialize the LightGBM booster to a string.
 * @param booster LightGBM booster
 * @return Model content as a string
 */
static std::string SaveModelToString(BoosterHandle booster) {
  return lightgbm::BoosterSaveModelToString(booster, 0, -1, C_API_FEATURE_IMPORTANCE_SPLIT);
}

/**
 * @brief Train a LightGBM model based on the model training configurations.
 * @param output_file The model file to save
 * @param num_threads Number of threads to use
 * @param num_rounds The maximum number of rounds of training.
 * @param var_group Type of variant to train
 * @param reduce_redundancy Reduce redundancy in training data
 * @param write_data Write training data to file
 */
void ModelTrainer::Train(const fs::path& output_file,
                         const size_t num_threads,
                         const u32 num_rounds,
                         const VariantGroup var_group,
                         bool reduce_redundancy,
                         const bool write_data) {
  // 1. Set up LightGBM scoring column names and categorical names based on the variant type.
  // 2. Set up a writer for the training data file.
  // 3. Determine over-represented truth labels in the training data if `reduce_redundancy` is true.
  // 4. Create LightGBM training dataset.
  // 5. Create LightGBM booster.
  // 6. Train the model for the specified number of rounds.
  // 7. Save the model to the output file.

  const bool is_germline =
      _config.workflow == Workflow::kGermline || _config.workflow == Workflow::kGermlineMultiSample;
  // set up scoring names, scoring columns, and categorical names
  vec<std::string> scoring_names;
  vec<UnifiedFeatureCols> scoring_cols;
  vec<std::string> cat_names;
  switch (var_group) {
    case VariantGroup::kSnvOnly:
      scoring_names = _config.snv_scoring_names;
      scoring_cols = _config.snv_scoring_cols;
      cat_names = _config.snv_categorical_names;
      break;
    case VariantGroup::kIndelOnly:
      scoring_names = _config.indel_scoring_names;
      scoring_cols = _config.indel_scoring_cols;
      cat_names = _config.indel_categorical_names;
      break;
    default:
      scoring_names = _config.scoring_names;
      scoring_cols = _config.scoring_cols;
      cat_names = _config.categorical_names;
      break;
  }

  // set up lightGBM column names
  vec<const char*> columns;
  columns.reserve(scoring_names.size());
  for (const auto& col : scoring_names) {
    columns.push_back(col.c_str());
  }
  // set up writer for training data file
  auto writer = write_data ? std::make_optional<LockedTsvWriter>(output_file.string() + ".data.tsv") : std::nullopt;
  vec<vec<std::string>> data_rows;
  if (write_data) {
    vec<std::string> header{"CHROM", "POS", "REF", "ALT", "LABEL"};
    header.insert(header.end(), scoring_names.begin(), scoring_names.end());
    writer->AppendRow(header);
  }
  // We can reduce redundancy in the training data by removing samples for variants with the same truth label.
  // Step 1:
  //   Determine over-represented truth labels in the training data.
  // Step 2:
  //   For a given variant with the same over-represented label in multiple samples, arbitrarily (i.e. not randomly)
  //   select one sample for the label. This ensures that the training data is reproducible for the same command line.
  StrUnorderedSet dominant_labels{};
  StrMap<u64> label_counts;
  if (reduce_redundancy) {
    // count each label in the training data and determine over-represented truth labels
    bool is_multi_sample = false;
    for (const auto& [chrom, pos_data] : _features) {
      for (const auto& [pos, vid_data] : pos_data) {
        for (const auto& [vid, feat_label_sid_vec] : vid_data) {
          if (var_group != VariantGroup::kAll) {
            const bool is_snv = vid.type == VariantType::kSNV;
            if ((var_group == VariantGroup::kSnvOnly && !is_snv) || (var_group == VariantGroup::kIndelOnly && is_snv)) {
              continue;
            }
          }
          for (const auto& [var_feat, label, sid] : feat_label_sid_vec) {
            ++label_counts[label];
          }
          if (!is_multi_sample && feat_label_sid_vec.size() > 1) {
            // at least one variant has more than one sample for one label
            is_multi_sample = true;
          }
        }
      }
    }
    if (is_multi_sample) {
      // find labels with higher-than-expected counts in the training data
      u64 total_count = 0;
      for (const auto& [label, count] : label_counts) {
        total_count += count;
      }
      const u64 balanced_count = total_count / label_counts.size();
      for (const auto& [label, count] : label_counts) {
        if (count >= balanced_count) {
          dominant_labels.emplace(label);
          if (is_germline) {
            // `label` is an integer in string format, e.g. "0"
            // convert to genotype in string format for logging, e.g. "0/0"
            auto genotype = GenotypeToString(IntToGenotype(std::stoul(label)));
          }
        }
      }
    } else {
      // turn off sampling of training data
      reduce_redundancy = false;
    }
    // reset the counts
    label_counts = {};
  }
  const bool has_dominant_labels = !dominant_labels.empty();
  // set up LightGBM training data
  vec<double> training_data;
  vec<float> labels;
  for (const auto& [chrom, pos_data] : _features) {
    for (const auto& [pos, vid_data] : pos_data) {
      for (const auto& [vid, feat_label_sid_vec] : vid_data) {
        if (var_group != VariantGroup::kAll) {
          const bool is_snv = vid.type == VariantType::kSNV;
          if ((var_group == VariantGroup::kSnvOnly && !is_snv) || (var_group == VariantGroup::kIndelOnly && is_snv)) {
            continue;
          }
        }
        StrUnorderedMap<u32> dominant_label_to_sid{};
        bool check_label = reduce_redundancy && has_dominant_labels && feat_label_sid_vec.size() > 1;
        if (check_label) {
          // extract sample IDs for the dominant truth labels only
          StrUnorderedMap<vec<u32>> dominant_label_to_sids;
          for (const auto& [feat, label, sid] : feat_label_sid_vec) {
            if (dominant_labels.contains(label)) {
              dominant_label_to_sids[label].emplace_back(sid);
            }
          }
          // check label afterward only if at least one dominant label has multi-sample data
          check_label = false;
          for (const auto& [label, sids] : dominant_label_to_sids) {
            if (sids.size() == 1) {
              dominant_label_to_sid[label] = sids.at(0);
            } else {
              // Select one arbitrary sample ID using the modulo operator to ensure reproducibility of training data.
              // Variant position is usually much larger than the number of samples, and it is sufficiently
              // unique for all variants in the training data. So, we can use it as a hash value in this context.
              // Contrary to random sampling, this helps reduce stochasticity in multi-sample models.
              dominant_label_to_sid[label] = sids.at(vid.pos % sids.size());
              check_label = true;
            }
          }
        }
        for (const auto& [var_feat, label, sid] : feat_label_sid_vec) {
          if (check_label && dominant_labels.contains(label) && dominant_label_to_sid.at(label) != sid) {
            // Skip training data from this sample
            continue;
          }
          UnifiedReferenceFeature ref_feature;
          try {
            ref_feature = _ref_features.at(chrom).at(sid).at(vid.GetRefFeaturePos());
          } catch ([[maybe_unused]] std::out_of_range& e) {
          }
          VcfFeature vcf_feat;
          bool vcf_feat_found = false;
          if (_vcf_features.contains(chrom) && _vcf_features.at(chrom).contains(pos) &&
              _vcf_features.at(chrom).at(pos).contains(vid)) {
            auto vcf_data = _vcf_features.at(chrom).at(pos).at(vid);
            for (const auto& [vcf_info, vcf_sample] : vcf_data) {
              if (vcf_sample == sid) {
                vcf_feat = vcf_info;
                vcf_feat_found = true;
                break;
              }
            }
          }
          if (!is_germline || vcf_feat_found) {
            // germline models require both BAM and VCF features
            std::optional<u32> normalize_target;
            if (!_normalize_targets.empty()) {
              try {
                normalize_target = _normalize_targets.at(sid).at(chrom);
              } catch ([[maybe_unused]] std::out_of_range& e) {
                // no normalize target for this chromosome
              }
            }
            const auto& feature_vec =
                GetFeatureVec(scoring_cols, vid, var_feat, ref_feature, vcf_feat, normalize_target);
            training_data.insert(training_data.end(), feature_vec.begin(), feature_vec.end());
            labels.emplace_back(std::stof(label));
            label_counts[label]++;
            if (writer) {
              // Format data and store in data_rows buffer.
              WriteDataRowHelper(data_rows, feature_vec, label, vid, is_germline);
              // Output when the data_rows buffer has reached a certain size.
              if (data_rows.size() > 10000) {
                FlushDataHelper(*writer, data_rows);
                data_rows.clear();
              }
            }
          }
        }
      }
    }
  }
  if (writer) {
    // Flush any remaining rows that may be in data_rows buffer
    FlushDataHelper(*writer, data_rows);
  }
  {
    // print the label counts of training data
    vec<std::string> items{};
    if (is_germline) {
      for (const auto& [label, count] : label_counts) {
        auto genotype = IntToGenotype(std::stoul(label));
        items.emplace_back(fmt::format("{}:{}", GenotypeToString(genotype), count));
      }
    } else {
      for (const auto& [label, count] : label_counts) {
        items.emplace_back(fmt::format("{}:{}", label, count));
      }
    }
    Logging::Info("Label counts: {}", string::Join(items, ", "));
    if (label_counts.size() != _config.n_classes) {
      WarnAsErrorIfSet(
          "Number of unique labels ({}) and classes ({}) are different!", label_counts.size(), _config.n_classes);
    }
  }
  auto lightgbm_config_str = GetLightGBMConfig(_workflow, cat_names, _config, var_group, num_threads);
  Logging::Info("Training with config: {}", lightgbm_config_str);

  auto num_rows = static_cast<s32>(labels.size());
  auto num_cols = static_cast<s32>(columns.size());
  DatasetHandle train_data = nullptr;
  lightgbm::DatasetCreateFromMat(
      training_data.data(), C_API_DTYPE_FLOAT64, num_rows, num_cols, 1, "", nullptr, &train_data);
  _train_data.reset(train_data);

  lightgbm::DatasetSetFeatureNames(_train_data.get(), columns.data(), num_cols);
  lightgbm::DatasetSetField(_train_data.get(), "label", labels.data(), num_rows, C_API_DTYPE_FLOAT32);

  BoosterHandle booster = nullptr;
  lightgbm::BoosterCreate(_train_data.get(), lightgbm_config_str, &booster);
  _booster.reset(booster);

  for (u32 train_round = 0; train_round < num_rounds; ++train_round) {
    if (0 == (train_round % 100)) {
      Logging::Info("Training round {}", train_round);
    }
    if (lightgbm::BoosterUpdateOneIter(_booster.get())) {
      Logging::Info("Training finished at round {}", train_round);
      break;
    }
  }
  Logging::Info("Saving model to {}", output_file);
  if (IsCompressed(output_file)) {
    Compress(SaveModelToString(_booster.get()), output_file);
  } else {
    lightgbm::BoosterSaveModel(_booster.get(), 0, -1, C_API_FEATURE_IMPORTANCE_SPLIT, output_file);
  }
}

}  // namespace xoos::svc
