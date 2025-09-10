#include "train-model-cli.h"

#ifdef SOMATIC_ENABLE
#include <fstream>
#include <vector>
#endif  // SOMATIC_ENABLE

#include <xoos/cli/cli.h>
#include <xoos/cli/enum-option-util.h>
#include <xoos/cli/thread-count-option-util.h>
#include <xoos/cli/validators/file-permission-validator.h>
#include <xoos/util/string-functions.h>

#include "util/cli-util.h"

namespace xoos::svc {

using BedRegion = io::BedRegion;

using cli::AddOptionalEnumOption;

/**
 * @brief Extract a list of file paths from a text file.
 * @param argument_name Name of the CLI argument for error messages.
 * @param list_file_name Text file containing a list of file paths, one per line.
 * @param feature_files Vector to store the extracted file paths.
 */
static void ParseFileList(const std::string& argument_name, const fs::path& list_file_name, vec<fs::path>& files) {
  const auto results = cli::CheckFileListReadable(list_file_name);
  if (std::holds_alternative<std::string>(results)) {
    throw cli::ValidationError("Invalid file list for '{}'; {}", argument_name, std::get<std::string>(results));
  }
  files = std::get<vec<fs::path>>(results);
}

#ifdef SOMATIC_ENABLED
/**
 * @brief Extract target regions from a BED file.
 * @param file_name Path to the BED file.
 * @param bed_regions Vector to store the extracted target regions.
 */
static void ProcessRegionsFile(const fs::path& file_name, std::vector<BedRegion>& bed_regions) {
  bed_regions.clear();
  std::string line;
  std::ifstream bed(file_name);
  while (getline(bed, line)) {
    string::Trim(line);
    bed_regions.push_back(io::ParseRegion(line));
  }
}
#endif  // SOMATIC_ENABLED

/**
 * @brief Pre-callback function for model training CLI.
 */
void TrainModelPreCallback(cli::ConstAppPtr app, const TrainModelParamPtr& params) {
  params->command_line = GetCommandLineInfo(app);
  using enum Workflow;
  switch (params->workflow) {
#ifdef SOMATIC_ENABLE
    case kSomatic: {
      break;
    }
#endif  // SOMATIC_ENABLE
    case kGermlineMultiSample:
    case kGermline: {
      if (app->count("--positive-vcf-features") == 0) {
        throw CLI::ValidationError(
            "positive-vcf-features is required for the germline or germline-multi-sample workflows");
      }
      break;
    }
    default: {
      throw CLI::ValidationError("Workflow not supported");
    }
  }
}

/**
 * @brief Load the JSON configuration file and update model training parameters.
 * @param params Model training parameters that will be updated.
 * @param config_json Path to the JSON configuration file.
 * @return SVCConfig object with loaded configuration.
 */
static SVCConfig GetConfig(const TrainModelParamPtr& params, const fs::path& config_json) {
  SVCConfig model_config = JsonToConfig(config_json, params->workflow);
  // Load the default values from the config. If the user overwrites a value it is done below
  params->iterations = model_config.iterations;
  params->snv_iterations = model_config.snv_iterations;
  params->indel_iterations = model_config.indel_iterations;
  params->normalize_features = model_config.normalize_features;
  return model_config;
}

void DefineOptionsTrainModel(cli::AppPtr app, const TrainModelParamPtr& params) {
  app->add_option_function<fs::path>(
         "--positive-features",
         [&params](const fs::path& value) { ParseFileList("--positive-features", value, params->positive_features); },
         "List of features files for positive samples expected to contain known variants")
      ->required()
      ->check(kCliNonEmptyFile);
  app->add_option_function<fs::path>(
         "--positive-vcf-features",
         [&params](const fs::path& value) {
           ParseFileList("--positive-vcf-features", value, params->positive_vcf_features);
           if (!params->positive_features.empty() &&
               params->positive_features.size() != params->positive_vcf_features.size()) {
             throw CLI::ValidationError(
                 "positive-vcf-features must contain the same number of VCFs as positive-features");
           }
         },
         "List of vcf features files for positive samples expected to contain known variants")
      ->check(kCliNonEmptyFile);
#ifdef SOMATIC_ENABLED
  app->add_option_function<fs::path>(
         "--negative-features",
         [&params](const fs::path& value) { ParseFileList("--negative-features", value, params->negative_features); },
         "List of features files for negative/normal/healthy samples")
      ->check(kCliNonEmptyFile);
  app->add_option_function<fs::path>(
         "--negative-vcf-features",
         [&params](const fs::path& value) {
           ParseFileList("--negative-vcf-features", value, params->negative_vcf_features);
           if (!params->negative_features.empty() &&
               params->negative_features.size() != params->negative_vcf_features.size()) {
             throw CLI::ValidationError(
                 "negative-vcf-features must contain the same number of vcfs as negative-features");
           }
         },
         "List of vcf features files for negative/normal/healthy samples")
      ->check(kCliNonEmptyFile);
  app->add_option_function<fs::path>(
         "--blocklist",
         [&params](const fs::path& value) { ProcessRegionsFile(value, params->blocklist); },
         "0-based blocklist BED file of regions/positions to exclude from model training")
      ->check(kCliNonEmptyFile);
  app->add_option("--max-score",
                  params->max_score,
                  "The maximum weighted score to use for a feature to be considered negative (inclusive)")
      ->default_val(kMaxScore);
  app->add_option("--train-iterations",
                  params->iterations,
                  "The maximum number of training rounds the model will run before stopping (inclusive)")
      ->check(CLI::PositiveNumber);
#endif  // SOMATIC_ENABLED
  app->add_option_function<fs::path>(
         "--truth-vcfs",
         [&params](const fs::path& value) {
           ParseFileList("--truth-vcfs", value, params->truth_vcfs);
           if (!params->positive_features.empty() && params->positive_features.size() != params->truth_vcfs.size()) {
             throw CLI::ValidationError("truth-vcfs must contain the same number of VCFs as positive-features");
           }
         },
         "List of truth set VCFs for positive samples containing known variants in those samples")
      ->required()
      ->check(kCliNonEmptyFile);
  app->add_option("--output-file", params->output_file, "Output model file(s)");
  AddOptionalEnumOption(app, "--workflow", params->workflow, "Compute features for the designated workflow")
      ->required();
  // The `--workflow` option must be defined before the `--config` option.
  // Otherwise, `params->workflow` will always have the default value (somatic) within the `GetConfig` function.
  app->add_option_function<fs::path>(
         "--config",
         [&params](const fs::path& value) { params->config = GetConfig(params, value); },
         "Path to config JSON file")
      ->force_callback();
  app->add_option("--train-snv-iterations",
                  params->snv_iterations,
                  "The maximum number of training rounds the germline SNV model will run before stopping (inclusive)")
      ->check(CLI::PositiveNumber);
  app->add_option("--train-indel-iterations",
                  params->indel_iterations,
                  "The maximum number of training rounds the germline indel model will run before stopping (inclusive)")
      ->check(CLI::PositiveNumber);
  cli::AddThreadCountOption(app, "--threads", params->threads);
  app->add_flag("--normalize,!--no-normalize", params->normalize_features, "Normalize read count features.");
  app->add_flag("--write-training-data-tsv,!--no-write-training-data-tsv",
                params->write_training_data_tsv,
                "Write training data for SNV and indel models to TSV file.");
  AddWarnAsErrorOption(app);
}

}  // namespace xoos::svc
