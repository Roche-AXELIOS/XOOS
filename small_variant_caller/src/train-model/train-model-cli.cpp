#include "train-model-cli.h"

#include <xoos/cli/cli.h>
#include <xoos/cli/enum-option-util.h>
#include <xoos/cli/thread-count-option-util.h>
#include <xoos/cli/validators/file-permission-validator.h>
#include <xoos/util/string-functions.h>

#include "core/cli-option-names.h"
#include "util/cli-util.h"
#include "util/region-util.h"

namespace xoos::svc {

using cli::AddOptionalEnumOption;
using cli::AddThreadCountOption;
using enum Workflow;

// "custom" workflow not supported in train_model
constexpr std::array kSupportedWorkflows = {
    kGermline, kGermlineMultiSample, kTumorOnlyTe, kTumorNormalWgs, kGermlineTagging};

// default values
constexpr u32 kDefaultMaxWeightedScore = 4;
constexpr auto* kDefaultOutputFile = "model.txt.gz";
constexpr auto* kDefaultSnvOutputFile = "snv_model.txt.gz";
constexpr auto* kDefaultIndelOutputFile = "indel_model.txt.gz";

/**
 * @brief Extract a list of file paths from a text file.
 * @param argument_name Name of the CLI argument for error messages.
 * @param list_file_name Text file containing a list of file paths, one per line.
 * @param files Vector to store the extracted file paths.
 */
static void ParseFileList(const std::string& argument_name, const fs::path& list_file_name, vec<fs::path>& files) {
  const auto results = cli::CheckFileListReadable(list_file_name);
  if (std::holds_alternative<std::string>(results)) {
    throw cli::ValidationError("Invalid file list for '{}'; {}", argument_name, std::get<std::string>(results));
  }
  files = std::get<vec<fs::path>>(results);
}

/**
 * @brief Helper function to define CLI options for the main application.
 * These options are either required or important in this submodule, or they are common to all other submodules in SVC.
 * @param app Main application pointer where CLI options will be defined.
 * @param params Shared pointer to CLI parameters to store parsed option values
 */
static void AddMainOptions(CLI::App* const app, const TrainModelParamPtr& params) {
  AddWarnAsErrorOption(app);

  AddThreadCountOption(app, cli_opt_name::kThreads, params->threads);

  app->add_option(cli_opt_name::kConfig, params->config_file, "Path to config JSON file");
}

/**
 * @brief Helper function to add core CLI options for subcommands. CLI options added here are common across all
 * subcommands.
 * @param sub Subcommand application pointer where CLI options are to be added
 * @param params Shared pointer to CLI parameters to store parsed option values
 */
static void AddCoreOptions(CLI::App* const sub, const TrainModelParamPtr& params) {
  auto* const opt = sub->add_option_function<fs::path>(
                           cli_opt_name::kPosBamFeatures,
                           [&params](const fs::path& value) {
                             ParseFileList(cli_opt_name::kPosBamFeatures, value, params->positive_features);
                           },
                           "List of features files for positive samples expected to contain known variants")
                        ->required();
  CheckNonEmptyFile(opt);
}

/**
 * @brief Add positive VCF features CLI option.
 * @param sub CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
static void AddPosVcfFeaturesOption(CLI::App* const sub, const TrainModelParamPtr& params) {
  auto* const opt = sub->add_option_function<fs::path>(
                           cli_opt_name::kPosVcfFeatures,
                           [&params](const fs::path& value) {
                             ParseFileList(cli_opt_name::kPosVcfFeatures, value, params->positive_vcf_features);
                             if (!params->positive_features.empty() &&
                                 params->positive_features.size() != params->positive_vcf_features.size()) {
                               throw CLI::ValidationError(fmt::format("{} and {} must contain the same number of files",
                                                                      cli_opt_name::kPosVcfFeatures,
                                                                      cli_opt_name::kPosBamFeatures));
                             }
                           },
                           "List of vcf features files for positive samples expected to contain known variants")
                        ->required();
  CheckNonEmptyFile(opt);
}

/**
 * @brief Add truth VCFs CLI option.
 * @param sub CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
static void AddTruthVcfsOption(CLI::App* const sub, const TrainModelParamPtr& params) {
  auto* const opt = sub->add_option_function<fs::path>(
                           cli_opt_name::kTruthVcfs,
                           [&params](const fs::path& value) {
                             ParseFileList(cli_opt_name::kTruthVcfs, value, params->truth_vcfs);
                             if (!params->positive_features.empty() &&
                                 params->positive_features.size() != params->truth_vcfs.size()) {
                               throw CLI::ValidationError(fmt::format("{} and {} must contain the same number of files",
                                                                      cli_opt_name::kTruthVcfs,
                                                                      cli_opt_name::kPosBamFeatures));
                             }
                           },
                           "List of truth set VCFs for positive samples containing known variants in those samples")
                        ->required();
  CheckNonEmptyFile(opt);
}

/**
 * @brief Add negative BAM features CLI option.
 * @param sub CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
static void AddNegBamFeaturesOption(CLI::App* const sub, const TrainModelParamPtr& params) {
  auto* const opt = sub->add_option_function<fs::path>(
      cli_opt_name::kNegBamFeatures,
      [&params](const fs::path& value) {
        ParseFileList(cli_opt_name::kNegBamFeatures, value, params->negative_features);
      },
      "List of features files for negative/normal/healthy samples");
  CheckNonEmptyFile(opt);
}

/**
 * @brief Add negative VCF features CLI option.
 * @param sub CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
static void AddNegVcfFeaturesOption(CLI::App* const sub, const TrainModelParamPtr& params) {
  auto* const opt = sub->add_option_function<fs::path>(
      cli_opt_name::kNegVcfFeatures,
      [&params](const fs::path& value) {
        ParseFileList(cli_opt_name::kNegVcfFeatures, value, params->negative_vcf_features);
        if (!params->negative_features.empty() &&
            params->negative_features.size() != params->negative_vcf_features.size()) {
          throw CLI::ValidationError(fmt::format("{} and {} must contain the same number of files",
                                                 cli_opt_name::kNegVcfFeatures,
                                                 cli_opt_name::kNegBamFeatures));
        }
      },
      "List of vcf features files for negative/normal/healthy samples");
  CheckNonEmptyFile(opt);
}

/**
 * @brief Add CLI options common to non-germline workflows.
 * @param sub CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 * @param defaults Workflow-specific default configuration
 */
static void AddNonGermlineSpecificOptions(CLI::App* const sub,
                                          const TrainModelParamPtr& params,
                                          const SVCConfig& defaults) {
  sub->add_option(cli_opt_name::kIterations,
                  params->iterations,
                  "The maximum number of training rounds the model will run before stopping (inclusive)")
      ->default_val(defaults.iterations)
      ->check(CLI::PositiveNumber);

  sub->add_option(cli_opt_name::kOutputFile, params->output_file, "Output path for model file")
      ->default_val(kDefaultOutputFile)
      ->check(CLI::NonexistentPath);

  sub->add_option(
         cli_opt_name::kOutputTrainingDataTsv, params->output_training_data, "Output path for model training data TSV")
      ->check(CLI::NonexistentPath);
}

/**
 * @brief Helper function to add germline and germline-multi-sample specific CLI options.
 * @param app CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
static void AddGermlineSpecificOptions(const CLI::App* const app, const TrainModelParamPtr& params) {
  for (const Workflow workflow : {kGermline, kGermlineMultiSample}) {
    CLI::App* const sub = app->get_subcommand(enum_util::FormatEnumName(workflow));
    const auto defaults = SVCConfig(workflow);
    AddPosVcfFeaturesOption(sub, params);
    AddTruthVcfsOption(sub, params);

    sub->add_option(cli_opt_name::kSnvIterations,
                    params->snv_iterations,
                    "The maximum number of training rounds the SNV model will run before stopping (inclusive)")
        ->default_val(defaults.snv_iterations)
        ->check(CLI::PositiveNumber);

    sub->add_option(cli_opt_name::kIndelIterations,
                    params->indel_iterations,
                    "The maximum number of training rounds the indel model will run before stopping (inclusive)")
        ->default_val(defaults.indel_iterations)
        ->check(CLI::PositiveNumber);

    sub->add_option(cli_opt_name::kSnvModelOutput, params->snv_output_file, "Output path for SNV model")
        ->default_val(kDefaultSnvOutputFile)
        ->check(CLI::NonexistentPath);

    sub->add_option(cli_opt_name::kIndelModelOutput, params->indel_output_file, "Output path for indel model")
        ->default_val(kDefaultIndelOutputFile)
        ->check(CLI::NonexistentPath);

    sub->add_option(cli_opt_name::kSnvOutputTrainingDataTsv,
                    params->snv_output_training_data,
                    "Output path for SNV model training data TSV")
        ->check(CLI::NonexistentPath);

    sub->add_option(cli_opt_name::kIndelOutputTrainingDataTsv,
                    params->indel_output_training_data,
                    "Output path for indel model training data TSV")
        ->check(CLI::NonexistentPath);
  }
}

/**
 * @brief Helper function to add tumor-only-te specific CLI options.
 * @param app CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
static void AddTumorOnlyTeSpecificOptions(const CLI::App* const app, const TrainModelParamPtr& params) {
  CLI::App* const sub = app->get_subcommand(enum_util::FormatEnumName(kTumorOnlyTe));
  const auto defaults = SVCConfig(kTumorOnlyTe);
  AddNegBamFeaturesOption(sub, params);
  AddNonGermlineSpecificOptions(sub, params, defaults);

  sub->add_option(cli_opt_name::kMaxWeightedScore,
                  params->max_score,
                  "The maximum weighted score to use for a feature to be considered negative (inclusive)")
      ->default_val(kDefaultMaxWeightedScore);

  auto* const blocklist_opt = sub->add_option_function<fs::path>(
      cli_opt_name::kBlocklistBed,
      [&params](const fs::path& value) { params->blocklist = GetChromIntervalMap(value); },
      "BED file of 0-based regions to exclude from model training");
  CheckBedFile(blocklist_opt);
}

/**
 * @brief Helper function to add tumor-normal-wgs specific CLI options.
 * @param app CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
static void AddTumorNormalWgsSpecificOptions(const CLI::App* const app, const TrainModelParamPtr& params) {
  CLI::App* const sub = app->get_subcommand(enum_util::FormatEnumName(kTumorNormalWgs));
  const auto defaults = SVCConfig(kTumorNormalWgs);
  AddPosVcfFeaturesOption(sub, params);
  AddNegBamFeaturesOption(sub, params);
  AddNegVcfFeaturesOption(sub, params);
  AddTruthVcfsOption(sub, params);
  AddNonGermlineSpecificOptions(sub, params, defaults);
}

/**
 * @brief Helper function to add germline-tagging specific CLI options.
 * @param app CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
static void AddGermlineTaggingSpecificOptions(const CLI::App* const app, const TrainModelParamPtr& params) {
  CLI::App* const sub = app->get_subcommand(enum_util::FormatEnumName(kGermlineTagging));
  const auto defaults = SVCConfig(kGermlineTagging);
  AddPosVcfFeaturesOption(sub, params);
  AddTruthVcfsOption(sub, params);
  AddNonGermlineSpecificOptions(sub, params, defaults);
}

void train_model::DefineOptions(CLI::App* const app, TrainModelParamPtr& params) {
  AddMainOptions(app, params);

  // Add a subcommand for each workflow
  for (const Workflow workflow : kSupportedWorkflows) {
    const std::string name = enum_util::FormatEnumName(workflow);
    const std::string desc = fmt::format("Train model for the {} workflow", name);
    CLI::App* const sub = app->add_subcommand(name, desc)->fallthrough();
    // Do not apply force_callback() to subcommand options to avoid overwriting params set by other subcommands
    AddCoreOptions(sub, params);
    const auto defaults = SVCConfig(workflow);
    cli::AddEnumOption(sub,
                       cli_opt_name::kNormalizeFeatures,
                       params->normalize_features,
                       "Normalize read depth related features values",
                       defaults.normalize_features);
  }
  app->require_subcommand(kMinSubcommands, kMaxSubcommands);

  // Add workflow-specific CLI options for each subcommand
  AddGermlineSpecificOptions(app, params);
  AddTumorOnlyTeSpecificOptions(app, params);
  AddTumorNormalWgsSpecificOptions(app, params);
  AddGermlineTaggingSpecificOptions(app, params);

  // hide the `germline-tagging` subcommand by assigning it to an empty string group
  app->get_subcommand(enum_util::FormatEnumName(kGermlineTagging))->group("");
}

/**
 * @brief Helper function to set default output model paths if not provided via CLI.
 * @param sub Subcommand CLI application pointer to check which options were specified via CLI
 * @param params CLI parameters shared pointer to apply defaults to
 */
static void SetDefaultOutputModelPaths(const CLI::App* const sub, const TrainModelParamPtr& params) {
  if (params->workflow == kGermline || params->workflow == kGermlineMultiSample) {
    if (sub->count(cli_opt_name::kSnvModelOutput) == 0) {
      params->snv_output_file = kDefaultSnvOutputFile;
    }
    if (sub->count(cli_opt_name::kIndelModelOutput) == 0) {
      params->indel_output_file = kDefaultIndelOutputFile;
    }
  } else {
    if (sub->count(cli_opt_name::kOutputFile) == 0) {
      params->output_file = kDefaultOutputFile;
    }
  }
}

/**
 * @brief Helper function to apply CLI defaults and workflow config presets to CLI parameters if the corresponding CLI
 * options were not specified.
 * @param sub Subcommand CLI application pointer to check which options were specified via CLI
 * @param params CLI parameters shared pointer to apply defaults to
 */
static void ApplyConfig(const CLI::App* const sub, const TrainModelParamPtr& params) {
  SetDefaultOutputModelPaths(sub, params);

  if (sub->count(cli_opt_name::kNormalizeFeatures) == 0) {
    params->normalize_features = params->config.normalize_features;
  }
  if (sub->get_option_no_throw(cli_opt_name::kIterations) == nullptr || sub->count(cli_opt_name::kIterations) == 0) {
    params->iterations = params->config.iterations;
  }
  if (sub->get_option_no_throw(cli_opt_name::kSnvIterations) == nullptr ||
      sub->count(cli_opt_name::kSnvIterations) == 0) {
    params->snv_iterations = params->config.snv_iterations;
  }
  if (sub->get_option_no_throw(cli_opt_name::kIndelIterations) == nullptr ||
      sub->count(cli_opt_name::kIndelIterations) == 0) {
    params->indel_iterations = params->config.indel_iterations;
  }
}

void train_model::PreCallback(const cli::ConstAppPtr app, const TrainModelParamPtr& params) {
  params->command_line = GetCommandLineInfo(app);

  // Check which subcommand was used, set workflow and config accordingly, and apply config defaults as needed
  for (const Workflow workflow : kSupportedWorkflows) {
    const std::string name = enum_util::FormatEnumName(workflow);
    if (app->got_subcommand(name)) {
      params->workflow = workflow;
      params->config = JsonToConfig(params->config_file.value_or(fs::path{}), params->workflow);
      ApplyConfig(app->get_subcommand(name), params);
      break;
    }
  }
}

}  // namespace xoos::svc
