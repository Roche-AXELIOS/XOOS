#include "train-model-cli.h"

#include <string>
#include <string_view>

#include <xoos/cli/cli.h>
#include <xoos/cli/enum-option-util.h>
#include <xoos/cli/thread-count-option-util.h>
#include <xoos/cli/validators/bed-validator.h>
#include <xoos/cli/validators/file-permission-validator.h>
#include <xoos/cli/validators/nonempty-file-validator.h>
#include <xoos/util/string-functions.h>

#include "core/cli-option-names.h"
#include "util/cli-util.h"
#include "util/region-util.h"

namespace xoos::svc {

using cli::AddThreadCountOption;
using enum Workflow;

// "custom" workflow not supported in train_model
constexpr std::array kSupportedWorkflows = {
    kGermlineWgs, kTumorOnlyTe, kTumorOnlyWgs, kTumorNormalWgs, kGermlineTagging};

// default values
constexpr u32 kDefaultMaxWeightedScore = 4;
constexpr auto* kDefaultOutputFile = "model.txt.gz";
constexpr auto* kDefaultSnvOutputFile = "snv_model.txt.gz";
constexpr auto* kDefaultIndelOutputFile = "indel_model.txt.gz";

namespace {
void SyncGermlineTrainingRenderedDefaults(CLI::App* sub, const TrainModelParamPtr& params);

/**
 * @brief Extract a list of file paths from a text file.
 * @param argument_name Name of the CLI argument for error messages.
 * @param list_file_name Text file containing a list of file paths, one per line.
 * @param files Vector to store the extracted file paths.
 */
void ParseFileList(const std::string& argument_name, const fs::path& list_file_name, vec<fs::path>& files) {
  const auto results = cli::CheckFileListReadable(list_file_name);
  if (std::holds_alternative<std::string>(results)) {
    throw cli::ValidationError("Invalid file list for '{}'; {}", argument_name, std::get<std::string>(results));
  }
  files = std::get<vec<fs::path>>(results);
}

/**
 * @brief Helper function to add core CLI options for subcommands. CLI options added here are common across all
 * subcommands.
 * @param sub Subcommand application pointer where CLI options are to be added
 * @param params Shared pointer to CLI parameters to store parsed option values
 */
void AddCommonOptions(CLI::App* const sub, const TrainModelParamPtr& params) {
  AddWarnAsErrorOption(sub);

  AddThreadCountOption(sub, cli_opt_name::kThreads, params->threads);

  sub->add_option(cli_opt_name::kConfig, params->config_file, "Path to config JSON file");

  auto* const opt = sub->add_option_function<fs::path>(
                           cli_opt_name::kPosBamFeatures,
                           [&params](const fs::path& value) {
                             ParseFileList(cli_opt_name::kPosBamFeatures, value, params->positive_features);
                           },
                           "List of features files for positive samples expected to contain known variants")
                        ->required();
  opt->check(cli::NonEmptyFileValidator());
}

/**
 * @brief Add positive VCF features CLI option.
 * @param sub CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
void AddPosVcfFeaturesOption(CLI::App* const sub, const TrainModelParamPtr& params) {
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
  opt->check(cli::NonEmptyFileValidator());
}

/**
 * @brief Add truth VCFs CLI option.
 * @param sub CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
void AddTruthVcfsOption(CLI::App* const sub, const TrainModelParamPtr& params) {
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
  opt->check(cli::NonEmptyFileValidator());
}

/**
 * @brief Add negative BAM features CLI option.
 * @param sub CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
void AddNegBamFeaturesOption(CLI::App* const sub, const TrainModelParamPtr& params) {
  auto* const opt = sub->add_option_function<fs::path>(
      cli_opt_name::kNegBamFeatures,
      [&params](const fs::path& value) {
        ParseFileList(cli_opt_name::kNegBamFeatures, value, params->negative_features);
      },
      "List of features files for negative/normal/healthy samples");
  opt->check(cli::NonEmptyFileValidator());
}

/**
 * @brief Add negative VCF features CLI option.
 * @param sub CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
void AddNegVcfFeaturesOption(CLI::App* const sub, const TrainModelParamPtr& params) {
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
  opt->check(cli::NonEmptyFileValidator());
}

/**
 * @brief Add CLI options common to non-germline workflows.
 * @param sub CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 * @param defaults Workflow-specific default configuration
 */
void AddNonGermlineSpecificOptions(CLI::App* const sub, const TrainModelParamPtr& params, const SVCConfig& defaults) {
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
 * @brief Helper function to add germline WGS specific CLI options.
 * @param app CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
void AddGermlineSpecificOptions(const CLI::App* const app, const TrainModelParamPtr& params) {
  for (const Workflow workflow : {kGermlineWgs}) {
    CLI::App* const sub = app->get_subcommand(enum_util::FormatEnumName(workflow));
    // CLI11 needs values while options are registered. The duplex profile is only a display placeholder; PreCallback
    // replaces it with the profile selected by the required concrete library type before training starts.
    const auto defaults = SVCConfig(ConfigProfile::kGermlineWgsDuplex);
    AddPosVcfFeaturesOption(sub, params);
    AddTruthVcfsOption(sub, params);

    cli::AddEnumOption(sub,
                       cli_opt_name::kModelType,
                       params->germline_model_type,
                       "Optional germline WGS model family whose training defaults should be used. Defaults to 'large'",
                       GermlineModelType::kLarge);

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

    sub->parse_complete_callback([sub, params] { SyncGermlineTrainingRenderedDefaults(sub, params); });
  }
}

/**
 * @brief Helper function to add tumor-only-te specific CLI options.
 * @param app CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
void AddTumorOnlyTeSpecificOptions(const CLI::App* const app, const TrainModelParamPtr& params) {
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
      "Uncompressed or gzip/BGZF-compressed (.gz) BED file of 0-based regions to exclude from model training");
  blocklist_opt->check(cli::BedFileValidator());
}

/**
 * @brief Add a truth-VCF list option for a specific class (somatic or germline) in tumor-only-wgs.
 * @param sub CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 * @param option_name CLI option name (somatic or germline truth VCFs)
 * @param truth_vcfs Destination vector for the parsed file list
 * @param description CLI help description
 */
void AddClassTruthVcfsOption(CLI::App* const sub,
                             const TrainModelParamPtr& params,
                             const std::string& option_name,
                             vec<fs::path>& truth_vcfs,
                             const std::string& description) {
  auto* const opt =
      sub->add_option_function<fs::path>(
             option_name,
             [&params, option_name, &truth_vcfs](const fs::path& value) {
               ParseFileList(option_name, value, truth_vcfs);
               if (!params->positive_features.empty() && params->positive_features.size() != truth_vcfs.size()) {
                 throw CLI::ValidationError(fmt::format(
                     "{} and {} must contain the same number of files", option_name, cli_opt_name::kPosBamFeatures));
               }
             },
             description)
          ->required();
  opt->check(cli::NonEmptyFileValidator());
}

/**
 * @brief Helper function to add tumor-only-wgs specific CLI options.
 *
 * Tumor-only-wgs trains a single model from single-sample features. Positives are drawn from two truth sets
 * (somatic and germline), each subsampled independently; negatives are everything not in the union of the truth
 * sets, subsampled by the negative ratio. A single-sample alt-support gate (`--min-support`) applies to both
 * classes. The Panel of Normals VCF is optional and is supplied to `compute_vcf_features`, not here; when absent,
 * PON features default to 0.
 *
 * @param app CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
void AddTumorOnlyWgsSpecificOptions(const CLI::App* const app, const TrainModelParamPtr& params) {
  CLI::App* const sub = app->get_subcommand(enum_util::FormatEnumName(kTumorOnlyWgs));
  const auto defaults = SVCConfig(kTumorOnlyWgs);
  AddPosVcfFeaturesOption(sub, params);
  AddNegBamFeaturesOption(sub, params);
  AddNegVcfFeaturesOption(sub, params);
  AddNonGermlineSpecificOptions(sub, params, defaults);

  AddClassTruthVcfsOption(sub,
                          params,
                          cli_opt_name::kSomaticTruthVcfs,
                          params->somatic_truth_vcfs,
                          "List of somatic truth VCFs for positive samples (one per positive sample)");
  AddClassTruthVcfsOption(sub,
                          params,
                          cli_opt_name::kGermlineTruthVcfs,
                          params->germline_truth_vcfs,
                          "List of germline truth VCFs for positive samples (one per positive sample)");

  sub->add_option(cli_opt_name::kSomaticPositiveRatio,
                  params->somatic_positive_ratio,
                  "Fraction of eligible somatic positives to keep (0.0-1.0)")
      ->default_val(params->somatic_positive_ratio)
      ->check(CLI::Range(0.0F, 1.0F));

  sub->add_option(cli_opt_name::kGermlinePositiveRatio,
                  params->germline_positive_ratio,
                  "Fraction of eligible germline positives to keep (0.0-1.0); tune down to balance classes")
      ->default_val(params->germline_positive_ratio)
      ->check(CLI::Range(0.0F, 1.0F));

  sub->add_option(
         cli_opt_name::kNegativeRatio, params->negative_ratio, "Fraction of eligible negatives to keep (0.0-1.0)")
      ->default_val(params->negative_ratio)
      ->check(CLI::Range(0.0F, 1.0F));

  sub->add_option(cli_opt_name::kMinSupport,
                  params->min_support,
                  "Minimum alt read support for a variant to be eligible for training")
      ->default_val(defaults.min_support);
}

/**
 * @brief Helper function to add tumor-normal-wgs specific CLI options.
 * @param app CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
void AddTumorNormalWgsSpecificOptions(const CLI::App* const app, const TrainModelParamPtr& params) {
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
void AddGermlineTaggingSpecificOptions(const CLI::App* const app, const TrainModelParamPtr& params) {
  CLI::App* const sub = app->get_subcommand(enum_util::FormatEnumName(kGermlineTagging));
  const auto defaults = SVCConfig(kGermlineTagging);
  AddPosVcfFeaturesOption(sub, params);
  AddTruthVcfsOption(sub, params);
  AddNonGermlineSpecificOptions(sub, params, defaults);
}
}  // namespace

void train_model::DefineOptions(CLI::App* const app, const TrainModelParamPtr& params) {
  // Add a subcommand for each workflow
  for (const Workflow workflow : kSupportedWorkflows) {
    const std::string name = enum_util::FormatEnumName(workflow);
    const std::string desc = fmt::format("Train model for the {} workflow", name);
    CLI::App* const sub = app->add_subcommand(name, desc)->fallthrough();
    // Do not apply force_callback() to subcommand options to avoid overwriting params set by other subcommands
    AddCommonOptions(sub, params);
    // Germline defaults cannot be final until --library-type is parsed, so use duplex only while registering options.
    const auto defaults = workflow == kGermlineWgs ? SVCConfig(ConfigProfile::kGermlineWgsDuplex) : SVCConfig(workflow);
    cli::AddEnumOption(sub,
                       cli_opt_name::kNormalizeFeatures,
                       params->normalize_features,
                       "Normalize read depth related features values",
                       defaults.normalize_features);
  }
  app->require_subcommand(kMinSubcommands, kMaxSubcommands);

  AddExplicitLibraryTypeOption(app->get_subcommand(enum_util::FormatEnumName(kGermlineWgs)),
                               params->library_type,
                               "Library preparation type (required because this command has no BAM input)")
      ->required();

  // Add workflow-specific CLI options for each subcommand
  AddGermlineSpecificOptions(app, params);
  AddTumorOnlyTeSpecificOptions(app, params);
  AddTumorOnlyWgsSpecificOptions(app, params);
  AddTumorNormalWgsSpecificOptions(app, params);
  AddGermlineTaggingSpecificOptions(app, params);

  // hide the `germline-tagging` subcommand by assigning it to an empty string group
  app->get_subcommand(enum_util::FormatEnumName(kGermlineTagging))->group("");
}

namespace {
/**
 * @brief Helper function to set default output model paths if not provided via CLI.
 * @param sub Subcommand CLI application pointer to check which options were specified via CLI
 * @param params CLI parameters shared pointer to apply defaults to
 */
void SetDefaultOutputModelPaths(const CLI::App* const sub, const TrainModelParamPtr& params) {
  if (params->workflow == kGermlineWgs) {
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
void ApplyConfig(const CLI::App* const sub, const TrainModelParamPtr& params) {
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
  if (sub->get_option_no_throw(cli_opt_name::kMinSupport) == nullptr || sub->count(cli_opt_name::kMinSupport) == 0) {
    params->min_support = params->config.min_support;
  }
}

/**
 * @brief Write a resolved value into an option's rendered default when the user did not provide it explicitly.
 * @details CLI11 stores option defaults at registration time. Germline WGS training defaults depend on both
 * `--library-type` and `--model-type`, so they are only final after the pre-callback resolves the config profile.
 * @param sub Subcommand that owns the option.
 * @param option_name CLI option name to update.
 * @param value Resolved value to render in provenance output.
 */
void SetRenderedDefault(CLI::App* const sub, const std::string_view option_name, const std::string& value) {
  auto* const opt = sub->get_option_no_throw(std::string{option_name});
  if (opt != nullptr && opt->empty()) {
    opt->default_str(value);
  }
}

/**
 * @brief Resolve the germline WGS config selected by parsed CLI options.
 * @details This is used in two phases: parse completion, before the common "Running:" log renders defaults, and the
 * module pre-callback, before training starts. Keeping both phases on the same helper prevents provenance-only defaults
 * from drifting from the effective training configuration.
 * @param params Parsed training CLI parameters.
 * @return Config profile with selected germline model-family training defaults applied.
 */
SVCConfig ResolveGermlineTrainingConfig(const TrainModelParamPtr& params) {
  // The CLI accepts only concrete library types for training. The empty BAM list keeps programmatic callers from
  // silently resolving the auto-detect sentinel when no BAM input exists.
  params->library_type = ResolveLibraryType(params->library_type, {});
  auto config = JsonToConfig(params->config_file.value_or(fs::path{}), GermlineWgsConfigProfile(params->library_type));
  config.ApplyGermlineModelTrainingDefaults(params->germline_model_type);
  return config;
}

/**
 * @brief Synchronize resolved germline WGS training defaults with rendered command-line provenance.
 * @param sub Active germline WGS subcommand.
 * @param params Training parameters after config application.
 */
void SyncGermlineTrainingOptions(CLI::App* const sub, const TrainModelParamPtr& params) {
  if (params->workflow != kGermlineWgs) {
    return;
  }
  SetRenderedDefault(sub, cli_opt_name::kModelType, enum_util::FormatEnumName(params->germline_model_type));
  SetRenderedDefault(sub, cli_opt_name::kSnvIterations, std::to_string(params->snv_iterations));
  SetRenderedDefault(sub, cli_opt_name::kIndelIterations, std::to_string(params->indel_iterations));
}

/**
 * @brief Synchronize germline WGS rendered defaults before common CLI logging.
 * @details The common root parse-complete callback emits the "Running:" log before module pre-callbacks run. This
 * subcommand parse-complete callback updates only CLI11's rendered default strings so that the log line reflects
 * model-type-specific iteration defaults. Explicit user-provided iteration values are left untouched.
 * @param sub Active germline WGS subcommand.
 * @param params Parsed training parameters.
 */
void SyncGermlineTrainingRenderedDefaults(CLI::App* const sub, const TrainModelParamPtr& params) {
  const auto config = ResolveGermlineTrainingConfig(params);
  SetRenderedDefault(sub, cli_opt_name::kModelType, enum_util::FormatEnumName(params->germline_model_type));
  SetRenderedDefault(sub, cli_opt_name::kSnvIterations, std::to_string(config.snv_iterations));
  SetRenderedDefault(sub, cli_opt_name::kIndelIterations, std::to_string(config.indel_iterations));
}
}  // namespace

void train_model::PreCallback(const CLI::App* const app, const TrainModelParamPtr& params) {
  // Check which subcommand was used, set workflow and config accordingly, and apply config defaults as needed
  for (const Workflow workflow : kSupportedWorkflows) {
    const std::string name = enum_util::FormatEnumName(workflow);
    if (app->got_subcommand(name)) {
      params->workflow = workflow;
      if (workflow == kGermlineWgs) {
        params->config = ResolveGermlineTrainingConfig(params);
      } else {
        params->config = JsonToConfig(params->config_file.value_or(fs::path{}), params->workflow);
      }
      auto* const sub = app->get_subcommand(name);
      ApplyConfig(sub, params);
      SyncGermlineTrainingOptions(sub, params);
      break;
    }
  }
  params->command_line = cli::GetCommandLineInfo(app);
}

}  // namespace xoos::svc
