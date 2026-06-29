#include "compute-bam-features-cli.h"

#include <xoos/cli/thread-count-option-util.h>

namespace xoos::svc {

using cli::AddOptionalEnumOption;
using cli::AddThreadCountOption;
using enum Workflow;

constexpr std::array kSupportedWorkflows = {
    kGermline, kGermlineMultiSample, kTumorOnlyTe, kTumorNormalWgs, kGermlineTagging, kCustom};

// CLI option default values
constexpr auto* kDefaultBamFeaturesFileName = "bam_features.txt";

/**
 * @brief Helper function to define CLI options for the main application.
 * These options are either required or important in this submodule, or they are common to all other submodules in SVC.
 * @param app Main application pointer where CLI options will be defined.
 * @param params Shared pointer to CLI parameters to store parsed option values
 */
static void AddMainOptions(CLI::App* const app, const ComputeBamFeaturesCliParamsPtr& params) {
  AddWarnAsErrorOption(app);

  AddThreadCountOption(app, cli_opt_name::kThreads, params->threads);

  app->add_option(cli_opt_name::kConfig, params->config_file, "Path to config JSON file");

  auto* const bam_opt = app->add_option(cli_opt_name::kBamInput,
                                        params->bam_input,
                                        "Path(s) to input BAM file(s), produced by GATK HaplotypeCaller/Mutect2")
                            ->required();
  CheckIndexedBamFile(bam_opt);

  auto* const genome_opt =
      app->add_option(cli_opt_name::kGenome, params->genome, "Path to indexed FASTA file for reference genome")
          ->required();
  CheckIndexedFastaFile(genome_opt);

  auto* const regions_opt = app->add_option_function<fs::path>(
      cli_opt_name::kTargetRegions,
      [&params](const fs::path& value) { params->bed_regions = GetBedRegions(value); },
      "Path to BED file for 0-based target regions");
  CheckBedFile(regions_opt);

  app->add_option(cli_opt_name::kOutputFile, params->output_file, "Path to output features file")
      ->default_val(kDefaultBamFeaturesFileName)
      ->check(CLI::NonexistentPath);
}

/**
 * @brief Helper function to add core CLI options for subcommands. CLI options added here are common across all
 * subcommands.
 * @param sub Subcommand application pointer where CLI options are to be added
 * @param params Shared pointer to CLI parameters to store parsed option values
 */
static void AddCoreOptions(CLI::App* const sub, const ComputeBamFeaturesCliParamsPtr& params) {
  sub->add_option(
         cli_opt_name::kMaxRegionSizePerThread, params->max_region_size_per_thread, "Maximum region size per thread")
      ->default_val(kDefaultMaxBamRegionSizePerThread)
      ->check(CLI::PositiveNumber);
}

/**
 * @brief Helper function to add tumor-normal-wgs specific CLI options.
 * @param app CLI application pointer where options are to be added
 * @param params CLI parameters shared pointer to store option values
 */
static void AddTumorNormalWgsSpecificOptions(const CLI::App* const app, const ComputeBamFeaturesCliParamsPtr& params) {
  CLI::App* const sub = app->get_subcommand(enum_util::FormatEnumName(kTumorNormalWgs));
  sub->add_option(cli_opt_name::kTumorSampleName, params->tumor_sample_name, "tumor sample name for read groups")
      ->required();
}

void compute_bam_features::DefineOptions(CLI::App* const app, ComputeBamFeaturesCliParamsPtr& params) {
  AddMainOptions(app, params);

  // Add a subcommand for each workflow
  for (const Workflow workflow : kSupportedWorkflows) {
    const std::string name = enum_util::FormatEnumName(workflow);
    const std::string desc = fmt::format("Compute BAM features for the {} workflow", name);
    CLI::App* const sub = app->add_subcommand(name, desc)->fallthrough();
    // Do not apply force_callback() to subcommand options to avoid overwriting params set by other subcommands
    AddCoreOptions(sub, params);
    const auto defaults = SVCConfig(workflow);
    AddSharedOptions(sub, params, defaults);
  }
  app->require_subcommand(kMinSubcommands, kMaxSubcommands);

  AddTumorNormalWgsSpecificOptions(app, params);

  // hide the `germline-tagging` subcommand by assigning it to an empty string group
  app->get_subcommand(enum_util::FormatEnumName(kGermlineTagging))->group("");
}

void compute_bam_features::PreCallback(const cli::ConstAppPtr app, const ComputeBamFeaturesCliParamsPtr& params) {
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

  if (params->decode_yc != YcDecodeMethod::kNone && !IsDuplexProtocol(params->sequencing_protocol)) {
    throw CLI::ValidationError("Duplex sequencing protocol must be used when decoding YC tags");
  }
}

}  // namespace xoos::svc
