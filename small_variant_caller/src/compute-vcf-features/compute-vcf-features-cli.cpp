#include "compute-vcf-features-cli.h"

#include <xoos/cli/enum-option-util.h>
#include <xoos/cli/thread-count-option-util.h>

namespace xoos::svc {

using cli::AddOptionalEnumOption;
using cli::AddThreadCountOption;
using enum Workflow;

constexpr std::array kSupportedWorkflows = {
    kGermline, kGermlineMultiSample, kTumorOnlyTe, kTumorNormalWgs, kGermlineTagging, kCustom};

// CLI option default values
constexpr auto* kDefaultVcfFeaturesFileName = "vcf_features.txt";

/**
 * @brief Helper function to define CLI options for the main application.
 * These options are either required or important in this submodule, or they are common to all other submodules in SVC.
 * @param app Main application pointer where CLI options will be defined.
 * @param params Shared pointer to CLI parameters to store parsed option values
 */
static void AddMainOptions(CLI::App* const app, const ComputeVcfFeaturesParamPtr& params) {
  AddWarnAsErrorOption(app);

  AddThreadCountOption(app, cli_opt_name::kThreads, params->threads);

  app->add_option(cli_opt_name::kConfig, params->config_file, "Path to config JSON file")->check(CLI::ExistingFile);

  auto* const vcf_opt = app->add_option(cli_opt_name::kVcfInput,
                                        params->vcf_file,
                                        "Path to input VCF file, produced by GATK Mutect2/HaplotypeCaller")
                            ->required();
  CheckIndexedVcfFile(vcf_opt);

  auto* const genome_opt =
      app->add_option(cli_opt_name::kGenome, params->genome, "Path to indexed FASTA file for reference genome")
          ->required();
  CheckIndexedFastaFile(genome_opt);

  auto* const regions_opt = app->add_option_function<fs::path>(
      cli_opt_name::kTargetRegions,
      [&params](const fs::path& value) { params->target_regions = GetChromIntervalMap(value); },
      "Path to BED file for 0-based target regions");
  CheckBedFile(regions_opt);

  app->add_option(cli_opt_name::kOutputFile, params->output_file, "Path to output features file")
      ->default_val(kDefaultVcfFeaturesFileName)
      ->check(CLI::NonexistentPath);
}

/**
 * @brief Helper function to add core CLI options for subcommands. CLI options added here are common across all
 * subcommands.
 * @param sub Subcommand application pointer where CLI options are to be added
 * @param params Shared pointer to CLI parameters to store parsed option values
 */
static void AddCoreOptions(CLI::App* const sub, const ComputeVcfFeaturesParamPtr& params) {
  // Note there are no workflow-specific defaults/options at this time
  sub->add_option(cli_opt_name::kOutputBed, params->output_bed, "Path to output BED file for extracted features")
      ->check(CLI::NonexistentPath);

  sub->add_option(
         cli_opt_name::kLeftPad, params->left_pad, "left-padding to variant start position for output BED file")
      ->default_val(kDefaultLeftPad)
      ->check(CLI::NonNegativeNumber);

  sub->add_option(
         cli_opt_name::kRightPad, params->right_pad, "right-padding to variant end position for output BED file")
      ->default_val(kDefaultRightPad)
      ->check(CLI::NonNegativeNumber);

  sub->add_option(cli_opt_name::kCollapseDistance,
                  params->collapse_dist,
                  "distance to collapse nearby variants after padding is applied for output BED file")
      ->default_val(kDefaultCollapseDistance)
      ->check(CLI::NonNegativeNumber);
}

void compute_vcf_features::DefineOptions(CLI::App* const app, ComputeVcfFeaturesParamPtr& params) {
  AddMainOptions(app, params);

  // Add a subcommand for each workflow
  for (const Workflow workflow : kSupportedWorkflows) {
    const std::string name = enum_util::FormatEnumName(workflow);
    const std::string desc = fmt::format("Compute VCF features for the {} workflow", name);
    CLI::App* const sub = app->add_subcommand(name, desc)->fallthrough();
    // Do not apply force_callback() to subcommand options to avoid overwriting params set by other subcommands
    AddCoreOptions(sub, params);
    AddSharedOptions(sub, params);
  }
  app->require_subcommand(kMinSubcommands, kMaxSubcommands);

  // hide the `germline-tagging` subcommand by assigning it to an empty string group
  app->get_subcommand(enum_util::FormatEnumName(kGermlineTagging))->group("");
}

void compute_vcf_features::PreCallback(const cli::ConstAppPtr app, const ComputeVcfFeaturesParamPtr& params) {
  params->command_line = GetCommandLineInfo(app);
  // Check which subcommand was used, set workflow and config accordingly
  for (const Workflow workflow : kSupportedWorkflows) {
    const auto workflow_name = enum_util::FormatEnumName(workflow);
    if (app->got_subcommand(workflow_name)) {
      params->workflow = workflow;
      params->config = JsonToConfig(params->config_file.value_or(fs::path{}), workflow);
      // No config defaults to apply here
      break;
    }
  }
}

}  // namespace xoos::svc
