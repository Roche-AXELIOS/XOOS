#include "compute-vcf-features-cli.h"

#include <xoos/cli/enum-option-util.h>
#include <xoos/cli/thread-count-option-util.h>
#include <xoos/cli/validators/bed-validator.h>
#include <xoos/cli/validators/indexed-fasta-validator.h>
#include <xoos/cli/validators/indexed-vcf-validator.h>

#include "util/cli-util.h"

namespace xoos::svc {

using cli::AddThreadCountOption;
using enum Workflow;

constexpr std::array kSupportedWorkflows = {
    kGermlineWgs, kGermlineMitochondria, kTumorOnlyTe, kTumorOnlyWgs, kTumorNormalWgs, kGermlineTagging, kCustom};

// CLI option default values
constexpr auto* kDefaultVcfFeaturesFileName = "vcf_features.txt";

namespace {
/**
 * @brief Helper function to add core CLI options for subcommands. CLI options added here are common across all
 * subcommands.
 * @param sub Subcommand application pointer where CLI options are to be added
 * @param params Shared pointer to CLI parameters to store parsed option values
 */
void AddCommonOptions(CLI::App* const sub, const ComputeVcfFeaturesParamPtr& params) {
  AddWarnAsErrorOption(sub);

  AddThreadCountOption(sub, cli_opt_name::kThreads, params->threads);

  sub->add_option(cli_opt_name::kConfig, params->config_file, "Path to config JSON file")->check(CLI::ExistingFile);

  auto* const vcf_opt = sub->add_option(cli_opt_name::kVcfInput,
                                        params->vcf_file,
                                        "Path to input VCF file, produced by GATK Mutect2/HaplotypeCaller")
                            ->required();
  vcf_opt->check(cli::IndexedVcfFileValidator());

  auto* const genome_opt =
      sub->add_option(cli_opt_name::kGenome, params->genome, "Path to indexed FASTA file for reference genome")
          ->required();
  genome_opt->check(cli::IndexedFastaFileValidator());

  auto* const regions_opt = sub->add_option_function<fs::path>(
      cli_opt_name::kTargetRegions,
      [&params](const fs::path& value) { params->target_regions = GetChromIntervalMap(value); },
      "Path to uncompressed or gzip/BGZF-compressed (.gz) BED file for 0-based target regions");
  regions_opt->check(cli::BedFileValidator());

  sub->add_option(cli_opt_name::kOutputFile, params->output_file, "Path to output features file")
      ->default_val(kDefaultVcfFeaturesFileName)
      ->check(CLI::NonexistentPath);

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
}  // namespace

void compute_vcf_features::DefineOptions(CLI::App* const app, ComputeVcfFeaturesParamPtr& params) {
  // Add a subcommand for each workflow
  for (const Workflow workflow : kSupportedWorkflows) {
    const std::string name = enum_util::FormatEnumName(workflow);
    const std::string desc = fmt::format("Compute VCF features for the {} workflow", name);
    CLI::App* const sub = app->add_subcommand(name, desc)->fallthrough();
    // Do not apply force_callback() to subcommand options to avoid overwriting params set by other subcommands
    AddCommonOptions(sub, params);
    AddSharedOptions(sub, params);
  }
  app->require_subcommand(kMinSubcommands, kMaxSubcommands);

  AddExplicitLibraryTypeOption(app->get_subcommand(enum_util::FormatEnumName(kGermlineWgs)),
                               params->library_type,
                               "Library preparation type (required because this command has no BAM input)")
      ->required();

  // Hide experimental subcommands by assigning them to an empty string group.
  app->get_subcommand(enum_util::FormatEnumName(kGermlineTagging))->group("");
  app->get_subcommand(enum_util::FormatEnumName(kCustom))->group("");
}

void compute_vcf_features::PreCallback(const CLI::App* const app, const ComputeVcfFeaturesParamPtr& params) {
  // Check which subcommand was used, set workflow and config accordingly
  for (const Workflow workflow : kSupportedWorkflows) {
    const auto workflow_name = enum_util::FormatEnumName(workflow);
    if (app->got_subcommand(workflow_name)) {
      params->workflow = workflow;
      if (workflow == kGermlineWgs) {
        // The CLI accepts only concrete library types for this command. The empty BAM list keeps programmatic callers
        // from silently resolving the auto-detect sentinel when no BAM input exists.
        params->library_type = ResolveLibraryType(params->library_type, {});
        params->config =
            JsonToConfig(params->config_file.value_or(fs::path{}), GermlineWgsConfigProfile(params->library_type));
      } else {
        params->config = JsonToConfig(params->config_file.value_or(fs::path{}), workflow);
      }
      // No config defaults to apply here
      break;
    }
  }
  params->command_line = cli::GetCommandLineInfo(app);
}

}  // namespace xoos::svc
