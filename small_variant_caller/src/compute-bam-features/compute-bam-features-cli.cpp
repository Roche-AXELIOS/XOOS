#include "compute-bam-features-cli.h"

#include <xoos/cli/thread-count-option-util.h>
#include <xoos/cli/validators/bed-validator.h>
#include <xoos/cli/validators/indexed-bam-validator.h>
#include <xoos/cli/validators/indexed-fasta-validator.h>

namespace xoos::svc {

using cli::AddThreadCountOption;
using enum Workflow;

constexpr std::array kSupportedWorkflows = {kGermlineWgs,
                                            kGermlineMitochondria,
                                            kTumorOnlyTe,
                                            kTumorOnlyWgs,
                                            kTumorNormalWgs,
                                            kGermlineTagging,
                                            kPanelOfNormals,
                                            kCustom};

// CLI option default values
constexpr auto* kDefaultBamFeaturesFileName = "bam_features.txt";

namespace {
/**
 * @brief Helper function to add core CLI options for subcommands. CLI options added here are common across all
 * subcommands.
 * @param sub Subcommand application pointer where CLI options are to be added
 * @param params Shared pointer to CLI parameters to store parsed option values
 */
void AddCommonOptions(CLI::App* const sub, const ComputeBamFeaturesCliParamsPtr& params) {
  AddWarnAsErrorOption(sub);

  AddThreadCountOption(sub, cli_opt_name::kThreads, params->threads);

  sub->add_option(cli_opt_name::kConfig, params->config_file, "Path to config JSON file");

  auto* const bam_opt = sub->add_option(cli_opt_name::kBamInput,
                                        params->bam_input,
                                        "Path(s) to input BAM file(s), produced by GATK HaplotypeCaller/Mutect2")
                            ->required();
  bam_opt->check(cli::IndexedBamFileValidator());

  auto* const genome_opt =
      sub->add_option(cli_opt_name::kGenome, params->genome, "Path to indexed FASTA file for reference genome")
          ->required();
  genome_opt->check(cli::IndexedFastaFileValidator());

  auto* const regions_opt = sub->add_option_function<fs::path>(
      cli_opt_name::kTargetRegions,
      [&params](const fs::path& value) { params->bed_regions = GetBedRegions(value); },
      "Path to uncompressed or gzip/BGZF-compressed (.gz) BED file for 0-based target regions");
  regions_opt->check(cli::BedFileValidator());

  sub->add_option(cli_opt_name::kOutputFile, params->output_file, "Path to output features file")
      ->default_val(kDefaultBamFeaturesFileName)
      ->check(CLI::NonexistentPath);

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
void AddTumorNormalWgsSpecificOptions(const CLI::App* const app, const ComputeBamFeaturesCliParamsPtr& params) {
  CLI::App* const sub = app->get_subcommand(enum_util::FormatEnumName(kTumorNormalWgs));
  sub->add_option(cli_opt_name::kTumorSampleName, params->tumor_sample_name, "tumor sample name for read groups")
      ->required();
}

/**
 * @brief Resolve and render germline WGS profile defaults before the common CLI logs the invocation.
 * @details The root parse-complete callback writes the `Running:` log before the module pre-callback commits runtime
 *          parameters. Resolve a local profile here so the rendered command line does not retain duplex placeholders
 *          when simplex is selected.
 */
void SyncGermlineWgsRenderedDefaults(CLI::App* const sub, const ComputeBamFeaturesCliParamsPtr& params) {
  params->library_type = ResolveLibraryType(params->library_type, params->bam_input);
  const auto config =
      JsonToConfig(params->config_file.value_or(fs::path{}), GermlineWgsConfigProfile(params->library_type));
  SetRenderedDefault(sub, cli_opt_name::kLibraryType, enum_util::FormatEnumName(params->library_type));
  compute_bam_features::SyncRenderedConfigDefaults(sub, config);
}
}  // namespace

void compute_bam_features::DefineOptions(CLI::App* const app, ComputeBamFeaturesCliParamsPtr& params) {
  // Add a subcommand for each workflow
  for (const Workflow workflow : kSupportedWorkflows) {
    const std::string name = enum_util::FormatEnumName(workflow);
    const std::string desc = fmt::format("Compute BAM features for the {} workflow", name);
    CLI::App* const sub = app->add_subcommand(name, desc)->fallthrough();
    // Do not apply force_callback() to subcommand options to avoid overwriting params set by other subcommands
    AddCommonOptions(sub, params);
    // CLI11 requires defaults when options are registered, before the germline library type can be resolved. Use the
    // duplex profile as a placeholder; PreCallback applies the resolved simplex or duplex profile to unset options.
    const auto defaults = workflow == kGermlineWgs ? SVCConfig(ConfigProfile::kGermlineWgsDuplex) : SVCConfig(workflow);
    AddSharedOptions(sub, params, defaults);
  }
  app->require_subcommand(kMinSubcommands, kMaxSubcommands);

  AddTumorNormalWgsSpecificOptions(app, params);

  cli::AddEnumOption(app->get_subcommand(enum_util::FormatEnumName(kGermlineWgs)),
                     cli_opt_name::kLibraryType,
                     params->library_type,
                     "Library preparation type; auto-detect inspects BAM metadata and read content",
                     LibraryType::kAutoDetect);
  auto* const germline_sub = app->get_subcommand(enum_util::FormatEnumName(kGermlineWgs));
  germline_sub->parse_complete_callback(
      [germline_sub, params] { SyncGermlineWgsRenderedDefaults(germline_sub, params); });

  // Hide experimental subcommands by assigning them to an empty string group.
  app->get_subcommand(enum_util::FormatEnumName(kGermlineTagging))->group("");
  app->get_subcommand(enum_util::FormatEnumName(kCustom))->group("");
}

void compute_bam_features::PreCallback(const CLI::App* const app, const ComputeBamFeaturesCliParamsPtr& params) {
  // Check which subcommand was used, set workflow and config accordingly, and apply config defaults as needed
  for (const Workflow workflow : kSupportedWorkflows) {
    const std::string name = enum_util::FormatEnumName(workflow);
    if (app->got_subcommand(name)) {
      params->workflow = workflow;
      auto* const sub = app->get_subcommand(name);
      if (workflow == kGermlineWgs) {
        params->library_type = ResolveLibraryType(params->library_type, params->bam_input);
        params->config =
            JsonToConfig(params->config_file.value_or(fs::path{}), GermlineWgsConfigProfile(params->library_type));
      } else {
        params->config = JsonToConfig(params->config_file.value_or(fs::path{}), params->workflow);
      }
      ApplyConfig(sub, params);
      if (workflow == kGermlineWgs) {
        SetRenderedDefault(sub, cli_opt_name::kLibraryType, enum_util::FormatEnumName(params->library_type));
        compute_bam_features::SyncRenderedConfigDefaults(sub, params->config);
      }
      break;
    }
  }

  if (params->decode_yc != YcDecodeMethod::kNone && !IsDuplexProtocol(params->sequencing_protocol)) {
    throw CLI::ValidationError("Duplex sequencing protocol must be used when decoding YC tags");
  }
  params->command_line = cli::GetCommandLineInfo(app);
  if (params->workflow == kGermlineWgs) {
    ReplaceOrAppendCommandLineOption(params->command_line->command_line,
                                     cli_opt_name::kLibraryType,
                                     enum_util::FormatEnumName(params->library_type));
  }
}

}  // namespace xoos::svc
