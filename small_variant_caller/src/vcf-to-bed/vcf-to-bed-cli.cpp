#include "vcf-to-bed-cli.h"

#include <string>

#include <xoos/cli/thread-count-option-util.h>

#include "core/cli-option-names.h"
#include "util/cli-util.h"

namespace xoos::svc {

using cli::AddThreadCountOption;

// Default values for CLI options
const std::string kDefaultOutputFile = "vcf-regions.bed";

void DefineOptionsVcfToBed(CLI::App* app, const VcfToBedParamPtr& params) {
  auto* const vcf_opt =
      app->add_option(cli_opt_name::kVcfInput, params->vcf_file, "Path to input VCF file")->required();
  CheckNonEmptyFile(vcf_opt);
  app->add_option(cli_opt_name::kOutputFile, params->output_bed_file, "Output BED file name")
      ->default_val(kDefaultOutputFile)
      ->check(CLI::NonexistentPath);
  app->add_option(cli_opt_name::kLeftPad, params->left_pad, "left-padding to variant start position")
      ->default_val(kDefaultLeftPad)
      ->check(CLI::NonNegativeNumber);
  app->add_option(cli_opt_name::kRightPad, params->right_pad, "right-padding to variant end position")
      ->default_val(kDefaultRightPad)
      ->check(CLI::NonNegativeNumber);
  app->add_option(cli_opt_name::kCollapseDistance,
                  params->collapse_dist,
                  "distance to collapse nearby variants after padding is applied")
      ->default_val(kDefaultCollapseDistance)
      ->check(CLI::NonNegativeNumber);
  auto* const regions_opt = app->add_option_function<std::string>(
      cli_opt_name::kTargetRegions,
      [&params](const std::string& value) { params->chrom_intervals_map = GetChromIntervalMap(value); },
      "Path to BED file for 0-based target regions");
  CheckBedFile(regions_opt);
  AddThreadCountOption(app, cli_opt_name::kThreads, params->threads);
  AddWarnAsErrorOption(app);
}

}  // namespace xoos::svc
