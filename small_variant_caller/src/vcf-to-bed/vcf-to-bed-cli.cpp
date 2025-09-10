#include "vcf-to-bed-cli.h"

#include <string>

#include <xoos/cli/thread-count-option-util.h>

#include "util/cli-util.h"

namespace xoos::svc {

using cli::AddThreadCountOption;

void DefineOptionsVcfToBed(cli::AppPtr app, const VcfToBedParamPtr& params) {
  app->add_option("--vcf-input", params->vcf_file, "Path to input VCF file")->required()->check(kCliNonEmptyFile);
  app->add_option("--output-file", params->output_bed_file, "Output BED file name")
      ->default_val(kDefaultVcfToBedFileName);
  app->add_option("--left-pad", params->left_pad, "left-padding to variant start position")
      ->default_val(kDefaultLeftPad)
      ->check(CLI::NonNegativeNumber);
  app->add_option("--right-pad", params->right_pad, "right-padding to variant end position")
      ->default_val(kDefaultRightPad)
      ->check(CLI::NonNegativeNumber);
  app->add_option(
         "--collapse-distance", params->collapse_dist, "distance to collapse nearby variants after padding is applied")
      ->default_val(kDefaultCollapsableDist)
      ->check(CLI::NonNegativeNumber);
  app->add_option_function<std::string>(
         "--region",
         [&params](const std::string& value) { params->chrom_intervals_map = GetChromIntervalMap(value); },
         "Path to a BED file of target regions")
      ->check(kCliNonEmptyFile);
  AddThreadCountOption(app, "--threads", params->threads);
  AddWarnAsErrorOption(app);
}

}  // namespace xoos::svc
