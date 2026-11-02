#include "rescue-secondary/rescue-secondary-cli.h"

#include <CLI/CLI.hpp>

#include <xoos/cli/cli.h>
#include <xoos/cli/thread-count-option-util.h>
#include <xoos/types/float.h>

#include "rescue-secondary/rescue-secondary-options.h"

namespace xoos::read_collapser::rescue_secondary {

const std::string kInputOptionsGroup = "Input Options";
const std::string kOutputOptionsGroup = "Output Options";
const std::string kRescueOptionsGroup = "Rescue Options";
const std::string kPerformanceOptionsGroup = "Performance Options";

const f64 kDefaultMinAlignmentScoreRatio = 0.8;

void SetRescueSecondaryCommandLineInfo(const CLI::App* const app, const RescueSecondaryOptionsPtr& options) {
  options->command_line_info = cli::GetCommandLineInfo(app);
}

void DefineRescueSecondaryOptions(CLI::App* const app, const RescueSecondaryOptionsPtr& options) {
  // Input Options
  app->add_option("--bam", options->bam, "Input BAM path or '-' for stdin.")
      ->default_val(kStdin)
      ->group(kInputOptionsGroup);

  // Output Options
  app->add_option("--output-dir", options->output_dir, "Directory for output file(s) and metrics.")
      ->default_val(".")
      ->group(kOutputOptionsGroup);

  app->add_option("--output-prefix",
                  options->output_prefix,
                  "Prefix prepended to every output filename, joined with an underscore "
                  "(e.g. 'SAMPLE' -> 'SAMPLE_<filename>'). Defaults to no prefix.")
      ->group(kOutputOptionsGroup);

  app->add_option<bool>("--overwrite", options->overwrite, "Overwrite existing output files.")
      ->expected(0)
      ->default_val("false")
      ->group(kOutputOptionsGroup);

  // Rescue Options
  app->add_option<bool>("--collated",
                        options->collated,
                        "Input is collated (records grouped by read name for single-pass processing).")
      ->expected(0)
      ->default_val("false")
      ->group(kRescueOptionsGroup);

  app->add_option("--min-alignment-score-ratio",
                  options->min_alignment_score_ratio,
                  "Minimum AS/read_length ratio for a secondary alignment to be eligible for rescue.")
      ->default_val(kDefaultMinAlignmentScoreRatio)
      ->check(CLI::NonNegativeNumber)
      ->group(kRescueOptionsGroup);

  // Performance Options
  cli::AddThreadCountOption(app, "--threads", options->threads)->group(kPerformanceOptionsGroup);
}

}  // namespace xoos::read_collapser::rescue_secondary
