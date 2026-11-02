#include "merge-and-label-bed/merge-and-label-bed-cli.h"

#include <CLI/CLI.hpp>

#include <xoos/cli/cli.h>
#include <xoos/cli/validators/bed-validator.h>

#include "merge-and-label-bed/merge-and-label-bed-options.h"

namespace xoos::alignment_metrics::merge_and_label_bed {

const std::string kInputOptionsGroup = "Input Options";
const std::string kOutputOptionsGroup = "Output Options";

void SetMergeAndLabelBedCommandLineInfo(const CLI::App* const app, const MergeAndLabelBedOptionsPtr& options) {
  options->command_line_info = cli::GetCommandLineInfo(app);
  // Pair each --bed with its --label by position. Enforced here because CLI11 cannot express the
  // cross-option "equal count" constraint declaratively.
  if (options->beds.size() != options->labels.size()) {
    throw cli::ValidationError("--bed and --label must be given the same number of times ({} bed(s) vs {} label(s))",
                               options->beds.size(),
                               options->labels.size());
  }
  options->inputs.clear();
  options->inputs.reserve(options->beds.size());
  for (size_t i = 0; i < options->beds.size(); ++i) {
    options->inputs.emplace_back(options->beds.at(i), options->labels.at(i));
  }
}

void DefineMergeAndLabelBedOptions(CLI::App* const app, const MergeAndLabelBedOptionsPtr& options) {
  // Input Options
  app->add_option("--bed",
                  options->beds,
                  "Input BED file(s). Repeat once per label; pairs with --label by position. Each file may be "
                  "gzip/bgzip compressed when named with a `.gz` extension.")
      ->required()
      ->check(CLI::ExistingFile)
      ->check(cli::BedFileValidator())
      ->group(kInputOptionsGroup);

  app->add_option("--label", options->labels, "Label assigned to the regions of the corresponding --bed file.")
      ->required()
      ->group(kInputOptionsGroup);

  // Output Options
  app->add_option("--output", options->output, "Output labeled BED file path.")->required()->group(kOutputOptionsGroup);

  app->add_option<bool>("--overwrite", options->overwrite, "Overwrite the output file if it already exists.")
      ->expected(0)
      ->default_val("false")
      ->group(kOutputOptionsGroup);
}

}  // namespace xoos::alignment_metrics::merge_and_label_bed
