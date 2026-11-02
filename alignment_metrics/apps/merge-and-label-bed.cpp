#include "merge-and-label-bed/merge-and-label-bed.h"

#include <memory>

#include <xoos/cli/cli.h>

#include "merge-and-label-bed/merge-and-label-bed-cli.h"
#include "merge-and-label-bed/merge-and-label-bed-options.h"

int main(int argc, char** argv) {
  namespace merge_and_label_bed = xoos::alignment_metrics::merge_and_label_bed;
  xoos::cli::StandardMainParam<merge_and_label_bed::MergeAndLabelBedOptions> param = {
      .program_name = PROGRAM_NAME,
      .version = VERSION,
      .cli_opts = std::make_shared<merge_and_label_bed::MergeAndLabelBedOptions>(),
      .define_options = merge_and_label_bed::DefineMergeAndLabelBedOptions,
      .main = merge_and_label_bed::MergeAndLabelBed,
      .pre_callback = merge_and_label_bed::SetMergeAndLabelBedCommandLineInfo,
  };
  return StandardMain(argc, argv, param);
}
