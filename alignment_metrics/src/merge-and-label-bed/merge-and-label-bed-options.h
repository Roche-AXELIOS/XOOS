#pragma once

#include <memory>
#include <string>
#include <utility>

#include <xoos/io/metadata-util.h>
#include <xoos/types/fs.h>
#include <xoos/types/vec.h>

namespace xoos::alignment_metrics::merge_and_label_bed {

using BedAndLabelPair = std::pair<fs::path, std::string>;

struct MergeAndLabelBedOptions {
  // Raw CLI inputs: `beds[i]` is assigned the label `labels[i]`. Paired into `inputs` during the
  // pre-callback after validating that the two lists have equal length.
  vec<fs::path> beds;
  vec<std::string> labels;

  // Each entry is a BED file paired with the label to assign to all of its regions.
  vec<BedAndLabelPair> inputs;

  // Output options
  fs::path output;
  bool overwrite{false};

  // Program metadata
  io::CommandLineInfo command_line_info{};
};

using MergeAndLabelBedOptionsPtr = std::shared_ptr<MergeAndLabelBedOptions>;

}  // namespace xoos::alignment_metrics::merge_and_label_bed
