#pragma once

#include <xoos/cli/cli.h>

#include "merge-and-label-bed/merge-and-label-bed-options.h"

namespace xoos::alignment_metrics::merge_and_label_bed {

// Called via DefineOptions<T>/PreCallback<T> std::function alias; SonarQube cannot trace indirect calls.
void SetMergeAndLabelBedCommandLineInfo(  // NOSONAR cpp:S5536
    const CLI::App* app,
    const MergeAndLabelBedOptionsPtr& options);

// Called via DefineOptions<T>/PreCallback<T> std::function alias; SonarQube cannot trace indirect calls.
void DefineMergeAndLabelBedOptions(  // NOSONAR cpp:S5536
    CLI::App* app,
    const MergeAndLabelBedOptionsPtr& options);
// SonarQube cannot trace calls through std::function
}  // namespace xoos::alignment_metrics::merge_and_label_bed
