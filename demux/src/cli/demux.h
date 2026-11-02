#pragma once
#include <memory>

#include <xoos/cli/cli.h>

#include "core/demux-and-trim-pipeline.h"

namespace xoos::demux {

using DemuxAndTrimParamPtr = std::shared_ptr<DemuxAndTrimParam>;

struct DemuxCliCallbacks {
  cli::DefineOptions<DemuxAndTrimParam> define_options;
  cli::PreCallback<DemuxAndTrimParam> pre_callback;
};

// Create CLI callbacks that share the state needed to process deprecated option names.
// Used by apps/demux.cpp and test utilities, which are outside Sonar's source scope.
DemuxCliCallbacks CreateCliCallbacks();  // NOSONAR cpp:S5536

/// Expands the input file list by recursively finding all sequence files in the directories specified in the input file
/// list. If an input file is not a directory, it is included in the expanded input file list as is.
std::vector<fs::path> ExpandInputFileList(const std::vector<fs::path>& input_files);
}  // namespace xoos::demux
