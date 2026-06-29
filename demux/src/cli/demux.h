#pragma once
#include <xoos/cli/cli.h>

#include <memory>

#include "core/demux-and-trim-pipeline.h"

namespace xoos::demux {

using DemuxAndTrimParamPtr = std::shared_ptr<DemuxAndTrimParam>;

/// Define the command line options and wire the up to the parameter object
void DefineOptions(cli::AppPtr app, DemuxAndTrimParamPtr& param);

/// Validate options after parsing
void ValidateOptions(const DemuxAndTrimParamPtr& param);

/// Pre-processing before main execution
void PreMainProcessing(cli::ConstAppPtr app, const DemuxAndTrimParamPtr& param);

cli::PreCallback<DemuxAndTrimParam> CreatePreCallback(const std::string& version);

void AddVersionAndCommandLineComment(const cli::ConstAppPtr& app, const DemuxAndTrimParamPtr& param,
                                     const std::string& version);

/// Expands the input file list by recursively finding all sequence files in the directories specified in the input file
/// list. If an input file is not a directory, it is included in the expanded input file list as is.
std::vector<fs::path> ExpandInputFileList(const std::vector<fs::path>& input_files);
}  // namespace xoos::demux
