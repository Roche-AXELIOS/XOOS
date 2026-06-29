#pragma once
#include <xoos/cli/cli.h>
#include <xoos/types/float.h>

#include "demux-strand-index/demux-strand-build-pipeline.h"

namespace xoos::demux::strand {

using StrandBuildParamPtr = std::shared_ptr<StrandBuildParam>;

constexpr size_t kDefaultStrandBuildIndexKmerSize = 19;
constexpr f32 kDefaultStrandBuildIndexMaxFalsePositiveRate = 0.015625f;

/// Define the command line options and wire the up to the parameter object
void DefineOptions(cli::AppPtr app, StrandBuildParamPtr& param);

}  // namespace xoos::demux::strand
