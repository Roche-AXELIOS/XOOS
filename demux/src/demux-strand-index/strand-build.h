#pragma once
#include <xoos/cli/cli.h>
#include <xoos/types/float.h>

#include "demux-strand-index/demux-strand-build-pipeline.h"

namespace xoos::demux::strand {

using StrandBuildParamPtr = std::shared_ptr<StrandBuildParam>;

constexpr size_t kDefaultStrandBuildIndexKmerSize = 19;
constexpr f32 kDefaultStrandBuildIndexMaxFalsePositiveRate = 0.015625F;

/// Define the command line options and wire the up to the parameter object
// Called via DefineOptions<T>/PreCallback<T> std::function alias; SonarQube cannot trace indirect calls.
void DefineOptions(  // NOSONAR cpp:S5536
    CLI::App* app,
    const StrandBuildParamPtr& param);
// std::function alias; SonarQube cannot trace calls through std::function
}  // namespace xoos::demux::strand
