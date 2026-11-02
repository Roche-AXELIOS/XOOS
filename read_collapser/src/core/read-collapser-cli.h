#pragma once

#include <htslib/sam.h>

#include <CLI/Formatter.hpp>

#include <xoos/cli/cli.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "core/read-collapser-options.h"

namespace xoos::read_collapser {

void SetCommandLineInfo(const CLI::App* app, const ReadCollapserOptionsPtr& options);

void DefineConsensusOptions(CLI::App* app, const ReadCollapserOptionsPtr& options);

void DefineMarkdupOptions(CLI::App* app, const ReadCollapserOptionsPtr& options);

// Called via DefineOptions<T>/PreCallback<T> std::function alias; SonarQube cannot trace indirect calls.
void DefineOptions(  // NOSONAR cpp:S5536
    CLI::App* app,
    ReadCollapserOptionsPtr& options);
}  // namespace xoos::read_collapser
