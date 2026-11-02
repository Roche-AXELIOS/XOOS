#pragma once

#include <xoos/cli/cli.h>

#include "rescue-secondary/rescue-secondary-options.h"

namespace xoos::read_collapser::rescue_secondary {

// Called via DefineOptions<T>/PreCallback<T> std::function alias; SonarQube cannot trace indirect calls.
void SetRescueSecondaryCommandLineInfo(  // NOSONAR cpp:S5536
    const CLI::App* app,
    const RescueSecondaryOptionsPtr& options);

// Called via DefineOptions<T>/PreCallback<T> std::function alias; SonarQube cannot trace indirect calls.
void DefineRescueSecondaryOptions(  // NOSONAR cpp:S5536
    CLI::App* app,
    const RescueSecondaryOptionsPtr& options);
}  // namespace xoos::read_collapser::rescue_secondary
