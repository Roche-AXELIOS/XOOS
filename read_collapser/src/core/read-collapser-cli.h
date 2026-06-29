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

void SetCommandLineInfo(cli::ConstAppPtr app, const ReadCollapserOptionsPtr& options);

void DefineConsensusOptions(cli::AppPtr app, const ReadCollapserOptionsPtr& options);

void DefineMarkdupOptions(cli::AppPtr app, const ReadCollapserOptionsPtr& options);

void DefineOptions(cli::AppPtr app, ReadCollapserOptionsPtr& options);

}  // namespace xoos::read_collapser
