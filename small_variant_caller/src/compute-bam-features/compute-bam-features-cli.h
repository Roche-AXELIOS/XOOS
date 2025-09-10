#pragma once

#include <memory>

#include <xoos/cli/cli.h>

#include "compute-bam-features.h"

namespace xoos::svc {
using ComputeBamFeaturesCliParamsPtr = std::shared_ptr<ComputeBamFeaturesCliParams>;
void DefineOptions(cli::AppPtr app, const ComputeBamFeaturesCliParamsPtr& params);
}  // namespace xoos::svc
