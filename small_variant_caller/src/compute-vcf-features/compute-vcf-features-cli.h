#pragma once

#include <memory>

#include <xoos/cli/cli.h>

#include "compute-vcf-features.h"

namespace xoos::svc {
using ComputeVcfFeaturesParamPtr = std::shared_ptr<ComputeVcfFeaturesParam>;
void DefineOptionsComputeVcfFeatures(cli::AppPtr app, const ComputeVcfFeaturesParamPtr& params);
}  // namespace xoos::svc
