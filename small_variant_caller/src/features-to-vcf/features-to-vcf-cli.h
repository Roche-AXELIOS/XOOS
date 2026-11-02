#pragma once

#include <memory>

#include <xoos/cli/cli.h>

#include "features-to-vcf.h"

namespace xoos::svc {

using FeaturesToVcfParamPtr = std::shared_ptr<FeaturesToVcfParam>;

void DefineOptionsFeaturesToVcf(CLI::App* app, const FeaturesToVcfParamPtr& params);

}  // namespace xoos::svc
