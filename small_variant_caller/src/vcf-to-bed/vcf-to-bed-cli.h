#pragma once

#include <memory>

#include <xoos/cli/cli.h>

#include "vcf-to-bed.h"

namespace xoos::svc {
using VcfToBedParamPtr = std::shared_ptr<VcfToBedParam>;
void DefineOptionsVcfToBed(cli::AppPtr app, const VcfToBedParamPtr& params);
}  // namespace xoos::svc
