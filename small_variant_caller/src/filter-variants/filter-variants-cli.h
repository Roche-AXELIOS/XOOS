#pragma once

#include <memory>

#include <xoos/cli/cli.h>

#include "filter-variants.h"

namespace xoos::svc {
using FilterVariantsParamPtr = std::shared_ptr<FilterVariantsParam>;

void DefineOptionsFilterVariants(cli::AppPtr app, const FilterVariantsParamPtr& params);
void FilterVariantsPreCallback(cli::ConstAppPtr app, const FilterVariantsParamPtr& params);

}  // namespace xoos::svc
