#pragma once

#include <memory>

#include <xoos/cli/cli.h>

#include "train-model.h"

namespace xoos::svc {

using TrainModelParamPtr = std::shared_ptr<TrainModelParam>;

void DefineOptionsTrainModel(cli::AppPtr app, const TrainModelParamPtr& params);
void TrainModelPreCallback(cli::ConstAppPtr app, const TrainModelParamPtr& params);

}  // namespace xoos::svc
