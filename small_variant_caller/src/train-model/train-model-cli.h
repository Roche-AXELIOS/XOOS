#pragma once

#include <memory>

#include <xoos/cli/cli.h>

#include "train-model.h"

namespace xoos::svc {
using TrainModelParamPtr = std::shared_ptr<TrainModelParam>;
}  // namespace xoos::svc

namespace xoos::svc::train_model {
/**
 * @brief Define CLI options for the main application.
 * @param app Main application pointer where CLI options will be defined.
 * @param params Shared pointer to CLI parameters to store parsed option values
 */
void DefineOptions(CLI::App* app, TrainModelParamPtr& params);

/**
 * @brief Preprocess the CLI options and validate the parameters before running the main application.
 * @param app Main application pointer where the callback is triggered
 * @param params Shared pointer to CLI parameters to be validated
 * @throws CLI::ValidationError if the parameters are invalid
 */
void PreCallback(cli::ConstAppPtr app, const TrainModelParamPtr& params);
}  // namespace xoos::svc::train_model
