#pragma once

#include <memory>

#include <xoos/cli/cli.h>

#include "generate-pon.h"

namespace xoos::svc {
using GeneratePonParamPtr = std::shared_ptr<GeneratePonParam>;

/**
 * @brief Define CLI options for the generate_pon subcommand.
 * @param app CLI11 application to add options to.
 * @param params Shared pointer to the parameter struct to populate.
 */
void DefineOptionsGeneratePon(CLI::App* app, const GeneratePonParamPtr& params);
}  // namespace xoos::svc

namespace xoos::svc::generate_pon {
/**
 * @brief Capture command line metadata before running the main function.
 * @param app CLI application pointer
 * @param params Shared pointer to CLI parameters
 */
void PreCallback(const CLI::App* app, const GeneratePonParamPtr& params);
}  // namespace xoos::svc::generate_pon
