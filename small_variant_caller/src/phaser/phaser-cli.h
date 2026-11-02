#pragma once

#include <memory>

#include <xoos/cli/cli.h>

#include "params.h"

namespace xoos::svc {
using PhaserParamPtr = std::shared_ptr<phaser::Params>;
}  // namespace xoos::svc

namespace xoos::svc::phaser {
/**
 * @brief Define CLI options for the phase_variants application.
 * @param app Main application pointer where CLI options will be defined.
 * @param params Shared pointer to CLI parameters to store parsed option values.
 */
// Invoked via the .define_options function pointer in apps/phase-variants.cpp; Sonar analyses only the src tree.
void DefineOptions(CLI::App* app, const PhaserParamPtr& params);  // NOSONAR(cpp:S5536)

/**
 * @brief Validate the parsed CLI parameters before running.
 * @param app Main application pointer where the callback is triggered.
 * @param params Shared pointer to CLI parameters to be validated.
 * @throws CLI::ValidationError if the parameters are invalid or the output VCF already exists.
 */
// Registered as a CLI pre-run callback in apps/phase-variants.cpp; Sonar analyses only the src tree.
void PreCallback(const CLI::App* app, const PhaserParamPtr& params);  // NOSONAR(cpp:S5536)
}  // namespace xoos::svc::phaser
