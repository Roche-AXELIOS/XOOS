#pragma once

#include <CLI/Validators.hpp>

#include <xoos/types/fs.h>

namespace xoos::cli {

/**
 * @brief Create a custom validator function based on a provided function.
 * @param validate A function that validates input filesystem path.
 * @return A function that takes a string and returns an error message if validation fails.
 */
std::function<std::string(const std::string&)> CreateValidatorFunction(
    const std::function<void(const fs::path&)>& validate);

}  // namespace xoos::cli
