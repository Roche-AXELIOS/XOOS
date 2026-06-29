#pragma once

#include <CLI/Validators.hpp>

#include <xoos/types/fs.h>

namespace xoos::cli {

/**
 * @brief Validates that a file has a size larger than zero.
 * @param path Path to the file to validate.
 * @throws std::runtime_error File is empty.
 */
void ValidateNonEmptyFile(const fs::path& path);

/**
 * @brief Constructs a CLI validator for non-empty files.
 */
struct NonEmptyFileValidator : CLI::Validator {
  NonEmptyFileValidator();
};

}  // namespace xoos::cli
