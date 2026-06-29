#pragma once

#include <CLI/Validators.hpp>

#include <xoos/types/fs.h>

namespace xoos::cli {

/**
 * @brief Validates that a BED file is readable.
 * @param path Path to the file to validate.
 * @throws std::runtime_error Issues encountered during validation:
 *   1. File not readable by htslib
 *   2. Not a BED file as determined by htslib
 */
void ValidateBed(const fs::path& path);

/**
 * @brief Constructs a CLI validator for BED files.
 */
struct BedFileValidator : CLI::Validator {
  BedFileValidator();
};

}  // namespace xoos::cli
