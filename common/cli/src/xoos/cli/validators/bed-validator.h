#pragma once

#include <CLI/Validators.hpp>

#include <xoos/types/fs.h>

namespace xoos::cli {

/**
 * @brief Validates that a BED file is readable and holds at least one region.
 *
 * Leading UCSC track/browser and `#` comment lines and blank lines are skipped, and the first data
 * row is validated (numeric non-negative coordinates, start <= end, at least 3 columns). A file that
 * begins with header/comment lines is accepted as long as a valid region follows. Plain, gzip-, and
 * BGZF-compressed BED files are all supported (`.gz` selects decompression).
 *
 * @param path Path to the file to validate.
 * @throws std::runtime_error Issues encountered during validation:
 *   1. File cannot be opened or decompressed
 *   2. The first data row is malformed (bad coordinates, too few columns, non-tabular content)
 *   3. The file contains no valid regions (empty, header-only, or blank-only)
 */
void ValidateBed(const fs::path& path);

/**
 * @brief Constructs a CLI validator for BED files.
 */
struct BedFileValidator : CLI::Validator {
  BedFileValidator();
};

}  // namespace xoos::cli
