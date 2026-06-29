#pragma once

#include <CLI/Validators.hpp>

#include <xoos/types/fs.h>

namespace xoos::cli {

/**
 * @brief Validates that a BAM file is readable, correctly formatted, non-empty, and has a valid index.
 * File format is checked using htslib.
 * @param path Path to the BAM file to validate.
 * @throws std::runtime_error Issues encountered during validation:
 *   1. File not readable by htslib
 *   2. Not a BAM file as determined by htslib
 *   3. BAM file is not bgzip-compressed
 *   4. BAM file may be truncated (EOF marker not found)
 *   5. BAM header cannot be read
 *   6. Cannot allocate memory for reading a BAM record
 *   7. BAM file has no records
 *   8. Error reading the BAM file
 *   9. Cannot find or open the corresponding index file
 */
void ValidateBamAndIndex(const fs::path& path);

/**
 * @brief Constructs a CLI validator for indexed BAM files.
 */
struct IndexedBamFileValidator : CLI::Validator {
  IndexedBamFileValidator();
};

}  // namespace xoos::cli
