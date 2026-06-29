#pragma once

#include <CLI/Validators.hpp>

#include <xoos/types/fs.h>

namespace xoos::cli {

/**
 * @brief Validates that a VCF file is readable, correctly formatted, non-empty, and has a valid index.
 * @param path Path to the VCF file to validate.
 * @throws std::runtime_error Issue encountered during validation:
 *   1. File not readable by htslib
 *   2. Not a VCF file as determined by htslib
 *   3. VCF file is not bgzip-compressed
 *   4. VCF file may be truncated (EOF marker not found)
 *   5. VCF header cannot be read
 *   6. Cannot allocate memory for reading a VCF record
 *   7. VCF file has no records
 *   8. Error reading the VCF file
 *   9. Cannot find or open the corresponding index file
 */
void ValidateVcfAndIndex(const fs::path& path);

/**
 * @brief Constructs a CLI validator for indexed VCF files.
 */
struct IndexedVcfFileValidator : CLI::Validator {
  IndexedVcfFileValidator();
};

}  // namespace xoos::cli
