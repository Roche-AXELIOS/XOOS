#pragma once

#include <CLI/Validators.hpp>

#include <xoos/types/fs.h>

namespace xoos::cli {

/**
 * @brief Validates that a FASTA file is readable and has an index.
 * @param path Path to the FASTA file to validate.
 * @throws std::runtime_error Issue encountered during validation:
 *   1. File not readable by htslib
 *   2. Not a FASTA file as determined by htslib
 *   3. Cannot find or open the corresponding index file
 */
void ValidateFastaAndIndex(const fs::path& path);

/**
 * @brief Constructs a CLI validator for indexed FASTA files.
 */
struct IndexedFastaFileValidator : CLI::Validator {
  IndexedFastaFileValidator();
};

}  // namespace xoos::cli
