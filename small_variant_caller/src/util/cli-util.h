#pragma once

#include <xoos/cli/cli.h>

#include "xoos/types/int.h"

namespace xoos::svc {

// Only one subcommand is allowed
constexpr s32 kMinSubcommands = 1;
constexpr s32 kMaxSubcommands = 1;

// CLI option default values shared by two or more applications
constexpr u32 kDefaultMaxBamRegionSizePerThread = 16384;
constexpr u32 kDefaultLeftPad = 0;
constexpr u32 kDefaultRightPad = 0;
constexpr u32 kDefaultCollapseDistance = 0;

// inclusive range for MAPQ in SAM spec
static const CLI::Range kCliRangeMapq(0, 255);

// inclusive range for QUAL in SAM spec
static const CLI::Range kCliRangeBaseq(0, 93);

// inclusive range for fractional values
static const CLI::Range kCliRangeFraction(0.0, 1.0);

/**
 * @brief Adds a command line option to treat warnings as errors.
 * @param app Pointer to the CLI application instance.
 * @return Pointer to the added CLI option.
 */
CLI::Option* AddWarnAsErrorOption(CLI::App* app);

// Each function below adds a series of validation checks to the given CLI option.
// Combined validators (i.e. via `&`) produce an error chaining multiple messages together with "AND", which is not
// user-friendly.
// By adding validators separately, we can ensure that only the first validation error is thrown, making it clearer to
// the user what went wrong.

/**
 * @brief Checks if the file exists and is non-empty.
 * @param opt Pointer to the CLI option to validate.
 * @return Pointer to the validated CLI option.
 */
CLI::Option* CheckNonEmptyFile(CLI::Option* opt);

/**
 * @brief Checks if the file exists, is non-empty, and is a valid indexed BAM file.
 * @param opt Pointer to the CLI option to validate.
 * @return Pointer to the validated CLI option.
 */
CLI::Option* CheckIndexedBamFile(CLI::Option* opt);

/**
 * @brief Checks if the file exists, is non-empty, and is a valid indexed VCF file.
 * @param opt Pointer to the CLI option to validate.
 * @return Pointer to the validated CLI option.
 */
CLI::Option* CheckIndexedVcfFile(CLI::Option* opt);

/**
 * @brief Checks if the file exists, is non-empty, and is a valid indexed FASTA file.
 * @param opt Pointer to the CLI option to validate.
 * @return Pointer to the validated CLI option.
 */
CLI::Option* CheckIndexedFastaFile(CLI::Option* opt);

/**
 * @brief Checks if the file exists, is non-empty, and is a valid BED file.
 * @param opt Pointer to the CLI option to validate.
 * @return Pointer to the validated CLI option.
 */
CLI::Option* CheckBedFile(CLI::Option* opt);

}  // namespace xoos::svc
