#pragma once

#include <string>
#include <string_view>

#include <xoos/cli/cli.h>

#include "core/library-type.h"
#include "xoos/types/int.h"

namespace xoos::svc {

// Only one subcommand is allowed
constexpr s32 kMinSubcommands = 1;
constexpr s32 kMaxSubcommands = 1;

// CLI option default values shared by two or more applications
constexpr u32 kDefaultMaxBamRegionSizePerThread = 16'384;
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

/**
 * @brief Add a `--library-type` option that accepts only concrete SVC library types.
 * @details This helper is for commands without BAM inputs. `auto-detect` depends on BAM header/read-content
 * inspection, so non-BAM commands must require users to choose the concrete library type explicitly.
 * @param app Pointer to the CLI application or subcommand instance.
 * @param library_type Parsed library type destination.
 * @param description User-facing option description.
 * @return Pointer to the added CLI option.
 */
CLI::Option* AddExplicitLibraryTypeOption(CLI::App* app, LibraryType& library_type, const std::string& description);

/**
 * @brief Write a resolved value into an option's rendered default when the user did not provide it explicitly.
 * @details CLI11 stores option defaults when options are registered. Some SVC defaults depend on a configuration
 *          profile resolved later in parsing, so this updates only the rendered default used for logs and metadata.
 * @param sub Subcommand that owns the option.
 * @param option_name CLI option name to update.
 * @param value Resolved value to render in provenance output.
 */
void SetRenderedDefault(CLI::App* sub, std::string_view option_name, const std::string& value);

/**
 * @brief Replace an option value in a rendered command line or append it when absent.
 * @details Matches the flag as a whole token so a shorter flag does not corrupt a longer one that shares its
 * prefix (e.g. `--sample-type` must not match inside `--sample-type-max-af`). Used after deferred option
 * resolution so command-line provenance reflects the resolved configuration.
 * @param command_line Rendered command line updated in place.
 * @param option_name Option flag to replace or append.
 * @param value Rendered option value.
 */
void ReplaceOrAppendCommandLineOption(std::string& command_line,
                                      std::string_view option_name,
                                      const std::string& value);

}  // namespace xoos::svc
