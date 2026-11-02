#pragma once

#include <optional>

#include <xoos/cli/cli.h>
#include <xoos/io/library-type.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "cli/copy-number-caller-options.h"

namespace xoos::cnc {

// Adds all copy-number-caller subcommands to the CLI app.
void AddSubcommands(CLI::App* app, CopyNumberCallerOptionsPtr& options);

/**
 * @brief Resolves the library type in place and reports whether it is SBX-S.
 *
 * When @p library_type is set (from an explicit --library-type) it is used as
 * is; otherwise the type is auto-detected from the BAM and written back to
 * @p library_type so downstream consumers see the resolved value. Emits an info
 * log for an explicit type and a warning when auto-detection is inconclusive
 * (in which case @p library_type is left unset). Shared so subcommands that
 * branch on SBX-S input resolve the library type consistently.
 *
 * @param library_type In/out: explicit --library-type value, or unset to
 *        auto-detect; populated with the auto-detected type on success.
 * @param bam_path BAM to auto-detect from when no explicit type is given.
 * @return true if the resolved library type is SBX-S; false otherwise (including
 *         when the type cannot be determined).
 */
bool IsSbxSLibraryType(std::optional<io::LibraryType>& library_type, const fs::path& bam_path);

void ConfigureSomaticLikelihoodOptions(CopyNumberCallerOptions& options);
void ConfigureGermlineLikelihoodOptions(CopyNumberCallerOptions& options);
void ConfigureAugmentIntervalsOptions(const CLI::App* app, CopyNumberCallerOptions& options);
void ConfigureSomaticTumorNormalWGSOptions(const CLI::App* app, CopyNumberCallerOptions& options);
void ConfigureSomaticTumorOnlyWGSOptions(const CLI::App* app, CopyNumberCallerOptions& options);
void ConfigurePurityPloidySearchOptions(const CLI::App* app, CopyNumberCallerOptions& options);
void ConfigureGermlineNormalWGSOptions(const CLI::App* app, CopyNumberCallerOptions& options);
void ConfigureSegmentationOptions(const CLI::App* app, CopyNumberCallerOptions& options);

/**
 * @brief Applies post-parse fixups: resolves output file paths relative to
 *        output_dir, resets conditional outputs (BAF, BW), and derives sample_id.
 *        Called automatically by CopyNumberCallerCliMain; expose here for tests.
 */
void CopyNumberCallerCliPreCallback(const CLI::App* app, CopyNumberCallerOptions& options);

s32 CopyNumberCallerCliMain(s32 argc, char** argv);
}  // namespace xoos::cnc
