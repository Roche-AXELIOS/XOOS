#include "core/read-collapser-cli.h"

#include <xoos/cli/cli.h>
#include <xoos/cli/enum-option-util.h>
#include <xoos/cli/validators/file-extension-validator.h>
#include <xoos/enum/enum-util.h>
#include <xoos/types/int.h>
#include <xoos/util/container-functions.h>
#include <xoos/util/file-functions.h>

#include "consensus/consensus.h"
#include "core/read-collapser-options.h"
#include "mark-duplicate/mark-duplicate.h"
#include "util/cli-option-util.h"

namespace xoos::read_collapser {

void SetCommandLineInfo(const CLI::App* const app, const ReadCollapserOptionsPtr& options) {
  options->command_line_info = cli::GetCommandLineInfo(app);
}

void DefineConsensusOptions(CLI::App* const app, const ReadCollapserOptionsPtr& options) {
  // Options are added in the standard help-menu section order:
  // Preset -> Input -> Output -> Region -> Read Filtering -> Cluster -> Consensus -> Performance.
  AddPresetOption(app, kConsensusPresetsMap, kOptGroupNamePresetOptions);
  AddCommonInputOptions(app, options, kOptGroupNameInputOptions);
  AddCommonOutputOptions(app, options, kOptGroupNameOutputOptions);
  AddConsensusOutputOptions(app, options, kOptGroupNameOutputOptions);
  AddCommonRegionOptions(app, options, kOptGroupNameRegionOptions);
  AddCommonReadFilteringOptions(app, options, kOptGroupNameReadFilterOptions);
  AddConsensusReadFilterOptions(app, options, kOptGroupNameReadFilterOptions);
  AddCommonClusterOptions(app, options, kOptGroupNameClusterOptions);
  AddConsensusClusterOptions(app, options, kOptGroupNameClusterOptions);
  AddConsensusOptions(app, options, kOptGroupNameConsensusOptions);
  AddCommonPerformanceOptions(app, options, kOptGroupNamePerformanceOptions);
  AddCommonHiddenOptions(app, options);
}

void DefineMarkdupOptions(CLI::App* const app, const ReadCollapserOptionsPtr& options) {
  // Options are added in the standard help-menu section order:
  // Preset -> Input -> Output -> Region -> Read Filtering -> Cluster -> Performance.
  AddPresetOption(app, kMarkdupPresetsMap, kOptGroupNamePresetOptions);
  AddCommonInputOptions(app, options, kOptGroupNameInputOptions);
  AddCommonOutputOptions(app, options, kOptGroupNameOutputOptions);
  AddMarkdupOutputOptions(app, options, kOptGroupNameOutputOptions);
  AddCommonRegionOptions(app, options, kOptGroupNameRegionOptions);
  AddCommonReadFilteringOptions(app, options, kOptGroupNameReadFilterOptions);
  AddMarkdupReadFilterOptions(app, options, kOptGroupNameReadFilterOptions);
  AddCommonClusterOptions(app, options, kOptGroupNameClusterOptions);
  AddCommonPerformanceOptions(app, options, kOptGroupNamePerformanceOptions);
  AddCommonHiddenOptions(app, options);
}

void DefineOptions(CLI::App* const app, ReadCollapserOptionsPtr& options) {
  // Add 'markdup' subcommand
  cli::AddSubcommand<ReadCollapserOptions>(app,
                                           "markdup",
                                           DefineMarkdupOptions,
                                           options,
                                           MarkDuplicateAndMergeOutput,
                                           "Identifies, marks, and optionally removes duplicate reads.",
                                           [](const CLI::App* const app, const ReadCollapserOptionsPtr& options) {
                                             ValidateMarkdupOptions(app, options);
                                             SetCommandLineInfo(app, options);
                                           });

  // Add 'consensus' subcommand
  cli::AddSubcommand<ReadCollapserOptions>(app,
                                           "consensus",
                                           DefineConsensusOptions,
                                           options,
                                           FastClusterAndConsensus,
                                           "Generates high-quality consensus sequences.",
                                           [](const CLI::App* const app, const ReadCollapserOptionsPtr& options) {
                                             ValidateConsensusOptions(app, options);
                                             SetCommandLineInfo(app, options);
                                           });

  // Require exactly one subcommand to be executed
  app->require_subcommand(1);
}

}  // namespace xoos::read_collapser
