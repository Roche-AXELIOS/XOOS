#pragma once

#include <string>

#include <CLI/App.hpp>

#include <xoos/cli/cli.h>

#include "cli/copy-number-caller-options.h"  // IWYU pragma: keep

namespace xoos::cnc {

// End-to-end pipeline options
void DefineSomaticTumorNormalWGSOptions(CLI::App* app, CopyNumberCallerOptionsPtr& options);
void DefineSomaticTumorOnlyWGSOptions(CLI::App* app, const CopyNumberCallerOptionsPtr& options);
void DefineGermlineNormalWGSOptions(CLI::App* app, CopyNumberCallerOptionsPtr& options);

// PredictSomaticCNA options
void DefinePredictSomaticCNASubcommandOptions(CLI::App* app, CopyNumberCallerOptionsPtr& options);

// Shared options
void DefineSomaticCommonOptions(CLI::App* app, CopyNumberCallerOptions& options, SegmentUndoMode undo_mode);
void DefineSomaticSegmentationOptions(CLI::App* app, SegmentationOptions& options, bool hide);
void DefineSegmentationAlgorithmOptions(CLI::App* app,
                                        SegmentationOptions& options,
                                        bool hide,
                                        SegmentUndoMode undo_mode);
void DefineBafFilteringOptions(CLI::App* app, CopyNumberCallerOptions& options, bool hide);
void DefineGermlineCallFlagOptions(CLI::App* app, LikelihoodOptions& options, bool hide);
void DefineAugmentIntervalsInputOptions(CLI::App* app,
                                        CopyNumberCallerOptions& options,
                                        bool hide,
                                        bool require_mappability);
void DefineIntervalOptions(
    CLI::App* app, AugmentIntervalsOptions& options, bool hide, size_t default_interval_size, f64 default_mappability);
void DefineIntervalOptions(CLI::App* app,
                           AugmentIntervalsOptions& options,
                           bool hide,
                           size_t default_interval_size,
                           f64 default_mappability,
                           const std::string& interval_size_description_suffix);
void DefineAdvancedIntervalOptions(CLI::App* app, AugmentIntervalsOptions& options);
void DefineCalculateCoverageAlgorithmOptions(CLI::App* app, CalculateCoverageOptions& options, bool hide);
void DefineGCCorrectAlgorithmOptions(CLI::App* app, GCCorrectOptions& options, bool hide);
void DefineCoverageOptions(CLI::App* app, CopyNumberCallerOptions& options);
void DefinePurityPloidySearchAlgorithmOptions(CLI::App* app, PurityPloidySearchOptions& options, bool hide);
void DefinePurityPloidySearchAlgorithmOptions(CLI::App* app,
                                              PurityPloidySearchOptions& options,
                                              bool hide,
                                              f64 baf_likelihood_weight_default);

// Hidden subcommand options
void DefineSegmentationSubcommandOptions(CLI::App* app, CopyNumberCallerOptionsPtr& options);
void DefinePredictGermlineCNVSubcommandOptions(CLI::App* app, CopyNumberCallerOptionsPtr& options);
void DefineAugmentIntervalsSubcommandOptions(CLI::App* app, CopyNumberCallerOptionsPtr& options);
void DefineCalculateCoverageSubcommandOptions(CLI::App* app, CopyNumberCallerOptionsPtr& options);
void DefineGCCorrectSubcommandOptions(CLI::App* app, CopyNumberCallerOptionsPtr& options);
void DefineOneSampleLogRSubcommandOptions(CLI::App* app, CopyNumberCallerOptionsPtr& options);
void DefineTwoSampleLogRSubcommandOptions(CLI::App* app, CopyNumberCallerOptionsPtr& options);
void DefinePurityPloidySearchSubcommandOptions(CLI::App* app, CopyNumberCallerOptionsPtr& options);
void DefineMergeSegmentsSubcommandOptions(CLI::App* app, CopyNumberCallerOptionsPtr& options);
void DefineSegToVcfSubcommandOptions(CLI::App* app, CopyNumberCallerOptionsPtr& options);
void DefinePredictSexSubcommandOptions(CLI::App* app, CopyNumberCallerOptionsPtr& options);
void DefineTumorOnlyBafsSubcommandOptions(CLI::App* app, CopyNumberCallerOptionsPtr& options);

// =============================================================================
// Add*Option helper functions
// =============================================================================

// Shared End-to-End options
CLI::Option* AddIntervalsOption(CLI::App* app, const std::string& group_name, std::optional<fs::path>& binding);

// Shared IO options
CLI::Option* AddSegmentsOption(CLI::App* app, const std::string& group_name, std::optional<fs::path>& binding);
CLI::Option* AddLogrsOption(CLI::App* app, const std::string& group_name, std::optional<fs::path>& binding);
CLI::Option* AddBafsOption(CLI::App* app, const std::string& group_name, std::optional<fs::path>& binding);
CLI::Option* AddMapqsOption(CLI::App* app, const std::string& group_name, std::optional<fs::path>& binding);
CLI::Option* AddVcfOption(CLI::App* app,
                          const std::string& group_name,
                          std::optional<fs::path>& binding,
                          const std::string& description);
CLI::Option* AddVcfOption(CLI::App* app,
                          const std::string& group_name,
                          std::optional<fs::path>& binding,
                          const std::string& description,
                          const std::string& option_name);
CLI::Option* AddOutputDirOption(CLI::App* app, const std::string& group_name, fs::path& binding);
CLI::Option* AddOutputPrefixOption(CLI::App* app, const std::string& group_name, std::string& binding);
CLI::Option* AddOverwriteOption(CLI::App* app, const std::string& group_name, bool& binding);

// VCF Parsing
CLI::Option* AddNormalSampleMinDepthOption(CLI::App* app,
                                           const std::string& group_name,
                                           s32& binding,
                                           s32 default_depth);
CLI::Option* AddNormalSampleMinBafOption(CLI::App* app, const std::string& group_name, f64& binding);
CLI::Option* AddNormalSampleMaxBafOption(CLI::App* app, const std::string& group_name, f64& binding);
CLI::Option* AddTumorSampleMinDepthOption(CLI::App* app, const std::string& group_name, s32& binding);
CLI::Option* AddTumorSampleMinBafOption(CLI::App* app, const std::string& group_name, f64& binding);
CLI::Option* AddTumorSampleMaxBafOption(CLI::App* app, const std::string& group_name, f64& binding);
CLI::Option* AddAdInOption(CLI::App* app, const std::string& group_name, std::optional<fs::path>& binding);
CLI::Option* AddForceEnableSomaticVariantParsingOption(CLI::App* app, const std::string& group_name, bool& binding);

// AugmentIntervals Inputs
CLI::Option* AddMappabilityBigwigOption(CLI::App* app, const std::string& group_name, fs::path& binding);
CLI::Option* AddBlocklistBedOption(CLI::App* app, const std::string& group_name, std::optional<fs::path>& binding);

// AugmentIntervals Other
CLI::Option* AddIncludeAltContigsOption(CLI::App* app, const std::string& group_name, bool& binding);

// AugmentIntervals Whole Genome
CLI::Option* AddWholeGenomeOption(CLI::App* app, const std::string& group_name, bool& binding);
CLI::Option* AddWholeGenomeIntervalSizeOption(CLI::App* app,
                                              const std::string& group_name,
                                              size_t& binding,
                                              size_t default_size);
CLI::Option* AddWholeGenomeIntervalSizeOption(CLI::App* app,
                                              const std::string& group_name,
                                              size_t& binding,
                                              size_t default_size,
                                              const std::string& description_suffix);
CLI::Option* AddMinGcContentOption(CLI::App* app, const std::string& group_name, f64& binding);
CLI::Option* AddMaxGcContentOption(CLI::App* app, const std::string& group_name, f64& binding);
CLI::Option* AddWholeGenomeMinMappabilityOption(CLI::App* app,
                                                const std::string& group_name,
                                                f64& binding,
                                                f64 default_mappability);

// AugmentIntervals Target Options
CLI::Option* AddBedOption(CLI::App* app, const std::string& group_name, std::optional<fs::path>& binding);
CLI::Option* AddNoOffTargetsOption(CLI::App* app, const std::string& group_name, bool& binding);
CLI::Option* AddOnTargetIntervalSizeOption(CLI::App* app, const std::string& group_name, size_t& binding);
CLI::Option* AddOffTargetMinWidthOption(CLI::App* app, const std::string& group_name, size_t& binding);
CLI::Option* AddOffTargetTrimLengthOption(CLI::App* app, const std::string& group_name, size_t& binding);
CLI::Option* AddOffTargetIntervalSizeOption(CLI::App* app, const std::string& group_name, size_t& binding);
CLI::Option* AddOnTargetMinMappabilityOption(CLI::App* app, const std::string& group_name, f64& binding);
CLI::Option* AddOffTargetMinMappabilityOption(CLI::App* app, const std::string& group_name, f64& binding);
CLI::Option* AddSexChromosomeTargetMinMappabilityOption(CLI::App* app, const std::string& group_name, f64& binding);

// CalculateCoverage
CLI::Option* AddCoverageExcludeFlagsOption(CLI::App* app, const std::string& group_name, std::string& binding);
CLI::Option* AddCoverageIgnoreDnOption(CLI::App* app, const std::string& group_name, bool& binding);

// GCCorrect
CLI::Option* AddGcCorrectFirstSpanOption(CLI::App* app, const std::string& group_name, f64& binding);

// TwoSampleLogR / Shared coverage
CLI::Option* AddTumorMinCoverageOption(CLI::App* app, const std::string& group_name, size_t& binding);
CLI::Option* AddNormalMinCoverageOption(CLI::App* app, const std::string& group_name, size_t& binding);
CLI::Option* AddTumorCoverageFileOption(CLI::App* app, const std::string& group_name, std::optional<fs::path>& binding);
CLI::Option* AddNormalCoverageFileOption(CLI::App* app,
                                         const std::string& group_name,
                                         std::optional<fs::path>& binding);
CLI::Option* AddCoverageFileOption(CLI::App* app, const std::string& group_name, std::optional<fs::path>& binding);

// Segmentation
CLI::Option* AddSeedSegmentsOption(CLI::App* app, const std::string& group_name, fs::path& binding);
CLI::Option* AddDisableHierarchicalPruningOption(CLI::App* app, const std::string& group_name, bool& binding);
CLI::Option* AddDisableMergingOption(CLI::App* app, const std::string& group_name, bool& binding);
CLI::Option* AddCbsMethodOption(CLI::App* app, const std::string& group_name, CbsMaxTMethod& binding);
CLI::Option* AddSegmentationModeOption(CLI::App* app, const std::string& group_name, SegmentationMode& binding);
CLI::Option* AddTruncateOutliersOption(CLI::App* app, const std::string& group_name, bool& binding);
CLI::Option* AddPruningClusteringParameterOption(CLI::App* app, const std::string& group_name, f64& binding);
CLI::Option* AddBafDistanceWeightOption(CLI::App* app, const std::string& group_name, f64& binding);
CLI::Option* AddMinTForAutomaticSegmentationOption(CLI::App* app, const std::string& group_name, f64& binding);
CLI::Option* AddUndoLogrSegmentsOption(CLI::App* app,
                                       const std::string& group_name,
                                       SegmentUndoMode& binding,
                                       SegmentUndoMode default_value);
CLI::Option* AddUndoLogrSegmentsSdFactorOption(CLI::App* app, const std::string& group_name, f64& binding);
CLI::Option* AddIncrementUndoLogrSegmentsSdFactorOption(CLI::App* app, const std::string& group_name, f64& binding);
CLI::Option* AddMaxNumSegmentsOption(CLI::App* app, const std::string& group_name, size_t& binding);
CLI::Option* AddUndoDhSegmentsOption(CLI::App* app, const std::string& group_name, bool& binding);
CLI::Option* AddUndoDhSegmentsSdFactorOption(CLI::App* app, const std::string& group_name, f64& binding);
CLI::Option* AddMinObservationsPerSegmentOption(CLI::App* app, const std::string& group_name, size_t& binding);

// PurityPloidySearch
CLI::Option* AddSegMaxLogrOption(CLI::App* app, const std::string& group_name, f64& binding);
CLI::Option* AddSegMinNumLogrsOption(CLI::App* app, const std::string& group_name, size_t& binding);
CLI::Option* AddSegMinNumSnpsOption(CLI::App* app, const std::string& group_name, size_t& binding);
CLI::Option* AddBafLikelihoodWeightForPurityPloidySearchRerunOption(CLI::App* app,
                                                                    const std::string& group_name,
                                                                    f64& binding,
                                                                    f64 default_value);
CLI::Option* AddMinLohFractionForWgdOption(CLI::App* app, const std::string& group_name, f64& binding);
CLI::Option* AddTumorOnlyOption(CLI::App* app,
                                const std::string& group_name,
                                bool& binding,
                                const std::string& description);

// Sample Metadata
CLI::Option* AddSexOption(CLI::App* app, const std::string& group_name, std::optional<sex_predict::Sex>& binding);
CLI::Option* AddSexSelectionOption(CLI::App* app, const std::string& group_name, sex_predict::SexSelection& binding);

CLI::Option* AddLibraryTypeOption(CLI::App* app,
                                  const std::string& group_name,
                                  std::optional<io::LibraryType>& binding);
CLI::Option* AddSampleNameOption(CLI::App* app, const std::string& group_name, std::optional<std::string>& binding);
CLI::Option* AddTumorPurityOption(CLI::App* app, const std::string& group_name, std::optional<f64>& binding);
CLI::Option* AddTumorPloidyOption(CLI::App* app, const std::string& group_name, std::optional<f64>& binding);
CLI::Option* AddSampleIdOption(CLI::App* app, const std::string& group_name, std::string& binding);
CLI::Option* AddNormalSampleNameOption(CLI::App* app,
                                       const std::string& group_name,
                                       std::optional<std::string>& binding);
CLI::Option* AddTumorSampleNameOption(CLI::App* app,
                                      const std::string& group_name,
                                      std::optional<std::string>& binding);

// Likelihood
CLI::Option* AddCnvLengthFlagMinSizeOption(CLI::App* app, const std::string& group_name, size_t& binding);
CLI::Option* AddMapqCutoffToFlagCallsOption(CLI::App* app, const std::string& group_name, s32& binding);

// Save-flag helpers (new hidden advanced options)
CLI::Option* AddSaveTaskflowGraphOption(CLI::App* app, const std::string& group_name, bool& binding);
CLI::Option* AddSaveLogRatioSegmentsOption(CLI::App* app, const std::string& group_name, bool& binding);
CLI::Option* AddSaveBafSegmentsOption(CLI::App* app, const std::string& group_name, bool& binding);
CLI::Option* AddWriteSegmentLikelihoodOption(CLI::App* app, const std::string& group_name, bool& binding);

// MergeSegments
CLI::Option* AddMinMapqOption(CLI::App* app, const std::string& group_name, f64& binding);
CLI::Option* AddUseMapqsFromFileOption(CLI::App* app, const std::string& group_name, bool& binding);
CLI::Option* AddRecalculatePerSegmentDataOption(CLI::App* app, const std::string& group_name, bool& binding);

}  // namespace xoos::cnc
