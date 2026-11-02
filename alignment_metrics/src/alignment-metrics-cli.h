#pragma once
#include <memory>
#include <string>

#include <xoos/cli/cli.h>

#include "alignment-metrics-options.h"

namespace xoos::alignment_metrics {

const std::string kOptGroupNameInputOptions = "Input Options";
const std::string kOptGroupNameOutputOptions = "Output Options";
const std::string kOptGroupNameRegionOptions = "Region Options";
const std::string kOptGroupNameReadFilteringOptions = "Read Filtering Options";
const std::string kOptGroupNameBaseFilteringOptions = "Base Filtering Options";
const std::string kOptGroupNameHomopolymerMetricOptions = "Homopolymer Metric Options";
const std::string kOptGroupNameHpMaskingOptions = "Discordant Homopolymer Masking Options";
const std::string kOptGroupNameAccuracyMetricsOptions = "Accuracy Metrics Options";
const std::string kOptGroupNameCoverageMetricsOptions = "Coverage Metrics Options";
const std::string kOptGroupNameGcBiasMetricsOptions = "GC Bias Metrics Options";
const std::string kOptGroupNameReadMetricsOptions = "Read Metrics Options";
const std::string kOptGroupNameTargetEnrichmentMetricsOptions = "Target Enrichment Metrics Options";
const std::string kOptGroupNameGeneralMetricsOptions = "General Metrics Options";
const std::string kOptGroupNamePerformanceOptions = "Performance Options";

const std::string kSubcommandNameCoverageMetrics = "coverage";
const std::string kSubcommandNameAccuracyMetrics = "accuracy";
const std::string kSubcommandNameReadMetrics = "read";
const std::string kSubcommandNameGcBiasMetrics = "gcbias";
const std::string kSubcommandNameAlignmentMetrics = "all";

const std::string kSubcommandDescriptionCoverageMetrics = "Calculate coverage metrics";
const std::string kSubcommandDescriptionAccuracyMetrics = "Calculate accuracy metrics";
const std::string kSubcommandDescriptionReadMetrics = "Calculate read metrics";
const std::string kSubcommandDescriptionGcBiasMetrics = "Calculate GC bias metrics";
const std::string kSubcommandDescriptionAlignmentMetrics = "Calculate all metrics";

using AlignmentMetricsOptsPtr = std::shared_ptr<AlignmentMetricsOptions>;

// Define the CLI options necessary for coverage metrics.
void DefineCoverageMetricsOptions(CLI::App* app, const AlignmentMetricsOptsPtr& opts);
// Define the CLI options necessary for accuracy metrics.
void DefineAccuracyMetricsOptions(CLI::App* app, const AlignmentMetricsOptsPtr& opts);
// Define the CLI options necessary for read metrics.
void DefineReadMetricsOptions(CLI::App* app, const AlignmentMetricsOptsPtr& opts);
// Define the CLI options necessary for GC bias metrics.
void DefineGcBiasMetricsOptions(CLI::App* app, const AlignmentMetricsOptsPtr& opts);
// Define the CLI options necessary for all metrics.
void DefineAllMetricsOptions(CLI::App* app, const AlignmentMetricsOptsPtr& opts);

void AddInputOptions(CLI::App* app, const AlignmentMetricsOptsPtr& opts);
void AddOutputOptions(CLI::App* app, const AlignmentMetricsOptsPtr& opts);
void AddRegionOptions(CLI::App* app, const AlignmentMetricsOptsPtr& opts);
void AddReadFilteringOptions(CLI::App* app, const AlignmentMetricsOptsPtr& opts);
void AddBaseFilteringOptions(CLI::App* app, const AlignmentMetricsOptsPtr& opts);
void AddHomopolymerMetricOptions(CLI::App* app, const AlignmentMetricsOptsPtr& opts);
void AddHpMaskingOption(CLI::App* app, const AlignmentMetricsOptsPtr& opts);
void AddAccuracyMetricsOptions(CLI::App* app, const AlignmentMetricsOptsPtr& opts);
void AddCoverageMetricsOptions(CLI::App* app, const AlignmentMetricsOptsPtr& opts);
void AddGcBiasOptions(CLI::App* app, const AlignmentMetricsOptsPtr& opts);
void AddTargetEnrichmentMetricsOptions(CLI::App* app, const AlignmentMetricsOptsPtr& opts);
void AddReadMetricsOptions(CLI::App* app, const AlignmentMetricsOptsPtr& opts);
void AddGeneralMetricsOptions(CLI::App* app, const AlignmentMetricsOptsPtr& opts);
void AddPerformanceOptions(CLI::App* app, const AlignmentMetricsOptsPtr& opts, bool ref_required);
void AddHiddenOptions(CLI::App* app, const AlignmentMetricsOptsPtr& opts, bool include_hp_subsampling_options);

// Run Application for CLI with subcommands for coverage metrics.
int RunAlignmentMetricsApp(int argc, char** argv, const std::string& program_name, const std::string& version);
}  // namespace xoos::alignment_metrics
