#pragma once
#include <memory>
#include <string>

#include <xoos/cli/cli.h>

#include "alignment-metrics-options.h"

namespace xoos::alignment_metrics {

const std::string kOptGroupNameInputOptions = "Input Options";
const std::string kOptGroupNameOutputOptions = "Output Options";
const std::string kOptGroupNameReadFilterOptions = "Read-level Filter Options";
const std::string kOptGroupNameBaseFilterOptions = "Base-level Filter Options";
const std::string kOptGroupNameReadTrimmingOptions = "Read Trimming Options";
const std::string kOptGroupNameHpMetricsOptions = "Homopolymer Metrics Options";
const std::string kOptGroupNameHpMaskingOptions = "Discordant Homopolymer Masking Options";
const std::string kOptGroupNameAccuracyMetricsOptions = "Accuracy Metrics Options";
const std::string kOptGroupNameCoverageMetricsOptions = "Coverage Metrics Options";
const std::string kOptGroupNameReadMetricsOptions = "Read Metrics Options";
const std::string kOptGroupNameSummaryStatsOptions = "Summary Statistics Options";
const std::string kOptGroupNamePerformanceOptions = "Performance Options";
const std::string kOptGroupNameTEMetricsOptions = "Target Enrichment Metrics Options";

const std::string kSubcommandNameCoverageMetrics = "coverage";
const std::string kSubcommandNameAccuracyMetrics = "accuracy";
const std::string kSubcommandNameReadMetrics = "read";
const std::string kSubcommandNameAlignmentMetrics = "all";

const std::string kSubcommandDescriptionCoverageMetrics = "Calculate coverage metrics";
const std::string kSubcommandDescriptionAccuracyMetrics = "Calculate accuracy metrics";
const std::string kSubcommandDescriptionReadMetrics = "Calculate read metrics";
const std::string kSubcommandDescriptionAlignmentMetrics = "Calculate all metrics";

using CoverageMetricsOptsPtr = std::shared_ptr<AlignmentMetricsOptions>;

// Define the CLI options necessary for coverage metrics.
void DefineCoverageMetricsOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts);
// Define the CLI options necessary for accuracy metrics.
void DefineAccuracyMetricsOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts);
// Define the CLI options necessary for read metrics.
void DefineReadMetricsOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts);
// Define the CLI options necessary for all metrics.
void DefineAllMetricsOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts);

void AddInputOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts, bool requires_reference);
void AddOutputOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts);
void AddReadFilterOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts);
void AddBaseFilterOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts);
void AddHpMetricsOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts);
void AddReadTrimmingOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts);
void AddHpMaskingOption(cli::AppPtr app, const CoverageMetricsOptsPtr& opts);
void AddAccuracyMetricsOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts);
void AddReadMetricsOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts);
void AddPerformanceOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts, bool ref_required);
void AddTEMetricsOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts);
void AddLegacyOptions(cli::AppPtr app, const CoverageMetricsOptsPtr& opts);

// Run Application for CLI with subcommands for coverage metrics.
int RunAlignmentMetricsApp(int argc, char** argv, const std::string& program_name, const std::string& version);
}  // namespace xoos::alignment_metrics
