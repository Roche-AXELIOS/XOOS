#pragma once

#include <string>

namespace xoos::alignment_metrics {

// Default file names for metrics output files
// accuracy metrics
const std::string kDefaultErrorByClusterSize = "errors_by_cluster_size.tsv";
const std::string kDefaultErrorBySubstitutionType = "errors_by_substitution_type.tsv";
const std::string kDefaultErrorByReadType = "errors_by_read_type.tsv";
const std::string kDefaultBaseLevelAccuracySummary = "base_level_accuracy_summary.tsv";
const std::string kDefaultQscoreStats = "qscore_stats.tsv";

// coverage metrics
const std::string kDefaultCoverageHistograms = "coverage_histograms.tsv";
const std::string kDefaultHpCoverageHistograms = "hp_coverage_histograms.tsv";
const std::string kDefaultHpCoverageStats = "hp_coverage_stats.tsv";
const std::string kDefaultHpCoverageDistributionSummary = "hp_coverage_distribution_summary.tsv";
const std::string kDefaultCoverageStats = "coverage_stats.tsv";
const std::string kDefaultCoverageDistributionSummary = "coverage_distribution_summary.tsv";
// coverage uniformity metrics
const std::string kDefaultCoverageUniformitySummary = "coverage_uniformity_summary.tsv";
const std::string kDefaultMeanCoverageHistogram = "mean_coverage_histogram.tsv";

// read metrics
const std::string kDefaultReadCountsByClusterSize = "read_counts_by_cluster_size.tsv";
const std::string kDefaultReadMetricsSummary = "read_metrics_summary.tsv";
const std::string kDefaultReadLengthHistograms = "read_length_histograms.tsv";
const std::string kDefaultReadLengthSummary = "read_length_summary.tsv";
const std::string kDefaultTeReadMetrics = "read_metrics_te.tsv";

// hp error output
const std::string kDefaultHpErrors = "hp_errors.tsv";

// metrics directory names
const std::string kDefaultAccuracyDir = "accuracy";
const std::string kDefaultCoverageDir = "coverage";
const std::string kDefaultReadMetricsDir = "read_metrics";

}  // namespace xoos::alignment_metrics
