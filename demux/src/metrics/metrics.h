#pragma once

namespace xoos::demux {

// Common metric name constants for both Duplex and Simplex workflows

// ============================================================================
// Run Metrics - Common
// ============================================================================
constexpr auto kNumExpectedSids = "num_expected_sids";
constexpr auto kNumSids = "num_sids";
constexpr auto kTotalReads = "total_reads";
constexpr auto kAssignedReads = "assigned_reads";
constexpr auto kPassingReads = "passing_reads";
constexpr auto kFailedAssignedReads = "failed_assigned_reads";
constexpr auto kTooShortTrimmedReads = "too_short_trimmed_reads";
constexpr auto kTooShortReads = "too_short_reads";
constexpr auto kTooManyErrorsReads = "too_many_errors_reads";
constexpr auto kUnassignedReads = "unassigned_reads";

// ============================================================================
// Sample Metrics - Common
// ============================================================================
constexpr auto kMetric = "metric";
constexpr auto kIndexSequence = "index_sequence";
constexpr auto kMetricsName = "metrics_name";
constexpr auto kCount = "count";
constexpr auto kPercentage = "percentage";

}  // namespace xoos::demux
