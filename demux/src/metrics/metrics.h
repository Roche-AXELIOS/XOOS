#pragma once

#include <string>
#include <string_view>

#include <xoos/types/float.h>
#include <xoos/types/int.h>

namespace xoos::demux {

// Common metric name constants for both Duplex and Simplex workflows
// ============================================================================
// Run Metrics - Common
// ============================================================================
constexpr std::string_view kNumExpectedSids = "num_expected_sids";
constexpr std::string_view kNumSids = "num_sids";
constexpr std::string_view kTotalReads = "total_reads";
constexpr std::string_view kAssignedReads = "assigned_reads";
constexpr std::string_view kPassingReads = "passing_reads";
constexpr std::string_view kFailedAssignedReads = "failed_assigned_reads";
constexpr std::string_view kTooShortTrimmedReads = "too_short_trimmed_reads";
constexpr std::string_view kSidDiscordantDiscardedReads = "sid_discordant_discarded_reads";
constexpr std::string_view kFullOnlyDiscardedReads = "full_only_discarded_reads";
constexpr std::string_view kTooShortReads = "too_short_reads";
constexpr std::string_view kTooLongReads = "too_long_reads";
constexpr std::string_view kUnassignedReads = "unassigned_reads";

// ============================================================================
// Base Metrics - Common
// ============================================================================
constexpr std::string_view kTotalBases = "total_bases";
constexpr std::string_view kRawBases = "raw_bases";
constexpr std::string_view kTrimmedBases = "trimmed_bases";
constexpr std::string_view kUnassignedBases = "unassigned_bases";
constexpr std::string_view kFailedAssignedBases = "failed_assigned_bases";
constexpr std::string_view kMeanPassingReadLength = "mean_passing_read_length";

// ============================================================================
// Sample Metrics - Common
// ============================================================================
constexpr std::string_view kMetric = "metric";
constexpr std::string_view kIndexSequence = "index_sequence";
constexpr std::string_view kAssignedReadsPercentage = "assigned_reads_percentage";
constexpr std::string_view kPassingReadsPercentage = "passing_reads_percentage";
constexpr std::string_view kPassingReadsPercentageAllSamples = "passing_reads_percentage_all_samples";
constexpr std::string_view kTotalBasesPercentageAllSamples = "total_bases_percentage_all_samples";
constexpr std::string_view kTrimmedBasesPercentage = "trimmed_bases_percentage";

// Canonical Key-Value Summary header (see docs command-line-interface-best-practices.md, TSV section):
// metric_name / value / percentage / denominator. `denominator` names the metric used as the
// percentage base (not a number); it is `NA` when the percentage is not applicable.
constexpr std::string_view kMetricName = "metric_name";
constexpr std::string_view kValue = "value";
constexpr std::string_view kPercentage = "percentage";
constexpr std::string_view kDenominator = "denominator";

// ============================================================================
// Percentage helpers - shared by Duplex and Simplex writers
// ============================================================================

// Percentage of `numerator` over `denominator` in the [0, 100] range. Returns NaN when
// `denominator` is 0 so the caller can emit the not-applicable sentinel. The division is done in
// f64 to avoid precision loss before narrowing to f32.
f32 CalculatePercentage(u64 numerator, u64 denominator);

// Default number of decimal places for formatted percentages in the metrics files.
constexpr s32 kPercentagePrecision = 2;

// Format the percentage of `numerator` over `denominator` to `precision` decimal places. A zero
// denominator (e.g. 0/0) yields the canonical `NA` sentinel rather than the platform "nan" string.
std::string FormatPercentage(u64 numerator, u64 denominator, s32 precision);

// Format the percentage of `numerator` over `denominator` at the default precision.
std::string FormatPercentage(u64 numerator, u64 denominator);

}  // namespace xoos::demux
