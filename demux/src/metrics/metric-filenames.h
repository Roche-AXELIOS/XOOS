#pragma once

namespace xoos::demux {

// metrics constants
constexpr auto kMetricsDirectory = "metrics";

// metrics file names
constexpr auto kRunMetricsFile = "run_stats.tsv";
constexpr auto kSampleMetricsFile = "sample_stats.tsv";

constexpr auto kPassingReadLengthDistr = "passing_read_length_distr.tsv";
constexpr auto kFullDuplexReadLengthDistr = "full_duplex_read_length_distr.tsv";
constexpr auto kPartialDuplexReadLengthDistr = "partial_duplex_read_length_distr.tsv";
constexpr auto kEndadapterPositionDistr = "endadapter_position_distr.tsv";
constexpr auto kTotalReadLengthDistr = "total_read_length_distr.tsv";
constexpr auto kUnassignedReadLengthDistr = "unassigned_read_length_distr.tsv";
constexpr auto kNoHairpinReadLengthDistr = "no_hairpin_read_length_distr.tsv";

// read length metrics file names
constexpr auto kUntrimmedFullReadLenDist = "untrimmed_full_read_len_dist.tsv";
constexpr auto kUntrimmedPartialReadLenDist = "untrimmed_partial_read_len_dist.tsv";
constexpr auto kTrimmedFullReadLenDist = "trimmed_full_read_len_dist.tsv";
constexpr auto kTrimmedPartialReadLenDist = "trimmed_partial_read_len_dist.tsv";

// metrics file names used in demux-and-trim pipeline
constexpr auto kSampleAssignmentMetrics = "sample_assignment_metrics.tsv";

}  // namespace xoos::demux
