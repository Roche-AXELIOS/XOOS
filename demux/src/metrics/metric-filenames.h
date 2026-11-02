#pragma once

#include <string>
#include <string_view>

namespace xoos::demux {

// metrics constants
constexpr auto kMetricsDirectory = "metrics";

// metrics file names
constexpr std::string_view kRunMetricsFile = "run_stats.tsv";
constexpr std::string_view kSampleMetricsFile = "sample_stats.tsv";

constexpr std::string_view kTotalEndadapterPositionDistr = "total_endadapter_position_distr.tsv";
constexpr std::string_view kSampleEndadapterPositionDistr = "sample_endadapter_position_distr.tsv";
constexpr std::string_view kTotalReadLengthDistr = "total_read_length_distr.tsv";
constexpr std::string_view kNoHairpinReadLengthDistr = "no_hairpin_read_length_distr.tsv";
constexpr std::string_view kRawFailedReadLengthDistr = "raw_failed_read_length_distr.tsv";

constexpr std::string_view kTotalUntrimmedSimplexPassingReadLengthDistr =
    "total_untrimmed_simplex_passing_read_length_distr.tsv";
constexpr std::string_view kTotalTrimmedSimplexPassingReadLengthDistr =
    "total_trimmed_simplex_passing_read_length_distr.tsv";
constexpr std::string_view kSampleUntrimmedSimplexPassingReadLengthDistr =
    "sample_untrimmed_simplex_passing_read_length_distr.tsv";
constexpr std::string_view kSampleTrimmedSimplexPassingReadLengthDistr =
    "sample_trimmed_simplex_passing_read_length_distr.tsv";
constexpr std::string_view kTotalUntrimmedDuplexPassingReadLengthDistr =
    "total_untrimmed_duplex_passing_read_length_distr.tsv";
constexpr std::string_view kTotalTrimmedDuplexPassingReadLengthDistr =
    "total_trimmed_duplex_passing_read_length_distr.tsv";
constexpr std::string_view kSampleUntrimmedDuplexPassingReadLengthDistr =
    "sample_untrimmed_duplex_passing_read_length_distr.tsv";
constexpr std::string_view kSampleTrimmedDuplexPassingReadLengthDistr =
    "sample_trimmed_duplex_passing_read_length_distr.tsv";
const constexpr char* kReadLengthDistributionDirectory = "read_length_distr";
constexpr std::string_view kSampleReadLengthDistributionSuffix = "_read_length_distr";

inline std::string SampleReadLengthDistributionDirectory(const std::string_view sample_name) {
  return std::string{sample_name} + std::string{kSampleReadLengthDistributionSuffix};
}

// metrics file names used in demux-and-trim pipeline
constexpr std::string_view kSampleAssignmentMetrics = "sample_assignment_metrics.tsv";

}  // namespace xoos::demux
