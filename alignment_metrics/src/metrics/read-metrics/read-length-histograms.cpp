#include "metrics/read-metrics/read-length-histograms.h"

#include <filesystem>
#include <tuple>
#include <vector>

#include <csv.hpp>

#include <xoos/types/int.h>

#include "metrics/metrics-names.h"

namespace xoos::alignment_metrics {

namespace {

// Percentiles always computed for read-length summaries (required for downstream ratio metrics).
constexpr u64 kPercentile10 = 10;
constexpr u64 kPercentile50 = 50;
constexpr u64 kPercentile90 = 90;

/**
 * A set of (all, partial, full) read length histograms together with their column names.
 *
 * The same struct is used for both the mapped (before read-level filters) and post-filter
 * (passing read-level filters) histogram groups.
 */
struct ReadLengthHistogramGroup {
  // This is a short-lived, non-owning view over histograms owned by the caller, used only to pass a
  // group of histograms to WriteReadLengthHistograms; the references are intentional.
  // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
  const histogram::Histogram<u64>& all_histogram;
  const histogram::Histogram<u64>& partial_histogram;
  const histogram::Histogram<u64>& full_histogram;
  // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)
  std::string all_name;
  std::string partial_name;
  std::string full_name;
};

/**
 * Write read length histograms to a TSV file.
 *
 * The mapped group columns are written first, followed by the post-filter group columns. For datasets without read
 * type information, only the "all" column of each group is written.
 *
 * @param mapped_group histograms and column names for mapped reads (before read-level filters)
 * @param post_filter_group histograms and column names for reads passing read-level filters
 * @param has_read_type_info whether partial/full columns should be included
 * @param binning_name the label for the first (x-axis) column
 */
void WriteReadLengthHistograms(const ReadLengthHistogramGroup& mapped_group,
                               const ReadLengthHistogramGroup& post_filter_group,
                               const bool has_read_type_info,
                               const std::string& binning_name,
                               const fs::path& output_path,
                               const io::TsvMetadata& tsv_metadata) {
  histogram::Histograms histograms{std::make_tuple(mapped_group.all_name, mapped_group.all_histogram)};
  if (has_read_type_info) {
    histograms.emplace_back(mapped_group.partial_name, mapped_group.partial_histogram);
    histograms.emplace_back(mapped_group.full_name, mapped_group.full_histogram);
  }
  histograms.emplace_back(post_filter_group.all_name, post_filter_group.all_histogram);
  if (has_read_type_info) {
    histograms.emplace_back(post_filter_group.partial_name, post_filter_group.partial_histogram);
    histograms.emplace_back(post_filter_group.full_name, post_filter_group.full_histogram);
  }
  histogram::WriteHistogramsToTsv(histograms,
                                  output_path,
                                  binning_name,
                                  histogram::HistogramBinOutput::kOmitFirstBinAndMaxLastBin,
                                  {},
                                  tsv_metadata);
}

/**
 * Append summary statistics for a histogram group to the given summary collection.
 *
 * @param summary_histograms the collection to append to
 * @param group histograms and column names to summarize
 * @param has_read_type_info whether partial/full columns should be included
 * @param required_percentiles percentiles to compute (a superset of the ones written to the output)
 */
void AppendReadLengthSummaries(histogram::HistogramSummaries<u64>& summary_histograms,
                               const ReadLengthHistogramGroup& group,
                               const bool has_read_type_info,
                               const std::vector<u64>& required_percentiles) {
  summary_histograms.emplace_back(group.all_name,
                                  histogram::CalculateHistogramSummary(group.all_histogram, {}, required_percentiles));
  if (has_read_type_info) {
    summary_histograms.emplace_back(
        group.partial_name, histogram::CalculateHistogramSummary(group.partial_histogram, {}, required_percentiles));
    summary_histograms.emplace_back(
        group.full_name, histogram::CalculateHistogramSummary(group.full_histogram, {}, required_percentiles));
  }
}

/**
 * Write a summary of read length histograms to a TSV file.
 *
 * The mapped group columns are written first, followed by the post-filter group columns. For datasets without read
 * type information, only the "all" column of each group is written.
 *
 * @param mapped_group histograms and column names for mapped reads (before read-level filters)
 * @param post_filter_group histograms and column names for reads passing read-level filters
 * @param has_read_type_info whether partial/full columns should be included
 * @param summary_stats_percentiles percentiles requested by the user for output
 */
void WriteReadLengthSummary(const ReadLengthHistogramGroup& mapped_group,
                            const ReadLengthHistogramGroup& post_filter_group,
                            const bool has_read_type_info,
                            const std::vector<u64>& summary_stats_percentiles,
                            const fs::path& output_path,
                            const io::TsvMetadata& tsv_metadata) {
  // Make sure we always calculate the 10th, 50th, and 90th percentile
  // because we need these values for the ratio calculations.
  // NOTE: These percentiles are not necessarily written to the output file
  // unless they are explicitly requested by the user.
  std::vector<u64> required_percentiles = {kPercentile10, kPercentile50, kPercentile90};
  // Add any additional percentiles specified in the summary_stats_percentiles
  for (const auto& percentile : summary_stats_percentiles) {
    if (std::ranges::find(required_percentiles, percentile) == required_percentiles.end()) {
      required_percentiles.emplace_back(percentile);
    }
  }
  // Calculate summary statistics for each histogram group.
  // Include required percentiles in addition to the summary_stats_percentiles.
  histogram::HistogramSummaries<u64> summary_histograms;
  AppendReadLengthSummaries(summary_histograms, mapped_group, has_read_type_info, required_percentiles);
  AppendReadLengthSummaries(summary_histograms, post_filter_group, has_read_type_info, required_percentiles);
  // When writing to the output, we can omit the required percentiles if they are not
  // specifically requested by the user. The values are still implicitly stored in
  // the HistogramSummary struct for computing percentile ratios.
  histogram::WriteSummaryHistogramsToTsv(
      summary_histograms, output_path, summary_stats_percentiles, {}, false, tsv_metadata);
}

}  // namespace

void ReadLengthHistograms::WriteTsv(const fs::path& output_path, const io::TsvMetadata& tsv_metadata) const {
  WriteReadLengthHistograms(ReadLengthHistogramGroup{.all_histogram = mapped_all_read_length_histogram,
                                                     .partial_histogram = mapped_partial_read_length_histogram,
                                                     .full_histogram = mapped_full_read_length_histogram,
                                                     .all_name = kNameMappedAllReadCount,
                                                     .partial_name = kNameMappedPartialReadCount,
                                                     .full_name = kNameMappedFullReadCount},
                            ReadLengthHistogramGroup{.all_histogram = post_filter_all_read_length_histogram,
                                                     .partial_histogram = post_filter_partial_read_length_histogram,
                                                     .full_histogram = post_filter_full_read_length_histogram,
                                                     .all_name = kNamePostFilterAllReadCount,
                                                     .partial_name = kNamePostFilterPartialReadCount,
                                                     .full_name = kNamePostFilterFullReadCount},
                            dataset_metadata.has_read_type_info,
                            kNameReadLength,
                            output_path,
                            tsv_metadata);
}

void ReadLengthHistograms::WriteExcludingSoftClipsTsv(const fs::path& output_path,
                                                      const io::TsvMetadata& tsv_metadata) const {
  WriteReadLengthHistograms(
      ReadLengthHistogramGroup{.all_histogram = mapped_all_read_length_excluding_softclips_histogram,
                               .partial_histogram = mapped_partial_read_length_excluding_softclips_histogram,
                               .full_histogram = mapped_full_read_length_excluding_softclips_histogram,
                               .all_name = kNameMappedAllReadCount,
                               .partial_name = kNameMappedPartialReadCount,
                               .full_name = kNameMappedFullReadCount},
      ReadLengthHistogramGroup{.all_histogram = post_filter_all_read_length_excluding_softclips_histogram,
                               .partial_histogram = post_filter_partial_read_length_excluding_softclips_histogram,
                               .full_histogram = post_filter_full_read_length_excluding_softclips_histogram,
                               .all_name = kNamePostFilterAllReadCount,
                               .partial_name = kNamePostFilterPartialReadCount,
                               .full_name = kNamePostFilterFullReadCount},
      dataset_metadata.has_read_type_info,
      kNameReadLengthExcludingSoftClips,
      output_path,
      tsv_metadata);
}

void ReadLengthHistograms::WriteSummaryTsv(const fs::path& output_path, const io::TsvMetadata& tsv_metadata) const {
  WriteReadLengthSummary(ReadLengthHistogramGroup{.all_histogram = mapped_all_read_length_histogram,
                                                  .partial_histogram = mapped_partial_read_length_histogram,
                                                  .full_histogram = mapped_full_read_length_histogram,
                                                  .all_name = kNameMappedAllReadsLength,
                                                  .partial_name = kNameMappedPartialReadsLength,
                                                  .full_name = kNameMappedFullReadsLength},
                         ReadLengthHistogramGroup{.all_histogram = post_filter_all_read_length_histogram,
                                                  .partial_histogram = post_filter_partial_read_length_histogram,
                                                  .full_histogram = post_filter_full_read_length_histogram,
                                                  .all_name = kNamePostFilterAllReadsLength,
                                                  .partial_name = kNamePostFilterPartialReadsLength,
                                                  .full_name = kNamePostFilterFullReadsLength},
                         dataset_metadata.has_read_type_info,
                         summary_stats_percentiles,
                         output_path,
                         tsv_metadata);
}

void ReadLengthHistograms::WriteExcludingSoftClipsSummaryTsv(const fs::path& output_path,
                                                             const io::TsvMetadata& tsv_metadata) const {
  WriteReadLengthSummary(
      ReadLengthHistogramGroup{.all_histogram = mapped_all_read_length_excluding_softclips_histogram,
                               .partial_histogram = mapped_partial_read_length_excluding_softclips_histogram,
                               .full_histogram = mapped_full_read_length_excluding_softclips_histogram,
                               .all_name = kNameMappedAllReadsLengthExcludingSoftClips,
                               .partial_name = kNameMappedPartialReadsLengthExcludingSoftClips,
                               .full_name = kNameMappedFullReadsLengthExcludingSoftClips},
      ReadLengthHistogramGroup{.all_histogram = post_filter_all_read_length_excluding_softclips_histogram,
                               .partial_histogram = post_filter_partial_read_length_excluding_softclips_histogram,
                               .full_histogram = post_filter_full_read_length_excluding_softclips_histogram,
                               .all_name = kNamePostFilterAllReadsLengthExcludingSoftClips,
                               .partial_name = kNamePostFilterPartialReadsLengthExcludingSoftClips,
                               .full_name = kNamePostFilterFullReadsLengthExcludingSoftClips},
      dataset_metadata.has_read_type_info,
      summary_stats_percentiles,
      output_path,
      tsv_metadata);
}

}  // namespace xoos::alignment_metrics
