#include "metrics/metrics.h"

#include <fstream>

#include <core/read-collapser-options.h>

#include <csv.hpp>

#include <xoos/io/metadata-util.h>
#include <xoos/types/float.h>
#include <xoos/util/file-functions.h>

#include "io/output-path.h"
#include "metrics/metrics-format-util.h"
#include "xoos/log/logging.h"

namespace xoos::read_collapser {

// Overload the += operator to add two Metrics objects together.
Metrics& operator+=(Metrics& lhs, const Metrics& rhs) {
  lhs.input_records += rhs.input_records;
  lhs.discarded_missing_umi_records += rhs.discarded_missing_umi_records;
  lhs.discarded_low_mapq_records += rhs.discarded_low_mapq_records;
  lhs.discarded_by_flags_records += rhs.discarded_by_flags_records;
  lhs.discarded_high_discordant_duplex_percentage_records += rhs.discarded_high_discordant_duplex_percentage_records;
  lhs.discarded_total_records += rhs.discarded_total_records;
  lhs.rescued_secondary_alignments += rhs.rescued_secondary_alignments;
  lhs.unmapped_reads += rhs.unmapped_reads;
  lhs.unclustered_partial_reads += rhs.unclustered_partial_reads;
  lhs.unclustered_supplementary_alignments += rhs.unclustered_supplementary_alignments;
  lhs.unclustered_secondary_alignments += rhs.unclustered_secondary_alignments;
  lhs.clustering_input_reads += rhs.clustering_input_reads;
  lhs.clustering_reads += rhs.clustering_reads;
  lhs.clustering_full_reads += rhs.clustering_full_reads;
  lhs.clustering_partial_reads += rhs.clustering_partial_reads;
  lhs.clustering_unclustered_partial_reads += rhs.clustering_unclustered_partial_reads;
  lhs.duplicate_reads += rhs.duplicate_reads;
  lhs.duplicate_supplementary_alignments += rhs.duplicate_supplementary_alignments;
  lhs.unclustered_reads += rhs.unclustered_reads;
  lhs.total_clusters += rhs.total_clusters;
  lhs.singleton_clusters += rhs.singleton_clusters;
  lhs.full_read_clusters += rhs.full_read_clusters;
  lhs.full_and_partial_read_clusters += rhs.full_and_partial_read_clusters;
  lhs.partial_read_clusters += rhs.partial_read_clusters;
  lhs.forward_strand_clusters += rhs.forward_strand_clusters;
  lhs.reverse_strand_clusters += rhs.reverse_strand_clusters;
  lhs.mixed_strand_clusters += rhs.mixed_strand_clusters;
  lhs.consensus_discarded_by_size_clusters += rhs.consensus_discarded_by_size_clusters;
  lhs.consensus_discarded_by_subsampling_reads += rhs.consensus_discarded_by_subsampling_reads;
  lhs.consensus_input_reads += rhs.consensus_input_reads;
  lhs.consensus_input_full_reads += rhs.consensus_input_full_reads;
  lhs.consensus_input_partial_reads += rhs.consensus_input_partial_reads;
  lhs.consensus_clusters += rhs.consensus_clusters;
  lhs.consensus_median_cluster_size += rhs.consensus_median_cluster_size;
  lhs.consensus_full_read_clusters += rhs.consensus_full_read_clusters;
  lhs.consensus_full_and_partial_read_clusters += rhs.consensus_full_and_partial_read_clusters;
  lhs.consensus_partial_read_clusters += rhs.consensus_partial_read_clusters;
  lhs.consensus_forward_strand_clusters += rhs.consensus_forward_strand_clusters;
  lhs.consensus_reverse_strand_clusters += rhs.consensus_reverse_strand_clusters;
  lhs.consensus_mixed_strand_clusters += rhs.consensus_mixed_strand_clusters;
  lhs.consensus_discarded_by_length_reads += rhs.consensus_discarded_by_length_reads;
  lhs.consensus_discarded_by_same_strand_cluster_size_reads +=
      rhs.consensus_discarded_by_same_strand_cluster_size_reads;
  lhs.consensus_discarded_by_mixed_strand_cluster_size_reads +=
      rhs.consensus_discarded_by_mixed_strand_cluster_size_reads;
  lhs.total_consensus_reads += rhs.total_consensus_reads;

  lhs.consensus_cluster_sizes += rhs.consensus_cluster_sizes;
  // for each distribution in cluster_sizes, add the counts for each bin from rhs to lhs
  lhs.cluster_sizes += rhs.cluster_sizes;
  return lhs;
}

// Update the cluster size histogram for a specific distribution name with the provided value.
void Metrics::UpdateClusterSizeHistogram(const u64 cluster_size,
                                         const std::string& distribution_name,
                                         const u64 value) {
  for (auto& [name, histogram] : cluster_sizes) {
    if (name == distribution_name) {
      histogram.AddCountToHistogram(cluster_size, value);
      return;
    }
  }
}

namespace {

// helper function to calculate a percentage histogram by dividing the cluster_sizes histograms specified by the
// denominator and numerator
histogram::Histogram<f64> CalculatePercentageHistogram(const histogram::Histograms<u64>& cluster_sizes,
                                                       const std::string& denominator,
                                                       const std::string& numerator) {
  if (cluster_sizes.empty()) {
    return histogram::Histogram<f64>{};
  }
  histogram::Histogram<f64> percentage_histogram(std::get<1>(cluster_sizes.at(0)).counts.size(), 0);
  const auto denominator_histogram = std::ranges::find_if(
      cluster_sizes, [&denominator](const auto& tuple) { return std::get<0>(tuple) == denominator; });
  const auto numerator_histogram =
      std::ranges::find_if(cluster_sizes, [&numerator](const auto& tuple) { return std::get<0>(tuple) == numerator; });

  // processing bin separately from outliers for efficiency
  for (u64 bin_index = 0; bin_index < std::get<1>(*denominator_histogram).counts.size(); ++bin_index) {
    if (std::get<1>(*denominator_histogram).counts.at(bin_index) > 0) {
      const auto numerator_bin_count = static_cast<f64>(std::get<1>(*numerator_histogram).counts.at(bin_index));
      const auto denominator_bin_count = static_cast<f64>(std::get<1>(*denominator_histogram).counts.at(bin_index));
      const f64 percentage = numerator_bin_count / denominator_bin_count * 100.0;
      percentage_histogram.AddCountToHistogram(bin_index, percentage);
    } else {
      percentage_histogram.AddCountToHistogram(bin_index, 0);
    }
  }

  u64 numerator_count{0};
  u64 denominator_count{0};
  // for the outliers, since we don't know the output bin index, we just store the counts, and add them to the end since
  // we KNOW that this method is ONLY used for the max bin
  for (const auto& [bin_index, denominator_bin_count] : std::get<1>(*denominator_histogram).outliers) {
    // in some cases the denominator histogram may not have the same outliers as the numerator histogram, so we
    if (std::get<1>(*numerator_histogram).outliers.contains(bin_index)) {
      numerator_count += std::get<1>(*numerator_histogram).outliers.at(bin_index);
    }
    denominator_count += denominator_bin_count;
  }
  // now we add it to the final bin for each histogram
  if (denominator_count > 0) {
    const f64 percentage = static_cast<f64>(numerator_count) / static_cast<f64>(denominator_count) * 100.0;
    percentage_histogram.AddCountToHistogram(percentage_histogram.counts.size(), percentage);
  }

  return percentage_histogram;
}

}  // namespace

// Write the summary statistics to a TSV file.
void Metrics::WriteSummaryStatsToTsv(const fs::path& output,
                                     const io::TsvMetadata& tsv_metadata,
                                     const ReadCollapserMode use_case,
                                     const ReadCollapserOptions& options) const {
  auto out = std::ofstream{output};
  auto writer = csv::make_tsv_writer(out).set_auto_flush(false);
  io::WriteTsvMetadata(out, tsv_metadata);

  // write header
  writer << vec<std::string>{"metric_name", "value", "percentage", "denominator"};

  // write data
  writer << FormatRow("input_records", input_records, 0, kNA);
  // If umi, discarded_missing_umi_records metric value will be applicable
  writer << FormatRow("discarded_missing_umi_records",
                      discarded_missing_umi_records,
                      input_records,
                      "input_records",
                      options.cluster_by_umi);
  writer << FormatRow("discarded_low_mapq_records", discarded_low_mapq_records, input_records, "input_records");
  writer << FormatRow("discarded_by_flags_records", discarded_by_flags_records, input_records, "input_records");
  writer << FormatRow("discarded_high_discordant_duplex_percentage_records",
                      discarded_high_discordant_duplex_percentage_records,
                      input_records,
                      "input_records");
  writer << FormatRow("discarded_total_records", discarded_total_records, input_records, "input_records");
  writer << FormatRow("rescued_secondary_alignments",
                      rescued_secondary_alignments,
                      input_records,
                      "input_records",
                      use_case == ReadCollapserMode::kMarkDuplicate);
  writer << FormatRow("unmapped_reads", unmapped_reads, input_records, "input_records");
  writer << FormatRow("unclustered_partial_reads", unclustered_partial_reads, input_records, "input_records");
  writer << FormatRow(
      "unclustered_supplementary_alignments", unclustered_supplementary_alignments, input_records, "input_records");
  writer << FormatRow(
      "unclustered_secondary_alignments", unclustered_secondary_alignments, input_records, "input_records");
  writer << FormatRow("clustering_input_reads", clustering_input_reads, input_records, "input_records");
  writer << FormatRow("clustering_reads", clustering_reads, 0, kNA);
  writer << FormatRow("clustering_full_reads", clustering_full_reads, clustering_reads, "clustering_reads");
  writer << FormatRow("clustering_partial_reads", clustering_partial_reads, clustering_reads, "clustering_reads");
  writer << FormatRow("clustering_unclustered_partial_reads",
                      clustering_unclustered_partial_reads,
                      clustering_reads,
                      "clustering_reads");
  writer << FormatRow("duplicate_reads", duplicate_reads, clustering_input_reads, "clustering_input_reads");
  writer << FormatRow(
      "duplicate_supplementary_alignments",
      duplicate_supplementary_alignments,
      unclustered_supplementary_alignments,
      "unclustered_supplementary_alignments",
      use_case == ReadCollapserMode::kMarkDuplicate && options.markdup_options.mark_supplementary_alignments);
  writer << FormatRow("unclustered_reads", unclustered_reads, input_records, "input_records");
  writer << FormatRow("total_clusters", total_clusters, 0, kNA);
  writer << FormatRow("singleton_clusters", singleton_clusters, total_clusters, "total_clusters");
  writer << FormatRow("full_read_clusters", full_read_clusters, total_clusters, "total_clusters");
  writer << FormatRow(
      "full_and_partial_read_clusters", full_and_partial_read_clusters, total_clusters, "total_clusters");
  writer << FormatRow("partial_read_clusters", partial_read_clusters, total_clusters, "total_clusters");
  writer << FormatRow("forward_strand_clusters", forward_strand_clusters, total_clusters, "total_clusters");
  writer << FormatRow("reverse_strand_clusters", reverse_strand_clusters, total_clusters, "total_clusters");
  writer << FormatRow("mixed_strand_clusters", mixed_strand_clusters, total_clusters, "total_clusters");
  // If not duplicate marking, consensus will be performed and the corresponding metric values will be applicable
  writer << FormatRow("consensus_discarded_by_size_clusters",
                      consensus_discarded_by_size_clusters,
                      total_clusters,
                      "total_clusters",
                      use_case == ReadCollapserMode::kConsensus);
  writer << FormatRow("consensus_discarded_by_subsampling_reads",
                      consensus_discarded_by_subsampling_reads,
                      clustering_reads,
                      "clustering_reads",
                      use_case == ReadCollapserMode::kConsensus);
  writer << FormatRow(
      "consensus_input_reads", consensus_input_reads, 0, kNA, use_case == ReadCollapserMode::kConsensus);
  writer << FormatRow("consensus_input_full_reads",
                      consensus_input_full_reads,
                      consensus_input_reads,
                      "consensus_input_reads",
                      use_case == ReadCollapserMode::kConsensus);
  writer << FormatRow("consensus_input_partial_reads",
                      consensus_input_partial_reads,
                      consensus_input_reads,
                      "consensus_input_reads",
                      use_case == ReadCollapserMode::kConsensus);
  writer << FormatRow("consensus_clusters", consensus_clusters, 0, kNA, use_case == ReadCollapserMode::kConsensus);
  // consensus_median_cluster_size is a statistic (median) rather than a count, and it has no named
  // denominator population for FormatRow's cascade to key on. A median over zero clusters is
  // undefined, so gate it explicitly on there being at least one consensus cluster; otherwise it is
  // reported as NA rather than a misleading 0.
  writer << FormatRow("consensus_median_cluster_size",
                      consensus_median_cluster_size,
                      0,
                      kNA,
                      use_case == ReadCollapserMode::kConsensus && consensus_clusters > 0);
  writer << FormatRow("consensus_full_read_clusters",
                      consensus_full_read_clusters,
                      consensus_clusters,
                      "consensus_clusters",
                      use_case == ReadCollapserMode::kConsensus);
  writer << FormatRow("consensus_full_and_partial_read_clusters",
                      consensus_full_and_partial_read_clusters,
                      consensus_clusters,
                      "consensus_clusters",
                      use_case == ReadCollapserMode::kConsensus);
  writer << FormatRow("consensus_partial_read_clusters",
                      consensus_partial_read_clusters,
                      consensus_clusters,
                      "consensus_clusters",
                      use_case == ReadCollapserMode::kConsensus);
  writer << FormatRow("consensus_forward_strand_clusters",
                      consensus_forward_strand_clusters,
                      consensus_clusters,
                      "consensus_clusters",
                      use_case == ReadCollapserMode::kConsensus);
  writer << FormatRow("consensus_reverse_strand_clusters",
                      consensus_reverse_strand_clusters,
                      consensus_clusters,
                      "consensus_clusters",
                      use_case == ReadCollapserMode::kConsensus);
  writer << FormatRow("consensus_mixed_strand_clusters",
                      consensus_mixed_strand_clusters,
                      consensus_clusters,
                      "consensus_clusters",
                      use_case == ReadCollapserMode::kConsensus);
  writer << FormatRow("consensus_discarded_by_length_reads",
                      consensus_discarded_by_length_reads,
                      0,
                      kNA,
                      use_case == ReadCollapserMode::kConsensus);
  writer << FormatRow("consensus_discarded_by_same_strand_cluster_size_reads",
                      consensus_discarded_by_same_strand_cluster_size_reads,
                      0,
                      kNA,
                      use_case == ReadCollapserMode::kConsensus);
  writer << FormatRow("consensus_discarded_by_mixed_strand_cluster_size_reads",
                      consensus_discarded_by_mixed_strand_cluster_size_reads,
                      0,
                      kNA,
                      use_case == ReadCollapserMode::kConsensus);
  writer << FormatRow(
      "total_consensus_reads", total_consensus_reads, 0, kNA, use_case == ReadCollapserMode::kConsensus);
}

// Write the cluster size histograms to a TSV file. Reformat the data to histograms of type double to permit percentage
// output.
void Metrics::WriteClusterSizeHistogramsToTsv(const fs::path& output, const io::TsvMetadata& tsv_metadata) const {
  vec<std::string> output_histogram_headers{"total_clusters",
                                            "forward_strand_clusters",
                                            "reverse_strand_clusters",
                                            "mixed_strand_clusters",
                                            "pct_mixed_strand_clusters",
                                            "pct_forward_strand_reads",
                                            "full_read_clusters",
                                            "full_and_partial_read_clusters",
                                            "partial_read_clusters"};
  histogram::Histograms<f64> output_histograms;
  for (auto& output_name : output_histogram_headers) {
    for (const auto& [name, histogram] : cluster_sizes) {
      if (output_name == "pct_mixed_strand_clusters") {
        output_histograms.emplace_back(
            output_name, CalculatePercentageHistogram(cluster_sizes, "total_clusters", "mixed_strand_clusters"));
        break;
      }
      if (output_name == "pct_forward_strand_reads") {
        output_histograms.emplace_back(
            output_name, CalculatePercentageHistogram(cluster_sizes, "total_reads", "forward_strand_reads"));
        break;
      }
      if (output_name == name) {
        // all other histograms must statically cast each element to f64
        histogram::Histogram<f64> histogram_double;
        for (const u64 count : histogram.counts) {
          histogram_double.counts.push_back(static_cast<f64>(count));
        }
        // copy the outliers next
        for (const auto& [key, value] : histogram.outliers) {
          histogram_double.outliers.emplace(key, static_cast<f64>(value));
        }
        output_histograms.emplace_back(output_name, histogram_double);
        break;
      }
    }
  }
  Logging::Info("Writing cluster size histograms to TSV file: {}", output.string());
  histogram::WriteHistogramsToTsv(output_histograms,
                                  output,
                                  "cluster_sizes",
                                  histogram::HistogramBinOutput::kOmitFirstBinAndMaxLastBinWithOutlier,
                                  kClusterSizesPrecision,
                                  tsv_metadata);
}

// Write the cluster size histogram summaries to a TSV file. This will calculate the summary statistics for each
// histogram that is not a read count.
void Metrics::WriteClusterSizeHistogramSummariesToTsv(const fs::path& output,
                                                      const io::TsvMetadata& tsv_metadata) const {
  histogram::HistogramSummaries<u64> summary_histograms;
  for (const auto& [name, histogram] : cluster_sizes) {
    // Only calculate summary for histograms that are stats of clusters
    if (name.find("_reads") != std::string::npos) {
      continue;
    }
    summary_histograms.emplace_back(name, histogram::CalculateHistogramSummary(histogram));
  }
  histogram::WriteSummaryHistogramsToTsv(
      summary_histograms, output, histogram::kDefaultPercentiles, {}, true, tsv_metadata);
}

bool Metrics::HasEmptyOutput(const ReadCollapserMode use_case) const {
  if (use_case == ReadCollapserMode::kMarkDuplicate) {
    // markdup writes retained input records to the output BAM. Output is header-only when input is
    // empty or when all input records are discarded by filtering.
    return input_records == 0 || discarded_total_records == input_records;
  }
  // consensus writes one record per consensus read to the output FASTQ.
  return total_consensus_reads == 0;
}

void Metrics::WriteAllMetricsToTsv(const ReadCollapserMode use_case,
                                   const ReadCollapserOptions& options,
                                   const io::TsvMetadata& tsv_metadata) const {
  // Warn (but exit successfully) when the run produced no output records.
  if (HasEmptyOutput(use_case)) {
    if (use_case == ReadCollapserMode::kMarkDuplicate) {
      Logging::Warn(
          "No records were written to output BAM; output BAM is header-only. "
          "This can happen when the input BAM has zero records or all records were discarded by filters.");
    } else {
      Logging::Warn(
          "No consensus reads were produced; output FASTQ is empty. "
          "This can happen when the input BAM has zero records or all clusters were filtered out.");
    }
  }

  // Metrics are written under output-dir/metrics/. The BAM/FASTQ writers create their own parent
  // directories, but the metrics writers open std::ofstream directly, so create the subdirectory here.
  file::CreateWritableDirectory(options.output_dir / "metrics");
  WriteSummaryStatsToTsv(MetricsOutputPath(options, "summary_stats.tsv"), tsv_metadata, use_case, options);
  WriteClusterSizeHistogramSummariesToTsv(MetricsOutputPath(options, "cluster_size_distribution_summary.tsv"),
                                          tsv_metadata);
  WriteClusterSizeHistogramsToTsv(MetricsOutputPath(options, "cluster_size_distributions.tsv"), tsv_metadata);
}

Metrics SumMetrics(const vec<Metrics>& metrics_subsets) {
  auto total = Metrics{};
  for (const auto& subset : metrics_subsets) {
    total += subset;
  }
  return total;
}

}  // namespace xoos::read_collapser
