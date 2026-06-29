#include "metrics/metrics.h"

#include <fstream>

#include <csv.hpp>

#include <xoos/io/metadata-util.h>
#include <xoos/types/float.h>

#include "xoos/log/logging.h"

namespace xoos::read_collapser {

// Overload the += operator to add two Metrics objects together.
Metrics& operator+=(Metrics& lhs, const Metrics& rhs) {
  lhs.input_reads += rhs.input_reads;
  lhs.discarded_missing_umi_reads += rhs.discarded_missing_umi_reads;
  lhs.discarded_low_mapq_reads += rhs.discarded_low_mapq_reads;
  lhs.discarded_by_flags_reads += rhs.discarded_by_flags_reads;
  lhs.discarded_high_discordant_duplex_percentage_reads += rhs.discarded_high_discordant_duplex_percentage_reads;
  lhs.discarded_total_reads += rhs.discarded_total_reads;
  lhs.unmapped_reads += rhs.unmapped_reads;
  lhs.unclustered_partial_reads += rhs.unclustered_partial_reads;
  lhs.unclustered_supplementary_reads += rhs.unclustered_supplementary_reads;
  lhs.unclustered_secondary_reads += rhs.unclustered_secondary_reads;
  lhs.clustering_input_reads += rhs.clustering_input_reads;
  lhs.clustering_reads += rhs.clustering_reads;
  lhs.clustering_full_reads += rhs.clustering_full_reads;
  lhs.clustering_partial_reads += rhs.clustering_partial_reads;
  lhs.clustering_unclustered_partial_reads += rhs.clustering_unclustered_partial_reads;
  lhs.duplicate_reads += rhs.duplicate_reads;
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

// Reset the metrics to their initial state.
void Metrics::Reset() {
  input_reads = 0;
  discarded_missing_umi_reads = 0;
  discarded_low_mapq_reads = 0;
  discarded_by_flags_reads = 0;
  discarded_high_discordant_duplex_percentage_reads = 0;
  discarded_total_reads = 0;
  unmapped_reads = 0;
  unclustered_partial_reads = 0;
  unclustered_supplementary_reads = 0;
  unclustered_secondary_reads = 0;
  clustering_input_reads = 0;
  clustering_reads = 0;
  clustering_full_reads = 0;
  clustering_partial_reads = 0;
  clustering_unclustered_partial_reads = 0;
  duplicate_reads = 0;
  unclustered_reads = 0;
  total_clusters = 0;
  singleton_clusters = 0;
  full_read_clusters = 0;
  full_and_partial_read_clusters = 0;
  partial_read_clusters = 0;
  forward_strand_clusters = 0;
  reverse_strand_clusters = 0;
  mixed_strand_clusters = 0;
  consensus_discarded_by_size_clusters = 0;
  consensus_discarded_by_subsampling_reads = 0;
  consensus_input_reads = 0;
  consensus_input_full_reads = 0;
  consensus_input_partial_reads = 0;
  consensus_clusters = 0;
  consensus_median_cluster_size = 0;
  consensus_full_read_clusters = 0;
  consensus_full_and_partial_read_clusters = 0;
  consensus_partial_read_clusters = 0;
  consensus_forward_strand_clusters = 0;
  consensus_reverse_strand_clusters = 0;
  consensus_mixed_strand_clusters = 0;
  consensus_discarded_by_length_reads = 0;
  consensus_discarded_by_same_strand_cluster_size_reads = 0;
  consensus_discarded_by_mixed_strand_cluster_size_reads = 0;
  total_consensus_reads = 0;

  consensus_cluster_sizes.Reset();

  for (auto& [name, histogram] : cluster_sizes) {
    histogram.Reset();
  }
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

// helper function to calculate percentage in appropriate string format
// must use stringstream to get the correct precision
std::string CalculatePercentage(u64 value, u64 total, u8 precision = 2) {
  std::stringstream ss;
  ss << std::fixed << std::setprecision(precision) << static_cast<f64>(value) / static_cast<f64>(total) * 100.0;
  return ss.str();
}

// helper function to format the row for the tsv output
vec<std::string> FormatRow(const std::string& metric_name,
                           const u64 value,
                           const u64 denominator,
                           const std::string& denominator_name,
                           const bool metric_not_na = true) {
  if (metric_not_na) {
    vec<std::string> output{metric_name, std::to_string(value)};
    // add the percentage if there is a valid denominator, otherwise add kNA
    if (denominator_name != kNA && denominator != 0) {
      output.emplace_back(CalculatePercentage(value, denominator));
    } else {
      output.push_back(kNA);
    }
    // add the denominator name
    output.push_back(denominator_name);
    return output;
  }
  return {metric_name, kNA, kNA, denominator_name};
}

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

// Write the summary statistics to a TSV file.
void Metrics::WriteSummaryStatsToTsv(const fs::path& output,
                                     const io::Comments& comments,
                                     const ReadCollapserMode use_case,
                                     const bool cluster_by_umi) const {
  auto out = std::ofstream{output};
  auto writer = csv::make_tsv_writer_buffered(out);
  io::WriteTsvComments(writer, comments);

  // write header
  writer << vec<std::string>{"metric_name", "value", "percentage", "denominator"};

  // write data
  writer << FormatRow("input_reads", input_reads, 0, kNA);
  // If umi, discarded_missing_umi_reads metric value will be applicable
  writer << FormatRow(
      "discarded_missing_umi_reads", discarded_missing_umi_reads, input_reads, "input_reads", cluster_by_umi);
  writer << FormatRow("discarded_low_mapq_reads", discarded_low_mapq_reads, input_reads, "input_reads");
  writer << FormatRow("discarded_by_flags_reads", discarded_by_flags_reads, input_reads, "input_reads");
  writer << FormatRow("discarded_high_discordant_duplex_percentage_reads",
                      discarded_high_discordant_duplex_percentage_reads,
                      input_reads,
                      "input_reads");
  writer << FormatRow("discarded_total_reads", discarded_total_reads, input_reads, "input_reads");
  writer << FormatRow("unmapped_reads", unmapped_reads, input_reads, "input_reads");
  writer << FormatRow("unclustered_partial_reads", unclustered_partial_reads, input_reads, "input_reads");
  writer << FormatRow("unclustered_supplementary_reads", unclustered_supplementary_reads, input_reads, "input_reads");
  writer << FormatRow("unclustered_secondary_reads", unclustered_secondary_reads, input_reads, "input_reads");
  writer << FormatRow("clustering_input_reads", clustering_input_reads, input_reads, "input_reads");
  writer << FormatRow("clustering_reads", clustering_reads, 0, kNA);
  writer << FormatRow("clustering_full_reads", clustering_full_reads, clustering_reads, "clustering_reads");
  writer << FormatRow("clustering_partial_reads", clustering_partial_reads, clustering_reads, "clustering_reads");
  writer << FormatRow("clustering_unclustered_partial_reads",
                      clustering_unclustered_partial_reads,
                      clustering_reads,
                      "clustering_reads");
  writer << FormatRow("duplicate_reads", duplicate_reads, clustering_input_reads, "clustering_input_reads");
  writer << FormatRow("unclustered_reads", unclustered_reads, input_reads, "input_reads");
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
  writer << FormatRow("consensus_median_cluster_size",
                      consensus_median_cluster_size,
                      0,
                      kNA,
                      use_case == ReadCollapserMode::kConsensus);
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
void Metrics::WriteClusterSizeHistogramsToTsv(const fs::path& output, const io::Comments& comments) const {
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
    for (auto& [name, histogram] : cluster_sizes) {
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
  histogram::WriteHistograms(output_histograms,
                             output,
                             "cluster_sizes",
                             histogram::HistogramBinOutput::kOmitFirstBinAndMaxLastBinWithOutlier,
                             kClusterSizesPrecision,
                             comments);
}

// Write the cluster size histogram summaries to a TSV file. This will calculate the summary statistics for each
// histogram that is not a read count.
void Metrics::WriteClusterSizeHistogramSummariesToTsv(const fs::path& output, const io::Comments& comments) const {
  histogram::HistogramSummaries<u64> summary_histograms;
  for (auto& [name, histogram] : cluster_sizes) {
    // Only calculate summary for histograms that are stats of clusters
    if (name.find("_reads") != std::string::npos) {
      continue;
    }
    summary_histograms.emplace_back(name, histogram::CalculateHistogramSummary(histogram));
  }
  histogram::WriteSummaryHistograms(summary_histograms, output, histogram::kDefaultPercentiles, {}, true, comments);
}

void Metrics::WriteAllMetricsToTsv(const fs::path& output_dir,
                                   const ReadCollapserMode use_case,
                                   const bool cluster_by_umi,
                                   const std::string& version,
                                   const std::string& command_line) const {
  io::Comments comments;
  io::AddVersionAndCommandLineComment(comments, version, command_line);
  WriteSummaryStatsToTsv(output_dir / "summary_stats.tsv", comments, use_case, cluster_by_umi);
  WriteClusterSizeHistogramSummariesToTsv(output_dir / "cluster_size_distribution_summary.tsv", comments);
  WriteClusterSizeHistogramsToTsv(output_dir / "cluster_size_distributions.tsv", comments);
}

thread_local concurrent::EnumerableThreadLocal<Metrics> ConcurrentMetrics::metrics{std::make_shared<Metrics>()};

Metrics& ConcurrentMetrics::Get() {
  return *metrics.Local();
}

Metrics ConcurrentMetrics::GetTotal() {
  auto total = Metrics{};
  metrics.ForEach([&total](const Metrics& item) { total += item; });
  return total;
}

void ConcurrentMetrics::Reset() {
  metrics.ForEachNonConst([](Metrics& item) { item.Reset(); });
}

}  // namespace xoos::read_collapser
