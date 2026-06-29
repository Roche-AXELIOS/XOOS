#pragma once

#include <xoos/concurrent/enumerable-thread-local.h>
#include <xoos/histogram/histogram-summary.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

namespace xoos::read_collapser {

// kNA represents metric fields that are not applicable (e.g. consensus metrics will be NA for duplicate marking mode).
constexpr std::string kNA = "NA";
constexpr u16 kDefaultMaxClusterSize = 50;
const vec<u8> kClusterSizesPrecision{
    0,  // total_clusters
    0,  // forward_strand_clusters
    0,  // reverse_strand_clusters
    0,  // mixed_strand_clusters
    2,  // pct_mixed_strand_clusters
    2,  // pct_forward_strand_reads
    0,  // full_read_clusters
    0,  // full_and_partial_read_clusters
    0,  // partial_read_clusters
};

// This enum defines the mode of operation for read collapser, specifically whether it is handling duplicate marking or
// consensus, so that it can be specified when generating and outputting the metrics. The metric output differs between
// modes (e.g. consensus metrics will be NA for duplicate marking mode).
enum class ReadCollapserMode {
  kMarkDuplicate,
  kConsensus
};

// This struct is used to store metrics related to the read collapser operations, such as clustering and consensus.
struct Metrics {
  // output metrics
  u64 input_reads{};
  u64 discarded_missing_umi_reads{};
  u64 discarded_low_mapq_reads{};
  u64 discarded_by_flags_reads{};
  u64 discarded_high_discordant_duplex_percentage_reads{};
  u64 discarded_total_reads{};
  u64 unmapped_reads{};
  u64 unclustered_partial_reads{};
  u64 unclustered_supplementary_reads{};
  u64 unclustered_secondary_reads{};
  u64 clustering_input_reads{};
  u64 clustering_reads{};
  u64 clustering_full_reads{};
  u64 clustering_partial_reads{};
  u64 clustering_unclustered_partial_reads{};
  u64 duplicate_reads{};
  u64 unclustered_reads{};
  u64 total_clusters{};
  u64 singleton_clusters{};
  u64 full_read_clusters{};
  u64 full_and_partial_read_clusters{};
  u64 partial_read_clusters{};
  u64 forward_strand_clusters{};
  u64 reverse_strand_clusters{};
  u64 mixed_strand_clusters{};
  u64 consensus_discarded_by_size_clusters{};
  u64 consensus_discarded_by_subsampling_reads{};
  u64 consensus_input_reads{};
  u64 consensus_input_full_reads{};
  u64 consensus_input_partial_reads{};
  u64 consensus_clusters{};
  u64 consensus_median_cluster_size{};
  u64 consensus_full_read_clusters{};
  u64 consensus_full_and_partial_read_clusters{};
  u64 consensus_partial_read_clusters{};
  u64 consensus_forward_strand_clusters{};
  u64 consensus_reverse_strand_clusters{};
  u64 consensus_mixed_strand_clusters{};
  u64 consensus_discarded_by_length_reads{};
  u64 consensus_discarded_by_same_strand_cluster_size_reads{};
  u64 consensus_discarded_by_mixed_strand_cluster_size_reads{};
  u64 total_consensus_reads{};

  // metric to help calculate consensus_median_cluster_size
  histogram::Histogram<u64> consensus_cluster_sizes{histogram::Histogram<u64>(kDefaultMaxClusterSize, 0)};

  // histogram::Histograms are used to store the cluster size distributions
  histogram::Histograms<u64> cluster_sizes{
      std::make_tuple("total_clusters", histogram::Histogram{histogram::Histogram<u64>(kDefaultMaxClusterSize, 0)}),
      std::make_tuple("forward_strand_clusters",
                      histogram::Histogram{histogram::Histogram<u64>(kDefaultMaxClusterSize, 0)}),
      std::make_tuple("reverse_strand_clusters",
                      histogram::Histogram{histogram::Histogram<u64>(kDefaultMaxClusterSize, 0)}),
      std::make_tuple("mixed_strand_clusters",
                      histogram::Histogram{histogram::Histogram<u64>(kDefaultMaxClusterSize, 0)}),
      std::make_tuple("full_read_clusters", histogram::Histogram{histogram::Histogram<u64>(kDefaultMaxClusterSize, 0)}),
      std::make_tuple("full_and_partial_read_clusters",
                      histogram::Histogram{histogram::Histogram<u64>(kDefaultMaxClusterSize, 0)}),
      std::make_tuple("partial_read_clusters",
                      histogram::Histogram{histogram::Histogram<u64>(kDefaultMaxClusterSize, 0)}),
      std::make_tuple("forward_strand_reads",
                      histogram::Histogram{histogram::Histogram<u64>(kDefaultMaxClusterSize, 0)}),
      std::make_tuple("total_reads", histogram::Histogram<u64>{histogram::Histogram<u64>(kDefaultMaxClusterSize, 0)}),
  };

  void Reset();

  void UpdateClusterSizeHistogram(u64 cluster_size, const std::string& distribution_name, u64 value);

  void WriteSummaryStatsToTsv(const fs::path& output,
                              const io::Comments& comments,
                              ReadCollapserMode use_case,
                              bool cluster_by_umi) const;

  void WriteClusterSizeHistogramsToTsv(const fs::path& output, const io::Comments& comments) const;

  void WriteClusterSizeHistogramSummariesToTsv(const fs::path& output, const io::Comments& comments) const;

  void WriteAllMetricsToTsv(const fs::path& output_dir,
                            ReadCollapserMode use_case,
                            bool cluster_by_umi,
                            const std::string& version,
                            const std::string& command_line) const;
};

Metrics& operator+=(Metrics& lhs, const Metrics& rhs);

// This class is used to store metrics in a thread-local manner.
class ConcurrentMetrics {
 public:
  // Get the thread-local instance of Metrics.
  static Metrics& Get();
  // Get the total metrics across all threads.
  static Metrics GetTotal();
  // Reset the metrics for the current thread.
  static void Reset();

 private:
  static thread_local concurrent::EnumerableThreadLocal<Metrics> metrics;
};

}  // namespace xoos::read_collapser
