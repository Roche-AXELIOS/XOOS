#include "metrics/metrics.h"

#include <filesystem>
#include <optional>

#include <xoos/io/metadata-util.h>
#include <xoos/types/fs.h>
#include <xoos/util/file-functions.h>

#include "alignment-metrics-options.h"
#include "io/metrics-filenames.h"
#include "metadata/dataset-metadata.h"

namespace xoos::alignment_metrics {

Metrics::Metrics(const AlignmentMetricsOptions& options, const DatasetMetadata& dataset_metadata)
    : accuracy_metrics(options.metric_types.has_accuracy_metrics
                           ? std::make_optional(AccuracyMetrics(options, dataset_metadata))
                           : std::nullopt),
      coverage_metrics(options.metric_types.has_coverage_metrics
                           ? std::make_optional(CoverageMetrics(options, dataset_metadata))
                           : std::nullopt),
      read_metrics(options.metric_types.has_read_metrics ? std::make_optional(ReadMetrics(options, dataset_metadata))
                                                         : std::nullopt),
      gc_bias_metrics(options.metric_types.has_gc_bias_metrics ? std::make_optional(GcBiasMetrics()) : std::nullopt) {
}

void Metrics::WriteMetrics(const fs::path& output_dir,
                           const io::TsvMetadata& tsv_metadata,
                           const std::string& output_prefix) const {
  using file::PrefixedOutputPath;
  const fs::path absolute_out_dir = fs::absolute(output_dir);
  // Create output directories if it doesn't exist
  // These calls to `create_directories` will not throw even if the directories already exist and will simply return
  // false and do nothing.
  fs::create_directories(absolute_out_dir);
  const auto accuracy_dir = output_dir / kDefaultAccuracyDir;
  const auto coverage_dir = output_dir / kDefaultCoverageDir;
  const auto read_metrics_dir = output_dir / kDefaultReadMetricsDir;
  if (accuracy_metrics.has_value()) {
    create_directories(absolute_out_dir / kDefaultAccuracyDir);
    accuracy_metrics->errors_by_substitution_type.WriteTsv(
        PrefixedOutputPath(accuracy_dir, output_prefix, kDefaultErrorBySubstitutionType), tsv_metadata);
    accuracy_metrics->base_level_accuracy_summary.WriteTsv(
        PrefixedOutputPath(accuracy_dir, output_prefix, kDefaultBaseLevelAccuracySummary), tsv_metadata);
    accuracy_metrics->qscore_stats.WriteTsv(PrefixedOutputPath(accuracy_dir, output_prefix, kDefaultQscoreStats),
                                            tsv_metadata);
    if (accuracy_metrics->errors_by_cluster_size.has_value()) {
      accuracy_metrics->errors_by_cluster_size->WriteTsv(
          PrefixedOutputPath(accuracy_dir, output_prefix, kDefaultErrorByClusterSize), tsv_metadata);
    }
    if (accuracy_metrics->errors_by_read_type.has_value()) {
      accuracy_metrics->errors_by_read_type->WriteTsv(
          PrefixedOutputPath(accuracy_dir, output_prefix, kDefaultErrorByReadType), tsv_metadata);
    }
    if (accuracy_metrics->hp_stats.has_value()) {
      accuracy_metrics->hp_stats->WriteTsv(PrefixedOutputPath(accuracy_dir, output_prefix, kDefaultHpErrors),
                                           tsv_metadata);
    }
  }
  if (coverage_metrics.has_value()) {
    create_directories(absolute_out_dir / kDefaultCoverageDir);
    coverage_metrics->coverage_histograms.WriteTsv(
        PrefixedOutputPath(coverage_dir, output_prefix, kDefaultCoverageHistograms), tsv_metadata);
    coverage_metrics->coverage_histograms.WriteCoverageStatsTsv(
        PrefixedOutputPath(coverage_dir, output_prefix, kDefaultCoverageStats), tsv_metadata);
    coverage_metrics->coverage_histograms.WriteCoverageDistributionSummaryTsv(
        PrefixedOutputPath(coverage_dir, output_prefix, kDefaultCoverageDistributionSummary), tsv_metadata);
    if (coverage_metrics->hp_coverage_histograms.has_value()) {
      coverage_metrics->hp_coverage_histograms->WriteTsv(
          PrefixedOutputPath(coverage_dir, output_prefix, kDefaultHpCoverageHistograms), tsv_metadata);
      coverage_metrics->hp_coverage_histograms->WriteCoverageStatsTsv(
          PrefixedOutputPath(coverage_dir, output_prefix, kDefaultHpCoverageStats), tsv_metadata);
      coverage_metrics->hp_coverage_histograms->WriteCoverageDistributionSummaryTsv(
          PrefixedOutputPath(coverage_dir, output_prefix, kDefaultHpCoverageDistributionSummary), tsv_metadata);
    }
    if (coverage_metrics->coverage_uniformity_metrics.has_value()) {
      coverage_metrics->coverage_uniformity_metrics->WriteMeanCoverageHistogram(
          PrefixedOutputPath(coverage_dir, output_prefix, kDefaultMeanCoverageHistogram), tsv_metadata);
      coverage_metrics->coverage_uniformity_metrics->WriteCoverageUniformitySummaryTsv(
          PrefixedOutputPath(coverage_dir, output_prefix, kDefaultCoverageUniformitySummary),
          tsv_metadata,
          coverage_metrics->coverage_histograms.filtered_coverage_histogram);
    }
  }
  if (gc_bias_metrics.has_value()) {
    const auto gc_bias_dir = output_dir / kDefaultGcBiasDir;
    create_directories(absolute_out_dir / kDefaultGcBiasDir);
    gc_bias_metrics->WriteTsv(PrefixedOutputPath(gc_bias_dir, output_prefix, kDefaultGcBiasMetrics), tsv_metadata);
    gc_bias_metrics->WriteSummaryTsv(PrefixedOutputPath(gc_bias_dir, output_prefix, kDefaultGcBiasSummaryMetrics),
                                     tsv_metadata);
  }
  if (read_metrics.has_value()) {
    create_directories(absolute_out_dir / kDefaultReadMetricsDir);
    if (read_metrics->read_counts_by_cluster_size.has_value()) {
      read_metrics->read_counts_by_cluster_size->WriteTsv(
          PrefixedOutputPath(read_metrics_dir, output_prefix, kDefaultReadCountsByClusterSize), tsv_metadata);
    }
    if (read_metrics->te_read_metrics.has_value()) {
      read_metrics->te_read_metrics->WriteTsv(
          PrefixedOutputPath(read_metrics_dir, output_prefix, kDefaultTeReadMetrics), tsv_metadata);
    }
    read_metrics->read_metrics_summary.WriteTsv(
        PrefixedOutputPath(read_metrics_dir, output_prefix, kDefaultReadMetricsSummary), tsv_metadata);
    read_metrics->read_length_histograms.WriteTsv(
        PrefixedOutputPath(read_metrics_dir, output_prefix, kDefaultReadLengthHistograms), tsv_metadata);
    read_metrics->read_length_histograms.WriteSummaryTsv(
        PrefixedOutputPath(read_metrics_dir, output_prefix, kDefaultReadLengthSummary), tsv_metadata);
    read_metrics->read_length_histograms.WriteExcludingSoftClipsTsv(
        PrefixedOutputPath(read_metrics_dir, output_prefix, kDefaultReadLengthExcludingSoftClipsHistograms),
        tsv_metadata);
    read_metrics->read_length_histograms.WriteExcludingSoftClipsSummaryTsv(
        PrefixedOutputPath(read_metrics_dir, output_prefix, kDefaultReadLengthExcludingSoftClipsSummary), tsv_metadata);
  }
}

Metrics& Metrics::operator+=(const Metrics& other) {
  if (accuracy_metrics && other.accuracy_metrics) {
    *accuracy_metrics += *other.accuracy_metrics;
  }
  if (coverage_metrics && other.coverage_metrics) {
    *coverage_metrics += *other.coverage_metrics;
  }
  if (read_metrics && other.read_metrics) {
    *read_metrics += *other.read_metrics;
  }
  if (gc_bias_metrics && other.gc_bias_metrics) {
    *gc_bias_metrics += *other.gc_bias_metrics;
  }
  return *this;
}

}  // namespace xoos::alignment_metrics
