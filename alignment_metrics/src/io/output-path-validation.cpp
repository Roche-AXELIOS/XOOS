#include "io/output-path-validation.h"

#include <filesystem>
#include <string>

#include <xoos/error/error.h>
#include <xoos/util/file-functions.h>

#include "io/metrics-filenames.h"

namespace xoos::alignment_metrics {

using file::PrefixedOutputPath;

vec<fs::path> CollectOutputPaths(const AlignmentMetricsOptions& options,
                                 const DatasetMetadata& metadata,
                                 const fs::path& out_dir,
                                 const std::string& output_prefix) {
  vec<fs::path> paths;
  if (options.metric_types.has_accuracy_metrics) {
    const auto dir = out_dir / kDefaultAccuracyDir;
    paths.push_back(PrefixedOutputPath(dir, output_prefix, kDefaultErrorBySubstitutionType));
    paths.push_back(PrefixedOutputPath(dir, output_prefix, kDefaultBaseLevelAccuracySummary));
    paths.push_back(PrefixedOutputPath(dir, output_prefix, kDefaultQscoreStats));
    if (metadata.has_cluster_info) {
      paths.push_back(PrefixedOutputPath(dir, output_prefix, kDefaultErrorByClusterSize));
    }
    if (metadata.has_read_type_info || metadata.has_strand_info) {
      paths.push_back(PrefixedOutputPath(dir, output_prefix, kDefaultErrorByReadType));
    }
    if (options.calculate_hp_metrics) {
      paths.push_back(PrefixedOutputPath(dir, output_prefix, kDefaultHpErrors));
    }
  }
  if (options.metric_types.has_coverage_metrics) {
    const auto dir = out_dir / kDefaultCoverageDir;
    paths.push_back(PrefixedOutputPath(dir, output_prefix, kDefaultCoverageHistograms));
    paths.push_back(PrefixedOutputPath(dir, output_prefix, kDefaultCoverageStats));
    paths.push_back(PrefixedOutputPath(dir, output_prefix, kDefaultCoverageDistributionSummary));
    if (options.calculate_hp_metrics) {
      paths.push_back(PrefixedOutputPath(dir, output_prefix, kDefaultHpCoverageHistograms));
      paths.push_back(PrefixedOutputPath(dir, output_prefix, kDefaultHpCoverageStats));
      paths.push_back(PrefixedOutputPath(dir, output_prefix, kDefaultHpCoverageDistributionSummary));
    }
    if (options.enable_te_metrics) {
      paths.push_back(PrefixedOutputPath(dir, output_prefix, kDefaultMeanCoverageHistogram));
      paths.push_back(PrefixedOutputPath(dir, output_prefix, kDefaultCoverageUniformitySummary));
    }
  }
  if (options.metric_types.has_gc_bias_metrics) {
    const auto dir = out_dir / kDefaultGcBiasDir;
    paths.push_back(PrefixedOutputPath(dir, output_prefix, kDefaultGcBiasMetrics));
    paths.push_back(PrefixedOutputPath(dir, output_prefix, kDefaultGcBiasSummaryMetrics));
  }
  if (options.metric_types.has_read_metrics) {
    const auto dir = out_dir / kDefaultReadMetricsDir;
    paths.push_back(PrefixedOutputPath(dir, output_prefix, kDefaultReadMetricsSummary));
    paths.push_back(PrefixedOutputPath(dir, output_prefix, kDefaultReadLengthHistograms));
    paths.push_back(PrefixedOutputPath(dir, output_prefix, kDefaultReadLengthSummary));
    paths.push_back(PrefixedOutputPath(dir, output_prefix, kDefaultReadLengthExcludingSoftClipsHistograms));
    paths.push_back(PrefixedOutputPath(dir, output_prefix, kDefaultReadLengthExcludingSoftClipsSummary));
    if (metadata.has_cluster_info) {
      paths.push_back(PrefixedOutputPath(dir, output_prefix, kDefaultReadCountsByClusterSize));
    }
    if (options.enable_te_metrics) {
      paths.push_back(PrefixedOutputPath(dir, output_prefix, kDefaultTeReadMetrics));
    }
  }
  return paths;
}

void ValidateOutputFilesDoNotExist(const AlignmentMetricsOptions& options,
                                   const DatasetMetadata& metadata,
                                   const RegionLabelsSet& labels,
                                   const std::string& output_prefix) {
  vec<fs::path> all_paths;
  if (labels.Empty()) {
    all_paths = CollectOutputPaths(options, metadata, options.out_dir, output_prefix);
  } else {
    for (const auto& label : labels.Names()) {
      auto label_paths = CollectOutputPaths(options, metadata, options.out_dir / label, output_prefix);
      all_paths.insert(all_paths.end(), label_paths.begin(), label_paths.end());
    }
  }
  vec<fs::path> existing;
  for (const auto& path : all_paths) {
    if (fs::exists(path)) {
      existing.push_back(path);
    }
  }
  if (!existing.empty()) {
    std::string message = "Output file(s) already exist. Use --overwrite to overwrite them:\n";
    for (const auto& path : existing) {
      message += "  " + path.string() + "\n";
    }
    throw error::Error(message);
  }
}

}  // namespace xoos::alignment_metrics
