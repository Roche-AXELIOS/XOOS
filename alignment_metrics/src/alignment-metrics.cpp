#include "alignment-metrics.h"

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <htslib/sam.h>

#include <csv.hpp>

#include <taskflow/taskflow.hpp>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>
#include <xoos/util/file-functions.h>

#include "alignment-metrics-options.h"
#include "core/region-lookup.h"
#include "io/alignment-reader.h"
#include "metadata/dataset-metadata.h"
#include "metrics-worker.h"

namespace xoos::alignment_metrics {

/**
 * Constructor for CoverageMetrics class.
 *
 * During initialization, we detect the dataset type from the BAM file
 * and initialize the object to keep the final metrics. The dataset metadata
 * is used to determine what metrics to output and how to format them.
 **/
AlignmentMetrics::AlignmentMetrics(const AlignmentMetricsOptions& opts)
    : options(opts),
      dataset_metadata{GetDatasetMetadataFromBam(opts.bam_input)},
      final_metrics(options, dataset_metadata) {
}

// Create the output directory if it doesn't exist and confirm that it is a directory.
void CreateOutputDir(const AlignmentMetricsOptions& options) {
  // if it exists, we check to make sure it is writable
  file::CreateWritableDirectory(std::filesystem::absolute(options.out_dir));
}

// Merge the metrics from each worker thread into a single set of metrics on this local thread.
void AlignmentMetrics::MergeMetrics(const std::vector<std::shared_ptr<Metrics>>& per_thread_metrics) {
  for (const auto& metrics : per_thread_metrics) {
    final_metrics += *metrics;
  }
}

/**
 * Loop over regions and calculate alignment metrics for each region in parallel.
 */
void AlignmentMetrics::CalculateMetrics() {
  CreateOutputDir(options);

  // We do not partition regions when TE metrics are enabled to ensure correct calculation
  // of the total number of target regions and coverage uniformity metrics.
  const RegionLookupTable region_lookup_table(options.bam_input,
                                              options.bed_input,
                                              options.reference_input,
                                              options.region_size,
                                              options.calculate_hp_metrics,
                                              options.enable_te_metrics);

  // Check the largest super region size and warn if it exceeds recommended size
  for (const auto& super_region : region_lookup_table.GetSuperRegions()) {
    const auto super_region_size = static_cast<size_t>(super_region.end - super_region.start);
    if (super_region_size > kRecommendedMaxRegionSize) {
      if (options.enable_te_metrics) {
        XOOS_LOG_WARN_ONCE(
            "Region {}:{}-{} is larger than the recommended maximum size of {} bases. "
            "This will lead to high memory usage because --enable-te-metrics disables region partitioning. "
            "Check if the provided BED file is an actual target panel and set --region-size to a smaller value.",
            super_region.chromosome,
            super_region.start,
            super_region.end,
            kRecommendedMaxRegionSize);
      } else {
        XOOS_LOG_WARN_ONCE(
            "Region {}:{}-{} is larger than the recommended maximum region size of {} bases. "
            "This will lead to high memory usage. Consider setting --region-size to a smaller value.",
            super_region.chromosome,
            super_region.start,
            super_region.end,
            kRecommendedMaxRegionSize);
      }
      break;
    }
  }

  // Create NewMetricsWorker for each thread to process the regions in parallel
  std::vector<MetricsWorker> metrics_workers;
  metrics_workers.reserve(options.threads);
  for (u32 i = 0; i < options.threads; ++i) {
    metrics_workers.emplace_back(options, dataset_metadata, region_lookup_table);
  }

  Logging::Info("Calculating metrics for BAM file: {}", options.bam_input.string());
  Logging::Info("Dataset is {} dataset",
                dataset_metadata.has_cluster_info ? "an intermolecular consensus"
                                                  : (dataset_metadata.is_duplex_dataset ? "a duplex" : "a simplex"));
  Logging::Info("Running with {} threads", options.threads);

  tf::Executor executor(options.threads);
  if (executor.this_worker_id() != -1) {
    throw error::Error("Worker thread is the master thread, exiting");
  }

  tf::Taskflow taskflow;
  // For each super region with coverage, schedule a task to process the region and calculate metrics for the region
  for (size_t super_region_index = 0; super_region_index < region_lookup_table.GetSuperRegions().size();
       ++super_region_index) {
    auto super_region = region_lookup_table.GetSuperRegions().at(super_region_index);
    taskflow.emplace([&metrics_workers, super_region, super_region_index, &executor]() {
      metrics_workers.at(executor.this_worker_id()).ProcessRegion(super_region, super_region_index);
    });
  }
  // Schedule a tasks to calculate pileup data or add regions with no reads to the zero-coverage histogram bin.
  for (const auto& [chromosome, regions] : region_lookup_table.GetRegionsWithoutCoverage()) {
    taskflow.emplace([&executor, &metrics_workers, &regions]() {
      for (const auto& reg : regions) {
        metrics_workers.at(ToUnsigned(executor.this_worker_id())).ProcessNoCoverageRegion(reg);
      }
    });
  }
  // Add a separate task for calculating read-level metrics for unplaced unmapped reads
  taskflow.emplace([this]() {
    if (options.metric_types.has_read_metrics) {
      const auto reader = OpenAlignmentFile(options.bam_input);
      // Unmapped reads only pass filter if the user specified min-mapq 0 and did not exclude unmapped reads by flags
      const bool unmapped_reads_passed_filter = options.min_mapq == 0 && (options.exclude_flags & BAM_FUNMAP) == 0;
      final_metrics.read_metrics->UpdateUnmappedReadMetrics(
          reader, unmapped_reads_passed_filter, options.trim_leading_bases + options.trim_trailing_bases);
    }
  });
  // Run the taskflow to process all regions in parallel and wait for completion.
  executor.run(taskflow).get();

  // Merge the metrics from each worker into a single set of metrics on this local thread.
  std::vector<std::shared_ptr<Metrics>> per_thread_metrics;
  per_thread_metrics.reserve(options.threads);
  for (const auto& worker : metrics_workers) {
    per_thread_metrics.emplace_back(worker.GetMetrics());
  }
  MergeMetrics(per_thread_metrics);

  // Write out the metrics to a file
  Logging::Info("Writing metrics to output directory: {}", options.out_dir.string());
  final_metrics.WriteMetrics(options.out_dir, options.comments);

  Logging::Info("Finished writing metrics");
}

}  // namespace xoos::alignment_metrics
