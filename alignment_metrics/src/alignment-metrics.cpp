#include "alignment-metrics.h"

#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include <htslib/sam.h>

#include <csv.hpp>

#include <taskflow/taskflow.hpp>

#include <xoos/enum/enum-util.h>
#include <xoos/error/error.h>
#include <xoos/io/alignment-reader.h>
#include <xoos/io/bam/read-group.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>
#include <xoos/util/file-functions.h>

#include "alignment-metrics-options.h"
#include "core/region-labels-set.h"
#include "core/region-lookup.h"
#include "io/output-path-validation.h"
#include "metadata/dataset-metadata.h"
#include "metrics-worker.h"
#include "metrics/stratified-metrics.h"

namespace xoos::alignment_metrics {

namespace {

/**
 * @brief Validates the size of a super region against the recommended maximum size.
 *
 * Logs a warning if the super region exceeds the recommended maximum size. The warning message
 * differs depending on whether target enrichment metrics are enabled.
 *
 * @param super_region The super region to validate
 * @param options The alignment metrics options containing the target enrichment setting
 * @return True if the super region size is within the recommended limit, false otherwise
 */
bool ValidateSuperRegionSize(const SuperRegion& super_region, const AlignmentMetricsOptions& options) {
  const auto super_region_size = static_cast<size_t>(super_region.end - super_region.start);
  if (super_region_size > kRecommendedMaxRegionSize) {
    if (options.enable_te_metrics) {
      XOOS_LOG_WARN_ONCE(
          "Region {}:{}-{} is larger than the recommended maximum size of {} bases. "
          "This will lead to high memory usage because --enable-target-enrichment-metrics disables region "
          "partitioning. "
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
  }
  return super_region_size <= kRecommendedMaxRegionSize;
}

/**
 * @brief Validates that the output would not be empty based on input and super region status.
 *
 * Logs a warning if the input BAM file contains no alignment records or if no contigs with the requested
 * regions have any mapped reads. Indicates that empty metrics will be written with zero values and
 * metrics that cannot be computed are reported as NA.
 *
 * @param is_input_empty True if the input BAM file contains no alignment records
 * @param is_super_region_empty True if no contigs with the requested regions have any mapped reads
 * @return True if the output is non-empty, false otherwise
 */
bool ValidateNonEmptyOutput(const bool is_input_empty, const bool is_super_region_empty) {
  bool passed = true;
  if (is_input_empty) {
    Logging::Warn(
        "Input BAM file contains no alignment records. "
        "Writing empty metrics with zero values. "
        "Metrics that cannot be computed are reported as NA.");
    passed = false;
  } else if (is_super_region_empty) {
    Logging::Warn(
        "No contigs containing the requested regions have any mapped reads. "
        "Writing empty metrics with zero values. "
        "Metrics that cannot be computed are reported as NA.");
    passed = false;
  }
  return passed;
}

/**
 * @brief Validates the base quality and base type settings against the dataset metadata.
 *
 * Logs warnings if the specified minimum base quality or base type settings are inconsistent with the
 * characteristics of the dataset (e.g., duplex, simplex, or intermolecular consensus). This helps
 * prevent the accidental inclusion of low-confidence bases or exclusion of a significant portion of
 * the dataset when the base quality or base type settings are not appropriately configured.
 *
 * @param options The alignment metrics options containing base quality and base type settings
 * @param dataset_metadata Metadata describing the dataset, including cluster and duplex information
 */
void ValidateBaseQualitySettings(const AlignmentMetricsOptions& options, const DatasetMetadata& dataset_metadata) {
  constexpr u16 kDuplexBaseQualityWarningThreshold = 25;
  constexpr u16 kSimplexBaseQualityWarningThreshold = 0;
  constexpr u16 kPostConsensusBaseQualityWarningThreshold = 23;
  constexpr auto kDuplexRequiredBaseType = yc_decode::BaseType::kConcordant;
  constexpr auto kSimplexBaseTypeWarningThreshold = yc_decode::BaseType::kSimplex;
  constexpr auto kPostConsensusBaseTypeWarningThreshold = yc_decode::BaseType::kSimplex;

  if (!options.metric_types.has_accuracy_metrics && !options.metric_types.has_coverage_metrics) {
    return;
  }

  if (dataset_metadata.has_cluster_info) {
    if (options.min_baseq < kPostConsensusBaseQualityWarningThreshold) {
      Logging::Warn(
          "Dataset is intermolecular consensus but --min-base-quality is set to {}. This will cause low-confidence "
          "consensus bases to be included in coverage and accuracy calculation",
          options.min_baseq);
    }
    if (options.min_base_type > kPostConsensusBaseTypeWarningThreshold) {
      Logging::Warn(
          "Dataset is intermolecular consensus but --min-base-type is set to {}. This will cause all bases to be "
          "excluded from coverage and accuracy calculation.",
          enum_util::FormatEnumName(options.min_base_type));
    }
    return;
  }

  if (dataset_metadata.is_duplex_dataset) {
    if (options.min_baseq < kDuplexBaseQualityWarningThreshold) {
      Logging::Warn(
          "Dataset is duplex but --min-base-quality is set to {}. This will cause non-concordant duplex bases to be "
          "included in coverage and accuracy calculation",
          options.min_baseq);
    }
    if (options.min_base_type != kDuplexRequiredBaseType) {
      Logging::Warn(
          "Dataset is duplex but --min-base-type is set to {}. This will cause non-concordant duplex bases to be "
          "included in coverage and accuracy calculation",
          enum_util::FormatEnumName(options.min_base_type));
    }
    return;
  }

  if (options.min_baseq > kSimplexBaseQualityWarningThreshold) {
    Logging::Warn(
        "Dataset is simplex but --min-base-quality is set to {}. This may cause some bases to be excluded from "
        "coverage and accuracy calculation",
        options.min_baseq);
  }
  if (options.min_base_type > kSimplexBaseTypeWarningThreshold) {
    Logging::Warn(
        "Dataset is simplex but --min-base-type is set to {}. This will cause all bases to be excluded from coverage "
        "and accuracy calculation.",
        enum_util::FormatEnumName(options.min_base_type));
  }
}

}  // namespace

/**
 * Constructor for CoverageMetrics class.
 *
 * During initialization, we detect the dataset type from the BAM file
 * and initialize the object to keep the final metrics. The dataset metadata
 * is used to determine what metrics to output and how to format them.
 **/
AlignmentMetrics::AlignmentMetrics(const AlignmentMetricsOptions& opts)
    : options(opts), dataset_metadata{GetDatasetMetadataFromBam(opts.bam_input)} {
  // Initialize with a non-stratified (single-container) metrics set. CalculateMetrics re-initializes
  // this once the region labels set is known from the RegionLookupTable.
  final_metrics.emplace(options, dataset_metadata, RegionLabelsSet{});
}

// Create the output directory if it doesn't exist and confirm that it is a directory.
void CreateOutputDir(const AlignmentMetricsOptions& options) {
  // if it exists, we check to make sure it is writable
  file::CreateWritableDirectory(std::filesystem::absolute(options.out_dir));
}

// Merge the per-label metrics from each worker thread into a single set of metrics on this local thread.
void AlignmentMetrics::MergeMetrics(const std::vector<std::shared_ptr<StratifiedMetrics>>& per_thread_metrics) {
  for (const auto& metrics : per_thread_metrics) {
    *final_metrics += *metrics;
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
                                              options.enable_te_metrics,
                                              options.stratify_by_label);

  // Now that the region labels set is known, re-initialize the final per-label metrics container. In
  // non-stratified mode (no labels) this holds a single container written to the top-level output dir.
  final_metrics.emplace(options, dataset_metadata, region_lookup_table.GetLabels());

  if (!options.overwrite) {
    ValidateOutputFilesDoNotExist(options, dataset_metadata, region_lookup_table.GetLabels(), options.output_prefix);
  }

  // Check the largest super region size and warn if it exceeds recommended size
  for (const auto& super_region : region_lookup_table.GetSuperRegions()) {
    const bool is_super_region_size_valid = ValidateSuperRegionSize(super_region, options);
    if (!is_super_region_size_valid) {
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
  std::string_view dataset_description = "a simplex";
  if (dataset_metadata.has_cluster_info) {
    dataset_description = "an intermolecular consensus";
  } else if (dataset_metadata.is_duplex_dataset) {
    dataset_description = "a duplex";
  }
  Logging::Info("Dataset is {} dataset", dataset_description);
  ValidateBaseQualitySettings(options, dataset_metadata);
  Logging::Info("Running with {} threads", options.threads);

  ValidateNonEmptyOutput(dataset_metadata.is_empty, region_lookup_table.GetSuperRegions().empty());

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
      size_t worker_id = executor.this_worker_id() >= 0 ? ToUnsigned(executor.this_worker_id()) : 0;
      metrics_workers.at(worker_id).ProcessRegion(super_region, super_region_index);
    });
  }
  // Schedule a tasks to calculate pileup data or add regions with no reads to the zero-coverage histogram bin.
  for (const auto& [chromosome, regions] : region_lookup_table.GetRegionsWithoutCoverage()) {
    taskflow.emplace([&executor, &metrics_workers, &regions, chromosome_str = chromosome]() {
      size_t worker_id = executor.this_worker_id() >= 0 ? ToUnsigned(executor.this_worker_id()) : 0;
      for (const auto& region : regions) {
        metrics_workers.at(worker_id).ProcessNoCoverageRegion(chromosome_str, region);
      }
    });
  }
  // Add a separate task for calculating read-level metrics for unplaced unmapped reads. These reads
  // belong to no region, so they are counted in every label container (mirroring an independent
  // single-label run per label).
  taskflow.emplace([this]() {
    if (options.metric_types.has_read_metrics) {
      // Unmapped reads only pass filter if the user specified min-mapq 0 and did not exclude unmapped reads by flags
      const bool unmapped_reads_passed_filter = options.min_mapq == 0 && (options.exclude_flags & BAM_FUNMAP) == 0;
      // Unplaced unmapped reads belong to no region, so their contribution is identical for every
      // label container. Scan the BAM once into a scratch container, then fold that result into every
      // container; this avoids re-scanning the BAM once per label (up to kMaxLabels times).
      Metrics unmapped_scratch(options, dataset_metadata);
      const auto reader = io::OpenAlignmentReader(options.bam_input);
      unmapped_scratch.read_metrics->UpdateUnmappedReadMetrics(
          reader, unmapped_reads_passed_filter, options.trim_leading_bases + options.trim_trailing_bases);
      final_metrics->ForAllContainers(
          [&unmapped_scratch](Metrics& target) { *target.read_metrics += *unmapped_scratch.read_metrics; });
    }
  });
  // Run the taskflow to process all regions in parallel and wait for completion.
  executor.run(taskflow).get();

  // Merge the metrics from each worker into a single set of metrics on this local thread.
  std::vector<std::shared_ptr<StratifiedMetrics>> per_thread_metrics;
  per_thread_metrics.reserve(options.threads);
  for (const auto& worker : metrics_workers) {
    per_thread_metrics.emplace_back(worker.GetStratifiedMetrics());
  }
  MergeMetrics(per_thread_metrics);

  // Write out the metrics to a file. In stratified mode this creates one subdirectory per label.
  Logging::Info("Writing metrics to output directory: {}", options.out_dir.string());
  const io::TsvMetadata tsv_metadata{.command_line_info = options.command_line_info,
                                     .read_group = io::GetFirstReadGroupFromBam(options.bam_input)};
  final_metrics->WriteMetrics(options.out_dir, tsv_metadata, options.output_prefix);

  Logging::Info("Finished writing metrics");
}

}  // namespace xoos::alignment_metrics
