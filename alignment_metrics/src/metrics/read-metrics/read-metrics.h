#pragma once

#include <optional>

#include <xoos/io/alignment-reader.h>
#include <xoos/types/int.h>

#include "alignment-metrics-options.h"
#include "core/alignment.h"
#include "core/region-lookup.h"
#include "core/supplementary-alignment.h"
#include "metadata/alignment-metadata.h"
#include "metadata/dataset-metadata.h"
#include "metrics/read-metrics/read-counts-by-cluster-size.h"
#include "metrics/read-metrics/read-length-histograms.h"
#include "metrics/read-metrics/read-metrics-summary.h"
#include "metrics/read-metrics/read-metrics-te.h"

namespace xoos::alignment_metrics {

class Metrics;

/**
 * @brief Container for all read-level and alignment statistics.
 *
 * Aggregates read counts, alignment flags, mapping statistics, and read length distributions to provide comprehensive
 * summaries of sequencing run quality and read characteristics.
 */
class ReadMetrics {
  friend class Metrics;  // Allow Metrics to access private members for updating TE and unmapped metrics
 public:
  ReadMetrics() = default;
  explicit ReadMetrics(const AlignmentMetricsOptions& options, const DatasetMetadata& dataset_metadata);

  // ReadMetrics aggregates sub-metric holders that are updated in place as alignments are processed;
  // these fields are intentionally public.
  // NOLINTBEGIN(misc-non-private-member-variables-in-classes,cppcoreguidelines-non-private-member-variables-in-classes)
  ReadMetricsSummary read_metrics_summary;
  ReadLengthHistograms read_length_histograms;
  std::optional<ReadCountsByClusterSize> read_counts_by_cluster_size;
  std::optional<TargetEnrichmentReadMetrics> te_read_metrics;
  // NOLINTEND(misc-non-private-member-variables-in-classes,cppcoreguidelines-non-private-member-variables-in-classes)

  /**
   * Update read metrics using the given alignment. We assume that the calling code has already checked that the
   * alignment is not double counted
   *
   * This updates various read counts, read length distributions, and base counts (including
   * `aligned_bases` and `soft_clipped_bases`).
   *
   * @param alignment Alignment to process; used to extract BAM flags and read length
   * @param alignment_metadata Metadata for the alignment containing information such as full/partial status
   * @param passed_filter Whether the read passed the minimum mapping quality and BAM flag filters
   * @param regions_of_interest When built from a BED file, supplementary alignments only contribute to
   *        `aligned_bases` when their reference span overlaps a region of interest; when no BED file
   *        was supplied every supplementary alignment contributes
   */
  void AddRead(const Alignment& alignment,
               const AlignmentMetadata& alignment_metadata,
               bool passed_filter,
               const RegionLookupTable& regions_of_interest);
  /**
   * Update read metrics using the given alignment. We assume that the calling code has already checked that the
   * alignment is not double counted
   *
   * This updates various read counts, read length distributions, and base counts (including
   * `aligned_bases` and `soft_clipped_bases`).
   *
   * @param alignment Alignment to process; used to extract BAM flags and read length
   * @param alignment_metadata Metadata for the alignment containing information such as full/partial status
   * @param sa_tag_filter Filter criteria for supplementary alignments (min mapQ, exclude flags)
   * @param passed_filter Whether the read passed the minimum mapping quality and BAM flag filters
   * @param regions_of_interest When built from a BED file, supplementary alignments only contribute to
   *        `aligned_bases` when their reference span overlaps a region of interest; when no BED file
   *        was supplied every supplementary alignment contributes
   */
  void AddRead(const Alignment& alignment,
               const AlignmentMetadata& alignment_metadata,
               const SaTagFilter& sa_tag_filter,
               bool passed_filter,
               const RegionLookupTable& regions_of_interest);
  /**
   * Updates the unmapped read metrics in the final metrics object
   * by querying the BAM index for the number of unplaced unmapped reads (no coordinate).
   * Placed unmapped reads (reads with the unmapped flag 0x4 but has a coordinate) are
   * counted separately when querying the BAM file by regions.
   *
   * The metric object is updated in place. This function should only be called once
   * for the final metrics object (after merging per-thread metric objects) to avoid
   * double counting.
   */
  void UpdateUnmappedReadMetrics(const io::AlignmentReader& alignment_reader, bool passed_filter, u16 trim);
  ReadMetrics& operator+=(const ReadMetrics& other);

 private:
  // Per-record base counts derived from a primary alignment, aggregated into read_metrics_summary.
  struct ReadBaseCounts {
    u32 read_length{};
    u32 aligned_bases{};
    u32 aligned_bases_passing_filter{};
    u32 unique_aligned_bases{};
    u32 unique_aligned_bases_passing_filter{};
    u32 soft_clipped_bases{};
    u32 unmapped_soft_clipped_bases{};
    u32 gc_bases{};
  };

  // Records counters that apply to every alignment record (total, supplementary, secondary, duplicate).
  void RecordAlignmentRecordCounts(u32 flag);
  // Aggregates primary-alignment base counts into the summary (pre-filter and, when passed, post-filter).
  void RecordBaseCounts(const ReadBaseCounts& base_counts, bool passed_filter);
  // Records read length histograms for a mapped primary read (all + full/partial variants).
  void RecordReadLengthHistograms(u32 read_length,
                                  u32 read_length_without_softclips,
                                  ReadType read_type,
                                  bool passed_filter);
  // Records forward/reverse read counts for a mapped primary read.
  void RecordReadCountsByStrand(ReadStrand read_strand, bool passed_filter);
  // Records full/partial read counts for a primary read.
  void RecordReadCountsByType(ReadType read_type, bool passed_filter);
  // Records per-cluster read counts for a primary read (post-consensus datasets only).
  void RecordReadCountsByClusterType(const AlignmentMetadata& alignment_metadata, bool passed_filter);
  // Updates target-enrichment on-target counters from the current aggregate summary.
  void UpdateTargetEnrichmentReadMetrics();

  bool _unmapped_metrics_updated = false;
};

}  // namespace xoos::alignment_metrics
