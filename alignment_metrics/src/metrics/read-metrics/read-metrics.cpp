#include "metrics/read-metrics/read-metrics.h"

#include <cstddef>
#include <optional>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <xoos/histogram/histogram.h>
#include <xoos/io/alignment-reader.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/io/htslib-util/read-util.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/util/math.h>

#include "alignment-metrics-options.h"
#include "core/alignment.h"
#include "core/region-lookup.h"
#include "core/super-region.h"
#include "core/supplementary-alignment.h"
#include "metadata/alignment-metadata.h"
#include "metadata/dataset-metadata.h"

namespace xoos::alignment_metrics {

ReadMetrics::ReadMetrics(const AlignmentMetricsOptions& options, const DatasetMetadata& dataset_metadata)
    : read_metrics_summary(dataset_metadata),
      read_length_histograms(options.max_read_length_bin, options.summary_stats_percentiles, dataset_metadata),
      read_counts_by_cluster_size(dataset_metadata.has_cluster_info
                                      ? std::make_optional(ReadCountsByClusterSize(options.max_cluster_size_bin))
                                      : std::nullopt),
      te_read_metrics(options.enable_te_metrics ? std::make_optional(TargetEnrichmentReadMetrics(options.bam_input))
                                                : std::nullopt) {
}

ReadMetrics& ReadMetrics::operator+=(const ReadMetrics& other) {
  read_metrics_summary += other.read_metrics_summary;
  read_length_histograms += other.read_length_histograms;
  if (read_counts_by_cluster_size && other.read_counts_by_cluster_size) {
    read_counts_by_cluster_size->operator+=(*other.read_counts_by_cluster_size);
  }
  if (te_read_metrics && other.te_read_metrics) {
    te_read_metrics->operator+=(*other.te_read_metrics);
  }
  // If the other metrics object has already updated unmapped metrics,
  // we don't want to double count them.
  _unmapped_metrics_updated = _unmapped_metrics_updated || other._unmapped_metrics_updated;
  return *this;
}

/**
 * @brief Records read-level statistics for a single alignment.
 *
 * Updates read metrics counters based on alignment flags (mapped/unmapped, forward/reverse,
 * supplementary/secondary, duplicate), read type (full/partial), and cluster characteristics.
 * Increments read length histograms for reads passing quality filters. For post-consensus data,
 * stratifies reads by cluster size, strand composition, and completeness. This method is called
 * for every alignment read processed.
 *
 * @param alignment The alignment object containing BAM record and metadata
 * @param alignment_metadata Parsed metadata including read type, strand, and cluster info
 * @param passed_filter Whether the alignment passed quality filters
 */
void ReadMetrics::AddRead(const Alignment& alignment,
                          const AlignmentMetadata& alignment_metadata,
                          const bool passed_filter,
                          const RegionLookupTable& regions_of_interest) {
  // Call the more general AddRead with a default SaTagFilter that does not filter any supplementary alignments
  AddRead(alignment, alignment_metadata, SaTagFilter{}, passed_filter, regions_of_interest);
}

/**
 * @brief Records read-level statistics for a single alignment.
 *
 * Updates read metrics counters based on alignment flags (mapped/unmapped, forward/reverse,
 * supplementary/secondary, duplicate), read type (full/partial), and cluster characteristics.
 * Increments read length histograms for reads passing quality filters. For post-consensus data,
 * stratifies reads by cluster size, strand composition, and completeness. This method is called
 * for every alignment read processed.
 *
 * @param alignment The alignment object containing BAM record and metadata
 * @param alignment_metadata Parsed metadata including read type, strand, and cluster info
 * @param sa_tag_filter Filter criteria for supplementary alignments (min mapQ, exclude flags)
 * @param passed_filter Whether the alignment passed quality filters
 */
void ReadMetrics::AddRead(const Alignment& alignment,
                          const AlignmentMetadata& alignment_metadata,
                          const SaTagFilter& sa_tag_filter,
                          const bool passed_filter,
                          const RegionLookupTable& regions_of_interest) {
  const u32 read_length = alignment.read_length;
  // If a read has no length (e.g. after trimming), it is not counted
  if (read_length == 0) {
    return;
  }
  const u32 flag = alignment.alignment_ptr->core.flag;
  const bool mapped = (flag & u32{BAM_FUNMAP}) == 0;
  const bool primary = (flag & (u32{BAM_FSUPPLEMENTARY} | u32{BAM_FSECONDARY})) == 0;

  RecordAlignmentRecordCounts(flag);

  // Only primary alignments contribute read-level counts; secondary and supplementary records are
  // handled entirely by RecordAlignmentRecordCounts above.
  if (primary) {
    // aligned_bases spans the primary and its supplementary alignments. The primary's contribution is
    // its own aligned base count; the supplementary contribution is derived from the SA tag and
    // restricted to supplementary alignments overlapping a region of interest. A base aligned in both
    // the primary and a supplementary is counted in each, which is why aligned_bases has no percentage.
    const u32 primary_aligned_bases = alignment.CountAlignedBases();
    const u32 read_length_without_softclips = alignment.read_length_without_softclips;
    const u32 soft_clipped_bases = read_length - read_length_without_softclips;
    // Soft-clip-derived counts only make sense for a mapped primary alignment, which carries the full
    // read and the SA tag describing its supplementary alignments.
    const SoftClipBaseCounts softclip_counts =
        mapped ? alignment.CountSoftClippedBasesWithSupplementary(sa_tag_filter, regions_of_interest)
               : SoftClipBaseCounts{};
    const ReadBaseCounts base_counts{
        .read_length = read_length,
        .aligned_bases = primary_aligned_bases + softclip_counts.supplementary_aligned_bases,
        .aligned_bases_passing_filter =
            passed_filter ? primary_aligned_bases + softclip_counts.supplementary_aligned_bases_passing_filter : 0,
        // unique_aligned_bases = aligned bases in the primary alignment + soft-clipped primary bases that
        // are aligned in some supplementary alignment (i.e. soft-clipped bases minus those aligned nowhere).
        .unique_aligned_bases = primary_aligned_bases + softclip_counts.aligned_soft_clipped_bases,
        .unique_aligned_bases_passing_filter =
            passed_filter ? primary_aligned_bases + softclip_counts.aligned_soft_clipped_bases_passing_filter : 0,
        .soft_clipped_bases = soft_clipped_bases,
        .unmapped_soft_clipped_bases = soft_clipped_bases - softclip_counts.mapped_soft_clipped_bases,
        .gc_bases = io::CountGCBases(bam_get_seq(alignment.alignment_ptr), alignment.qpos, read_length)};

    ++read_metrics_summary.total_reads;
    RecordBaseCounts(base_counts, passed_filter);

    if (!mapped) {
      ++read_metrics_summary.unmapped_reads;
      read_metrics_summary.unmapped_bases += read_length;
    } else {
      // Mapping quality only makes sense for mapped reads
      if (alignment.alignment_ptr->core.qual == 0) {
        ++read_metrics_summary.mapq_zero_reads;
      }
      ++read_metrics_summary.mapped_reads;
      RecordReadLengthHistograms(
          read_length, read_length_without_softclips, alignment_metadata.read_type, passed_filter);
      RecordReadCountsByStrand(alignment_metadata.read_strand, passed_filter);
    }
    if (passed_filter) {
      ++read_metrics_summary.reads_passing_filter;
    }
    RecordReadCountsByType(alignment_metadata.read_type, passed_filter);
    RecordReadCountsByClusterType(alignment_metadata, passed_filter);
  }

  UpdateTargetEnrichmentReadMetrics();
}

void ReadMetrics::RecordAlignmentRecordCounts(const u32 flag) {
  // Total alignment records include all primary, secondary, supplementary alignments, and unmapped reads.
  // A read can have multiple alignment records associated with it.
  ++read_metrics_summary.alignments_plus_unmapped_reads;
  if ((flag & u32{BAM_FUNMAP}) == 0) {
    ++read_metrics_summary.alignment_records;
  }
  if ((flag & u32{BAM_FSUPPLEMENTARY}) != 0) {
    ++read_metrics_summary.supplementary_alignments;
  }
  if ((flag & u32{BAM_FSECONDARY}) != 0) {
    ++read_metrics_summary.secondary_alignments;
  }
  if ((flag & u32{BAM_FDUP}) != 0) {
    ++read_metrics_summary.duplicate_reads;
  }
}

void ReadMetrics::RecordBaseCounts(const ReadBaseCounts& base_counts, const bool passed_filter) {
  read_metrics_summary.total_bases += base_counts.read_length;
  read_metrics_summary.aligned_bases += base_counts.aligned_bases;
  read_metrics_summary.unique_aligned_bases += base_counts.unique_aligned_bases;
  read_metrics_summary.soft_clipped_bases += base_counts.soft_clipped_bases;
  read_metrics_summary.unmapped_soft_clipped_bases += base_counts.unmapped_soft_clipped_bases;
  read_metrics_summary.gc_bases += base_counts.gc_bases;
  if (passed_filter) {
    read_metrics_summary.total_bases_passing_filter += base_counts.read_length;
    read_metrics_summary.aligned_bases_passing_filter += base_counts.aligned_bases_passing_filter;
    read_metrics_summary.unique_aligned_bases_passing_filter += base_counts.unique_aligned_bases_passing_filter;
    read_metrics_summary.soft_clipped_bases_passing_filter += base_counts.soft_clipped_bases;
    read_metrics_summary.unmapped_soft_clipped_bases_passing_filter += base_counts.unmapped_soft_clipped_bases;
    read_metrics_summary.gc_bases_passing_filter += base_counts.gc_bases;
  }
}

void ReadMetrics::RecordReadLengthHistograms(const u32 read_length,
                                             const u32 read_length_without_softclips,
                                             const ReadType read_type,
                                             const bool passed_filter) {
  const auto update_read_length_histograms = [read_length, read_length_without_softclips, passed_filter](
                                                 histogram::Histogram<u64>& mapped_including_softclips,
                                                 histogram::Histogram<u64>& mapped_excluding_softclips,
                                                 histogram::Histogram<u64>& post_filter_including_softclips,
                                                 histogram::Histogram<u64>& post_filter_excluding_softclips) {
    mapped_including_softclips.AddCountToHistogram(read_length, u64{1});
    mapped_excluding_softclips.AddCountToHistogram(read_length_without_softclips, u64{1});
    if (passed_filter) {
      post_filter_including_softclips.AddCountToHistogram(read_length, u64{1});
      post_filter_excluding_softclips.AddCountToHistogram(read_length_without_softclips, u64{1});
    }
  };
  // Read length histograms count mapped reads (primary alignments) before read-level filters are applied.
  // Both variants are recorded: including soft-clipped bases (read_length) and excluding them
  // (read_length_without_softclips). The post_filter histograms count the subset that also passed filters.
  update_read_length_histograms(read_length_histograms.mapped_all_read_length_histogram,
                                read_length_histograms.mapped_all_read_length_excluding_softclips_histogram,
                                read_length_histograms.post_filter_all_read_length_histogram,
                                read_length_histograms.post_filter_all_read_length_excluding_softclips_histogram);
  // The full/partial histograms further distinguish mapped reads by read type.
  switch (read_type) {
    case ReadType::kFull:
      update_read_length_histograms(read_length_histograms.mapped_full_read_length_histogram,
                                    read_length_histograms.mapped_full_read_length_excluding_softclips_histogram,
                                    read_length_histograms.post_filter_full_read_length_histogram,
                                    read_length_histograms.post_filter_full_read_length_excluding_softclips_histogram);
      break;
    case ReadType::kPartial:
      update_read_length_histograms(
          read_length_histograms.mapped_partial_read_length_histogram,
          read_length_histograms.mapped_partial_read_length_excluding_softclips_histogram,
          read_length_histograms.post_filter_partial_read_length_histogram,
          read_length_histograms.post_filter_partial_read_length_excluding_softclips_histogram);
      break;
    case ReadType::kUnknown:
    default:
      break;
  }
}

void ReadMetrics::RecordReadCountsByStrand(const ReadStrand read_strand, const bool passed_filter) {
  // Increments the total counter and, when the read passed filters, the passing-filter counter.
  const auto record = [passed_filter](u64& total, u64& total_passing_filter) {
    ++total;
    if (passed_filter) {
      ++total_passing_filter;
    }
  };
  using enum ReadStrand;
  switch (read_strand) {
    case kForward:
      record(read_metrics_summary.forward_reads, read_metrics_summary.forward_reads_passing_filter);
      break;
    case kReverse:
      record(read_metrics_summary.reverse_reads, read_metrics_summary.reverse_reads_passing_filter);
      break;
    case kUnknown:
    default:
      break;
  }
}

void ReadMetrics::RecordReadCountsByType(const ReadType read_type, const bool passed_filter) {
  // Increments the total counter and, when the read passed filters, the passing-filter counter.
  const auto record = [passed_filter](std::optional<u64>& total, std::optional<u64>& total_passing_filter) {
    ++total;
    if (passed_filter) {
      ++total_passing_filter;
    }
  };
  using enum ReadType;
  switch (read_type) {
    case kFull:
      record(read_metrics_summary.full_length_reads, read_metrics_summary.full_length_reads_passing_filter);
      break;
    case kPartial:
      record(read_metrics_summary.partial_length_reads, read_metrics_summary.partial_length_reads_passing_filter);
      break;
    case kUnknown:
    default:
      break;
  }
}

void ReadMetrics::RecordReadCountsByClusterType(const AlignmentMetadata& alignment_metadata, const bool passed_filter) {
  if (!read_counts_by_cluster_size) {
    return;
  }
  auto& cluster_counts = read_counts_by_cluster_size.value();
  // For cluster sizes greater than the maximum cluster size, we use the last index in the vector.
  // (i.e. the max+ bin)
  const size_t cluster_size_index = alignment_metadata.cluster_size > cluster_counts.read_counts_by_cluster_size.size()
                                        ? cluster_counts.read_counts_by_cluster_size.size() - 1
                                        : math::SatSub(alignment_metadata.cluster_size, u32{1});
  auto& cluster_bin = cluster_counts[cluster_size_index];

  // Increments the total counter, and (when the read passed filters) the corresponding
  // passing-filter counters in both the summary and the per-cluster-size bin.
  const auto record =
      [passed_filter](std::optional<u64>& total, std::optional<u64>& total_passing_filter, u64& bin_passing_filter) {
        ++total;
        if (passed_filter) {
          ++total_passing_filter;
          ++bin_passing_filter;
        }
      };

  {
    using enum ClusterStrand;
    switch (alignment_metadata.cluster_strand) {
      case kForwardCluster:
        record(read_metrics_summary.forward_cluster_reads,
               read_metrics_summary.forward_cluster_reads_passing_filter,
               cluster_bin.forward_cluster_reads_passing_filter);
        break;
      case kReverseCluster:
        record(read_metrics_summary.reverse_cluster_reads,
               read_metrics_summary.reverse_cluster_reads_passing_filter,
               cluster_bin.reverse_cluster_reads_passing_filter);
        break;
      case kMixedStrandCluster:
        record(read_metrics_summary.mixed_strand_cluster_reads,
               read_metrics_summary.mixed_strand_cluster_reads_passing_filter,
               cluster_bin.mixed_strand_cluster_reads_passing_filter);
        break;
      case kUnknown:
      default:
        break;
    }
  }

  {
    using enum ClusterType;
    switch (alignment_metadata.cluster_type) {
      case kFull:
        record(read_metrics_summary.full_cluster_reads,
               read_metrics_summary.full_cluster_reads_passing_filter,
               cluster_bin.full_cluster_reads_passing_filter);
        break;
      case kPartial:
        record(read_metrics_summary.partial_cluster_reads,
               read_metrics_summary.partial_cluster_reads_passing_filter,
               cluster_bin.partial_cluster_reads_passing_filter);
        break;
      case kMixed:
        record(read_metrics_summary.mixed_full_and_partial_cluster_reads,
               read_metrics_summary.mixed_full_and_partial_cluster_reads_passing_filter,
               cluster_bin.mixed_full_and_partial_cluster_reads_passing_filter);
        break;
      case kUnknown:
      default:
        break;
    }
  }
}

void ReadMetrics::UpdateTargetEnrichmentReadMetrics() {
  if (!te_read_metrics.has_value()) {
    return;
  }
  // Use values from the merged final metrics to get on-target reads
  te_read_metrics->on_target_alignments = read_metrics_summary.mapped_reads +
                                          read_metrics_summary.secondary_alignments +
                                          read_metrics_summary.supplementary_alignments;
  te_read_metrics->on_target_reads_passing_filter = read_metrics_summary.reads_passing_filter;
}

void ReadMetrics::UpdateUnmappedReadMetrics(const io::AlignmentReader& alignment_reader,
                                            const bool passed_filter,
                                            const u16 trim) {
  if (_unmapped_metrics_updated) {
    return;
  }
  // Check if there is any unplaced unmapped reads before proceeding to avoid the
  // "EOF marker absent" warning
  if (hts_idx_get_n_no_coor(alignment_reader.idx.get()) == 0) {
    _unmapped_metrics_updated = true;
    return;
  }
  // Count bases in unmapped reads
  const RegionLookupTable empty_region_lookup_table{vec<SuperRegion>{}};
  auto itr = io::SamItrQueryI(alignment_reader.idx.get(), HTS_IDX_NOCOOR, 0, 0);
  const io::Bam1Ptr bam1_ptr(bam_init1());
  while (io::SamItrMultiNext(alignment_reader.bam_file.get(), itr.get(), bam1_ptr.get()) > 0) {
    Alignment alignment(bam1_ptr.get(), {.disable_hp_quality_modification = true, .disable_base_type_decoding = true});
    // Since unmapped reads have no CIGAR ops, it doesn't matter which side we trim from
    // as it would only affect the read length
    alignment.Trim(trim, 0);
    const AlignmentMetadata alignment_metadata = CreateAlignmentMetadata(alignment.alignment_ptr, false);
    // Unmapped reads carry no supplementary alignments, so the regions of interest are never
    // consulted; an empty dummy table is passed to satisfy the AddRead signature.
    AddRead(alignment, alignment_metadata, passed_filter, empty_region_lookup_table);
  }
  _unmapped_metrics_updated = true;
}

}  // namespace xoos::alignment_metrics
