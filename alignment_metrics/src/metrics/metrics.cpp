#include "metrics/metrics.h"

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <xoos/histogram/histogram-summary.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/io/htslib-util/read-util.h>
#include <xoos/types/int.h>
#include <xoos/util/math.h>

#include "alignment-metrics-options.h"
#include "core/alignment.h"
#include "io/alignment-reader.h"
#include "io/metrics-default-filenames.h"
#include "metadata/alignment-metadata.h"
#include "metadata/dataset-metadata.h"
#include "metrics/accuracy-metrics/qscore-stats.h"

namespace xoos::alignment_metrics {
AccuracyMetrics::AccuracyMetrics(const AlignmentMetricsOptions& options, const DatasetMetadata& dataset_metadata)
    : errors_by_cluster_size(dataset_metadata.has_cluster_info
                                 ? std::make_optional(ErrorsByClusterSize(options.max_cluster_size_bin))
                                 : std::nullopt),
      errors_by_substitution_type(ErrorsBySubstitutionType(dataset_metadata)),
      qscore_stats(QscoreStats(dataset_metadata)),
      errors_by_read_type((dataset_metadata.has_read_type_info || dataset_metadata.has_strand_info)
                              ? std::make_optional(ErrorsByReadType(dataset_metadata))
                              : std::nullopt),
      hp_stats(options.calculate_hp_metrics ? std::make_optional(HpStats(options.max_hp_length)) : std::nullopt) {
}

AccuracyMetrics& AccuracyMetrics::operator+=(const AccuracyMetrics& other) {
  base_level_accuracy_summary += other.base_level_accuracy_summary;
  if (errors_by_cluster_size && other.errors_by_cluster_size) {
    errors_by_cluster_size->operator+=(*other.errors_by_cluster_size);
  }
  errors_by_substitution_type += other.errors_by_substitution_type;
  qscore_stats += other.qscore_stats;
  if (errors_by_read_type && other.errors_by_read_type) {
    errors_by_read_type->operator+=(*other.errors_by_read_type);
  }
  if (hp_stats && other.hp_stats) {
    hp_stats->operator+=(*other.hp_stats);
  }
  return *this;
}

CoverageMetrics::CoverageMetrics(const AlignmentMetricsOptions& options, const DatasetMetadata& dataset_metadata)
    : coverage_histograms(CoverageHistograms(options.max_coverage_bin,
                                             options.coverage_cutoff,
                                             options.summary_stats_percentiles,
                                             dataset_metadata,
                                             !options.enable_te_metrics)),
      hp_coverage_histograms(options.calculate_hp_metrics
                                 ? std::make_optional(CoverageHistograms(options.max_coverage_bin,
                                                                         options.coverage_cutoff,
                                                                         options.summary_stats_percentiles,
                                                                         dataset_metadata,
                                                                         !options.enable_te_metrics))
                                 : std::nullopt),
      coverage_uniformity_metrics(options.enable_te_metrics ? std::make_optional(CoverageUniformityMetrics())
                                                            : std::nullopt) {
}

CoverageMetrics& CoverageMetrics::operator+=(const CoverageMetrics& other) {
  coverage_histograms.full_coverage_histogram += other.coverage_histograms.full_coverage_histogram;
  coverage_histograms.filtered_coverage_histogram += other.coverage_histograms.filtered_coverage_histogram;
  if (coverage_histograms.dataset_metadata.is_duplex_dataset &&
      other.coverage_histograms.dataset_metadata.is_duplex_dataset) {
    coverage_histograms.concordant_duplex_coverage_histogram +=
        other.coverage_histograms.concordant_duplex_coverage_histogram;
  }
  if (hp_coverage_histograms && other.hp_coverage_histograms) {
    hp_coverage_histograms->full_coverage_histogram += other.hp_coverage_histograms->full_coverage_histogram;
    hp_coverage_histograms->filtered_coverage_histogram += other.hp_coverage_histograms->filtered_coverage_histogram;
    if (hp_coverage_histograms->dataset_metadata.is_duplex_dataset &&
        other.hp_coverage_histograms->dataset_metadata.is_duplex_dataset) {
      hp_coverage_histograms->concordant_duplex_coverage_histogram +=
          other.hp_coverage_histograms->concordant_duplex_coverage_histogram;
    }
  }
  if (coverage_uniformity_metrics && other.coverage_uniformity_metrics) {
    coverage_uniformity_metrics->operator+=(*other.coverage_uniformity_metrics);
  }
  return *this;
}

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

Metrics::Metrics(const AlignmentMetricsOptions& options, const DatasetMetadata& dataset_metadata)
    : accuracy_metrics(options.metric_types.has_accuracy_metrics
                           ? std::make_optional(AccuracyMetrics(options, dataset_metadata))
                           : std::nullopt),
      coverage_metrics(options.metric_types.has_coverage_metrics
                           ? std::make_optional(CoverageMetrics(options, dataset_metadata))
                           : std::nullopt),
      read_metrics(options.metric_types.has_read_metrics ? std::make_optional(ReadMetrics(options, dataset_metadata))
                                                         : std::nullopt) {
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
                          const bool passed_filter) {
  const u32 read_length = alignment.read_length;
  const u32 read_length_without_softclips = alignment.read_length_without_softclips;
  const u32 soft_clipped_bases = read_length - read_length_without_softclips;
  const u32 gc_bases = io::CountGCBases(bam_get_seq(alignment.alignment_ptr), alignment.qpos, read_length);
  const u32 aligned_bases = alignment.CountAlignedBases();

  const bool mapped = (alignment.alignment_ptr->core.flag & BAM_FUNMAP) == 0;
  const bool primary = (alignment.alignment_ptr->core.flag & (BAM_FSUPPLEMENTARY | BAM_FSECONDARY)) == 0;

  // If a read has no length (e.g. after trimming), it is not counted
  if (read_length == 0) {
    return;
  }

  // Total alignment records include all primary, secondary, supplementary alignments, and unmapped reads.
  // A read can have multiple alignment records associated with it.
  ++read_metrics_summary.alignments_plus_unmapped_reads;
  if ((alignment.alignment_ptr->core.flag & BAM_FSUPPLEMENTARY) != 0) {
    ++read_metrics_summary.supplementary_alignments;
  }
  if ((alignment.alignment_ptr->core.flag & BAM_FSECONDARY) != 0) {
    ++read_metrics_summary.secondary_alignments;
  }
  if ((alignment.alignment_ptr->core.flag & BAM_FDUP) != 0) {
    ++read_metrics_summary.duplicate_reads;
  }
  // Count a read only once, ignoring secondary and supplementary alignments
  if (primary) {
    ++read_metrics_summary.total_reads;
    read_metrics_summary.total_bases += read_length;
    read_metrics_summary.aligned_bases += aligned_bases;
    read_metrics_summary.soft_clipped_bases += soft_clipped_bases;
    read_metrics_summary.gc_bases += gc_bases;
    if (!mapped) {
      ++read_metrics_summary.unmapped_reads;
      read_metrics_summary.unmapped_bases += read_length;
    } else {
      // Mapping quality only makes sense for mapped reads
      if (alignment.alignment_ptr->core.qual == 0) {
        ++read_metrics_summary.mapq_zero_reads;
      }
      ++read_metrics_summary.mapped_reads;
      // Read strand only makes sense for mapped reads
      switch (alignment_metadata.read_strand) {
        case ReadStrand::kForward:
          ++read_metrics_summary.forward_reads;
          if (passed_filter) {
            ++read_metrics_summary.forward_reads_passing_filter;
          }
          break;
        case ReadStrand::kReverse:
          ++read_metrics_summary.reverse_reads;
          if (passed_filter) {
            ++read_metrics_summary.reverse_reads_passing_filter;
          }
          break;
        case ReadStrand::kUnknown:
        default:
          break;
      }
    }
    if (passed_filter) {
      ++read_metrics_summary.reads_passing_filter;
      read_metrics_summary.total_bases_passing_filter += read_length;
      read_metrics_summary.aligned_bases_passing_filter += aligned_bases;
      read_metrics_summary.soft_clipped_bases_passing_filter += soft_clipped_bases;
      read_metrics_summary.gc_bases_passing_filter += gc_bases;
      if (mapped) {
        read_length_histograms.post_filter_all_read_length_histogram.AddCountToHistogram(read_length_without_softclips,
                                                                                         u64{1});
      }
    }
    switch (alignment_metadata.read_type) {
      case ReadType::kFull:
        ++read_metrics_summary.full_length_reads;
        if (passed_filter) {
          ++read_metrics_summary.full_length_reads_passing_filter;
          if (mapped) {
            read_length_histograms.post_filter_full_read_length_histogram.AddCountToHistogram(
                read_length_without_softclips, u64{1});
          }
        }
        break;
      case ReadType::kPartial:
        ++read_metrics_summary.partial_length_reads;
        if (passed_filter) {
          ++read_metrics_summary.partial_length_reads_passing_filter;
          if (mapped) {
            read_length_histograms.post_filter_partial_read_length_histogram.AddCountToHistogram(
                read_length_without_softclips, u64{1});
          }
        }
        break;
      case ReadType::kUnknown:
      default:
        break;
    }

    if (read_counts_by_cluster_size) {
      auto& read_counts_by_cluster_size = this->read_counts_by_cluster_size.value();
      // For cluster sizes greater than the maximum cluster size, we use the last index in the vector.
      // (i.e. the max+ bin)
      const size_t cluster_size_index =
          alignment_metadata.cluster_size > read_counts_by_cluster_size.read_counts_by_cluster_size.size()
              ? read_counts_by_cluster_size.read_counts_by_cluster_size.size() - 1
              : math::SatSub(alignment_metadata.cluster_size, u32{1});
      switch (alignment_metadata.cluster_strand) {
        case ClusterStrand::kForwardCluster:
          ++read_metrics_summary.forward_cluster_reads;
          if (passed_filter) {
            ++read_metrics_summary.forward_cluster_reads_passing_filter;
            ++read_counts_by_cluster_size[cluster_size_index].forward_cluster_reads_passing_filter;
          }
          break;
        case ClusterStrand::kReverseCluster:
          ++read_metrics_summary.reverse_cluster_reads;
          if (passed_filter) {
            ++read_metrics_summary.reverse_cluster_reads_passing_filter;
            ++read_counts_by_cluster_size[cluster_size_index].reverse_cluster_reads_passing_filter;
          }
          break;
        case ClusterStrand::kMixedStrandCluster:
          ++read_metrics_summary.mixed_strand_cluster_reads;
          if (passed_filter) {
            ++read_metrics_summary.mixed_strand_cluster_reads_passing_filter;
            ++read_counts_by_cluster_size[cluster_size_index].mixed_strand_cluster_reads_passing_filter;
          }
          break;
        case ClusterStrand::kUnknown:
        default:
          break;
      }
      switch (alignment_metadata.cluster_type) {
        case ClusterType::kFull:
          ++read_metrics_summary.full_cluster_reads;
          if (passed_filter) {
            ++read_metrics_summary.full_cluster_reads_passing_filter;
            ++read_counts_by_cluster_size[cluster_size_index].full_cluster_reads_passing_filter;
          }
          break;
        case ClusterType::kPartial:
          ++read_metrics_summary.partial_cluster_reads;
          if (passed_filter) {
            ++read_metrics_summary.partial_cluster_reads_passing_filter;
            ++read_counts_by_cluster_size[cluster_size_index].partial_cluster_reads_passing_filter;
          }
          break;
        case ClusterType::kMixed:
          ++read_metrics_summary.mixed_full_and_partial_cluster_reads;
          if (passed_filter) {
            ++read_metrics_summary.mixed_full_and_partial_cluster_reads_passing_filter;
            ++read_counts_by_cluster_size[cluster_size_index].mixed_full_and_partial_cluster_reads_passing_filter;
          }
          break;
        case ClusterType::kUnknown:
        default:
          break;
      }
    }
  }

  if (te_read_metrics.has_value()) {
    // Use values from the merged final metrics to get on-target reads
    te_read_metrics->on_target_alignments = read_metrics_summary.mapped_reads +
                                            read_metrics_summary.secondary_alignments +
                                            read_metrics_summary.supplementary_alignments;
    te_read_metrics->on_target_reads_passing_filter = read_metrics_summary.reads_passing_filter;
  }
}

void CoverageMetrics::AddEmptyHistogramData(const u32 num_bases) {
  coverage_histograms.full_coverage_histogram.AddCountToHistogram(0, u64{num_bases});
  coverage_histograms.filtered_coverage_histogram.AddCountToHistogram(0, u64{num_bases});
  if (coverage_histograms.dataset_metadata.is_duplex_dataset) {
    coverage_histograms.concordant_duplex_coverage_histogram.AddCountToHistogram(0, u64{num_bases});
  }
}

/**
 * @brief Adds coverage depth counts to coverage histograms for a single genomic position.
 *
 * @param all_depth Total depth including all bases (low-quality and discordant)
 * @param concordant_duplex_depth Depth of concordant duplex bases only
 * @param filtered_depth Depth of high-quality filtered bases
 * @param is_part_of_hp Whether this position is part of a homopolymer sequence
 */
void CoverageMetrics::AddHistogramData(const u32 all_depth,
                                       const u32 concordant_duplex_depth,
                                       const u32 filtered_depth,
                                       const bool is_part_of_hp) {
  if (coverage_histograms.dataset_metadata.is_duplex_dataset) {
    coverage_histograms.concordant_duplex_coverage_histogram.AddCountToHistogram(concordant_duplex_depth, u64{1});
    if (is_part_of_hp && hp_coverage_histograms) {
      hp_coverage_histograms->concordant_duplex_coverage_histogram.AddCountToHistogram(concordant_duplex_depth, u64{1});
    }
  }
  coverage_histograms.filtered_coverage_histogram.AddCountToHistogram(filtered_depth, u64{1});
  coverage_histograms.full_coverage_histogram.AddCountToHistogram(all_depth, u64{1});
  if (is_part_of_hp && hp_coverage_histograms) {
    hp_coverage_histograms->filtered_coverage_histogram.AddCountToHistogram(filtered_depth, u64{1});
    hp_coverage_histograms->full_coverage_histogram.AddCountToHistogram(all_depth, u64{1});
  }
}

/**
 * @brief Writes all metrics to TSV files in organized subdirectories.
 *
 * Creates output directory structure (accuracy/, coverage/, read_metrics/) and writes each
 * metric category to separate TSV files. Only writes metrics that were enabled and populated
 * during processing. Comments header is included in all output files for metadata tracking.
 *
 * @param output_dir Path to the output directory where metric files will be written
 * @param comments Header comments to include in all TSV files (e.g., command line, timestamps)
 */
void Metrics::WriteMetrics(const fs::path& output_dir, const io::Comments& comments) const {
  const fs::path absolute_out_dir = fs::absolute(output_dir);
  // Create output directories if it doesn't exist
  // These calls to `create_directories` will not throw even if the directories already exist and will simply return
  // false and do nothing.
  fs::create_directories(absolute_out_dir);
  if (accuracy_metrics.has_value()) {
    create_directories(absolute_out_dir / kDefaultAccuracyDir);
    accuracy_metrics->errors_by_substitution_type.WriteTsv(
        output_dir / kDefaultAccuracyDir / kDefaultErrorBySubstitutionType, comments);
    accuracy_metrics->base_level_accuracy_summary.WriteTsv(
        output_dir / kDefaultAccuracyDir / kDefaultBaseLevelAccuracySummary, comments);
    accuracy_metrics->qscore_stats.WriteTsv(output_dir / kDefaultAccuracyDir / kDefaultQscoreStats, comments);
    if (accuracy_metrics->errors_by_cluster_size.has_value()) {
      accuracy_metrics->errors_by_cluster_size->WriteTsv(output_dir / kDefaultAccuracyDir / kDefaultErrorByClusterSize,
                                                         comments);
    }
    if (accuracy_metrics->errors_by_read_type.has_value()) {
      accuracy_metrics->errors_by_read_type->WriteTsv(output_dir / kDefaultAccuracyDir / kDefaultErrorByReadType,
                                                      comments);
    }
    if (accuracy_metrics->hp_stats.has_value()) {
      accuracy_metrics->hp_stats->WriteTsv(output_dir / kDefaultAccuracyDir / kDefaultHpErrors, comments);
    }
  }
  if (coverage_metrics.has_value()) {
    create_directories(absolute_out_dir / kDefaultCoverageDir);
    coverage_metrics->coverage_histograms.WriteTsv(output_dir / kDefaultCoverageDir / kDefaultCoverageHistograms,
                                                   comments);
    coverage_metrics->coverage_histograms.WriteCoverageStatsTsv(
        output_dir / kDefaultCoverageDir / kDefaultCoverageStats, comments);
    coverage_metrics->coverage_histograms.WriteCoverageDistributionSummaryTsv(
        output_dir / kDefaultCoverageDir / kDefaultCoverageDistributionSummary, comments);
    if (coverage_metrics->hp_coverage_histograms.has_value()) {
      coverage_metrics->hp_coverage_histograms->WriteTsv(
          output_dir / kDefaultCoverageDir / kDefaultHpCoverageHistograms, comments);
      coverage_metrics->hp_coverage_histograms->WriteCoverageStatsTsv(
          output_dir / kDefaultCoverageDir / kDefaultHpCoverageStats, comments);
      coverage_metrics->hp_coverage_histograms->WriteCoverageDistributionSummaryTsv(
          output_dir / kDefaultCoverageDir / kDefaultHpCoverageDistributionSummary, comments);
    }
    if (coverage_metrics->coverage_uniformity_metrics.has_value()) {
      coverage_metrics->coverage_uniformity_metrics->WriteMeanCoverageHistogram(
          output_dir / kDefaultCoverageDir / kDefaultMeanCoverageHistogram, comments);
      coverage_metrics->coverage_uniformity_metrics->WriteCoverageUniformitySummaryTsv(
          output_dir / kDefaultCoverageDir / kDefaultCoverageUniformitySummary,
          comments,
          coverage_metrics->coverage_histograms.filtered_coverage_histogram);
    }
  }
  if (read_metrics.has_value()) {
    create_directories(absolute_out_dir / kDefaultReadMetricsDir);
    if (read_metrics->read_counts_by_cluster_size.has_value()) {
      read_metrics->read_counts_by_cluster_size->WriteTsv(
          output_dir / kDefaultReadMetricsDir / kDefaultReadCountsByClusterSize, comments);
    }
    if (read_metrics->te_read_metrics.has_value()) {
      read_metrics->te_read_metrics->WriteTsv(output_dir / kDefaultReadMetricsDir / kDefaultTeReadMetrics, comments);
    }
    read_metrics->read_metrics_summary.WriteTsv(output_dir / kDefaultReadMetricsDir / kDefaultReadMetricsSummary,
                                                comments);
    read_metrics->read_length_histograms.WriteTsv(output_dir / kDefaultReadMetricsDir / kDefaultReadLengthHistograms,
                                                  comments);
    read_metrics->read_length_histograms.WriteSummaryTsv(
        output_dir / kDefaultReadMetricsDir / kDefaultReadLengthSummary, comments);
  }
}

void ReadMetrics::UpdateUnmappedReadMetrics(const AlignmentReader& alignment_reader,
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
  auto itr = io::SamItrQueryI(alignment_reader.idx.get(), HTS_IDX_NOCOOR, 0, 0);
  const io::Bam1Ptr bam1_ptr(bam_init1());
  while (io::SamItrMultiNext(alignment_reader.bam.get(), itr.get(), bam1_ptr.get()) > 0) {
    Alignment alignment(bam1_ptr.get(), {.disable_hp_quality_modification = true, .disable_base_type_decoding = true});
    // Since unmapped reads have no CIGAR ops, it doesn't matter which side we trim from
    // as it would only affect the read length
    alignment.Trim(trim, 0);
    AlignmentMetadata alignment_metadata = CreateAlignmentMetadata(alignment.alignment_ptr, false);
    AddRead(alignment, alignment_metadata, passed_filter);
  }
  _unmapped_metrics_updated = true;
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
  return *this;
}

}  // namespace xoos::alignment_metrics
