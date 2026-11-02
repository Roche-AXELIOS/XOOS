#include "metrics/read-metrics/read-metrics-summary.h"

#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>

#include <csv.hpp>

#include <xoos/log/logging.h>
#include <xoos/types/int.h>
#include <xoos/util/format-functions.h>
#include <xoos/util/missing-value.h>

#include "metrics/metrics-names.h"

namespace xoos::alignment_metrics {

using util::kNA;
using util::ToPercentageWithPrecision;

namespace {

// Number of decimal places used for all percentages in the read metrics summary.
constexpr u8 kPercentagePrecision = 2;

// Writes the pre-filter read count rows (mapping, strand, cluster stratification).
template <typename Writer>
void WriteReadCounts(Writer& writer, const ReadMetricsSummary& summary) {
  writer << std::make_tuple(kNameTotalAlignmentsIncludingUnmapped, summary.alignments_plus_unmapped_reads, kNA, kNA);
  writer << std::make_tuple(
      kNameAlignmentRecords,
      summary.alignment_records,
      kNameTotalAlignmentsIncludingUnmapped,
      ToPercentageWithPrecision(
          summary.alignment_records, summary.alignments_plus_unmapped_reads, kPercentagePrecision));
  writer << std::make_tuple(
      kNameSecondaryAlignments,
      summary.secondary_alignments,
      kNameTotalAlignmentsIncludingUnmapped,
      ToPercentageWithPrecision(
          summary.secondary_alignments, summary.alignments_plus_unmapped_reads, kPercentagePrecision));
  writer << std::make_tuple(
      kNameSupplementaryAlignments,
      summary.supplementary_alignments,
      kNameTotalAlignmentsIncludingUnmapped,
      ToPercentageWithPrecision(
          summary.supplementary_alignments, summary.alignments_plus_unmapped_reads, kPercentagePrecision));
  // Pre-filter read counts
  writer << std::make_tuple(
      kNameTotalReads,
      summary.total_reads,
      kNameTotalAlignmentsIncludingUnmapped,
      ToPercentageWithPrecision(summary.total_reads, summary.alignments_plus_unmapped_reads, kPercentagePrecision));
  writer << std::make_tuple(
      kNameUnmappedReads,
      summary.unmapped_reads,
      kNameTotalReads,
      ToPercentageWithPrecision(summary.unmapped_reads, summary.total_reads, kPercentagePrecision));
  writer << std::make_tuple(kNameMappedReads,
                            summary.mapped_reads,
                            kNameTotalReads,
                            ToPercentageWithPrecision(summary.mapped_reads, summary.total_reads, kPercentagePrecision));
  writer << std::make_tuple(
      kNameForwardReads,
      summary.forward_reads,
      kNameMappedReads,
      ToPercentageWithPrecision(summary.forward_reads, summary.mapped_reads, kPercentagePrecision));
  writer << std::make_tuple(
      kNameReverseReads,
      summary.reverse_reads,
      kNameMappedReads,
      ToPercentageWithPrecision(summary.reverse_reads, summary.mapped_reads, kPercentagePrecision));
  writer << std::make_tuple(
      kNameDuplicateReads,
      summary.duplicate_reads,
      kNameMappedReads,
      ToPercentageWithPrecision(summary.duplicate_reads, summary.mapped_reads, kPercentagePrecision));
  writer << std::make_tuple(
      kNameMapQZeroReads,
      summary.mapq_zero_reads,
      kNameMappedReads,
      ToPercentageWithPrecision(summary.mapq_zero_reads, summary.mapped_reads, kPercentagePrecision));
  if (summary.full_length_reads.has_value()) {
    writer << std::make_tuple(
        kNameFullLengthReads,
        *summary.full_length_reads,
        kNameTotalReads,
        ToPercentageWithPrecision(*summary.full_length_reads, summary.total_reads, kPercentagePrecision));
  }
  if (summary.partial_length_reads.has_value()) {
    writer << std::make_tuple(
        kNamePartialLengthReads,
        *summary.partial_length_reads,
        kNameTotalReads,
        ToPercentageWithPrecision(*summary.partial_length_reads, summary.total_reads, kPercentagePrecision));
  }
  if (summary.mixed_strand_cluster_reads.has_value()) {
    writer << std::make_tuple(
        kNameMixedStrandClusterReads,
        *summary.mixed_strand_cluster_reads,
        kNameTotalReads,
        ToPercentageWithPrecision(*summary.mixed_strand_cluster_reads, summary.total_reads, kPercentagePrecision));
  }
  if (summary.forward_cluster_reads.has_value()) {
    writer << std::make_tuple(
        kNameForwardClusterReads,
        *summary.forward_cluster_reads,
        kNameTotalReads,
        ToPercentageWithPrecision(*summary.forward_cluster_reads, summary.total_reads, kPercentagePrecision));
  }
  if (summary.reverse_cluster_reads.has_value()) {
    writer << std::make_tuple(
        kNameReverseClusterReads,
        *summary.reverse_cluster_reads,
        kNameTotalReads,
        ToPercentageWithPrecision(*summary.reverse_cluster_reads, summary.total_reads, kPercentagePrecision));
  }
  if (summary.full_cluster_reads.has_value()) {
    writer << std::make_tuple(
        kNameFullClusterReads,
        *summary.full_cluster_reads,
        kNameTotalReads,
        ToPercentageWithPrecision(*summary.full_cluster_reads, summary.total_reads, kPercentagePrecision));
  }
  if (summary.mixed_full_and_partial_cluster_reads.has_value()) {
    writer << std::make_tuple(
        kNameMixedFullAndPartialClusterReads,
        *summary.mixed_full_and_partial_cluster_reads,
        kNameTotalReads,
        ToPercentageWithPrecision(
            *summary.mixed_full_and_partial_cluster_reads, summary.total_reads, kPercentagePrecision));
  }
  if (summary.partial_cluster_reads.has_value()) {
    writer << std::make_tuple(
        kNamePartialClusterReads,
        *summary.partial_cluster_reads,
        kNameTotalReads,
        ToPercentageWithPrecision(*summary.partial_cluster_reads, summary.total_reads, kPercentagePrecision));
  }
}

// Writes the post-filter read count rows.
template <typename Writer>
void WriteReadCountsPassingFilter(Writer& writer, const ReadMetricsSummary& summary) {
  writer << std::make_tuple(
      kNameReadsPassingFilter,
      summary.reads_passing_filter,
      kNameTotalReads,
      ToPercentageWithPrecision(summary.reads_passing_filter, summary.total_reads, kPercentagePrecision));
  writer << std::make_tuple(
      kNameForwardReadsPassingFilter,
      summary.forward_reads_passing_filter,
      kNameReadsPassingFilter,
      ToPercentageWithPrecision(
          summary.forward_reads_passing_filter, summary.reads_passing_filter, kPercentagePrecision));
  writer << std::make_tuple(
      kNameReverseReadsPassingFilter,
      summary.reverse_reads_passing_filter,
      kNameReadsPassingFilter,
      ToPercentageWithPrecision(
          summary.reverse_reads_passing_filter, summary.reads_passing_filter, kPercentagePrecision));
  if (summary.full_length_reads_passing_filter.has_value()) {
    writer << std::make_tuple(
        kNameFullLengthReadsPassingFilter,
        *summary.full_length_reads_passing_filter,
        kNameReadsPassingFilter,
        ToPercentageWithPrecision(
            *summary.full_length_reads_passing_filter, summary.reads_passing_filter, kPercentagePrecision));
  }
  if (summary.partial_length_reads_passing_filter.has_value()) {
    writer << std::make_tuple(
        kNamePartialLengthReadsPassingFilter,
        *summary.partial_length_reads_passing_filter,
        kNameReadsPassingFilter,
        ToPercentageWithPrecision(
            *summary.partial_length_reads_passing_filter, summary.reads_passing_filter, kPercentagePrecision));
  }
  if (summary.mixed_strand_cluster_reads_passing_filter.has_value()) {
    writer << std::make_tuple(
        kNameMixedStrandClusterReadsPassingFilter,
        *summary.mixed_strand_cluster_reads_passing_filter,
        kNameReadsPassingFilter,
        ToPercentageWithPrecision(
            *summary.mixed_strand_cluster_reads_passing_filter, summary.reads_passing_filter, kPercentagePrecision));
  }
  if (summary.forward_cluster_reads_passing_filter.has_value()) {
    writer << std::make_tuple(
        kNameForwardClusterReadsPassingFilter,
        *summary.forward_cluster_reads_passing_filter,
        kNameReadsPassingFilter,
        ToPercentageWithPrecision(
            *summary.forward_cluster_reads_passing_filter, summary.reads_passing_filter, kPercentagePrecision));
  }
  if (summary.reverse_cluster_reads_passing_filter.has_value()) {
    writer << std::make_tuple(
        kNameReverseClusterReadsPassingFilter,
        *summary.reverse_cluster_reads_passing_filter,
        kNameReadsPassingFilter,
        ToPercentageWithPrecision(
            *summary.reverse_cluster_reads_passing_filter, summary.reads_passing_filter, kPercentagePrecision));
  }
  if (summary.full_cluster_reads_passing_filter.has_value()) {
    writer << std::make_tuple(
        kNameFullClusterReadsPassingFilter,
        *summary.full_cluster_reads_passing_filter,
        kNameReadsPassingFilter,
        ToPercentageWithPrecision(
            *summary.full_cluster_reads_passing_filter, summary.reads_passing_filter, kPercentagePrecision));
  }
  if (summary.mixed_full_and_partial_cluster_reads_passing_filter.has_value()) {
    writer << std::make_tuple(kNameMixedFullAndPartialClusterReadsPassingFilter,
                              *summary.mixed_full_and_partial_cluster_reads_passing_filter,
                              kNameReadsPassingFilter,
                              ToPercentageWithPrecision(*summary.mixed_full_and_partial_cluster_reads_passing_filter,
                                                        summary.reads_passing_filter,
                                                        kPercentagePrecision));
  }
  if (summary.partial_cluster_reads_passing_filter.has_value()) {
    writer << std::make_tuple(
        kNamePartialClusterReadsPassingFilter,
        *summary.partial_cluster_reads_passing_filter,
        kNameReadsPassingFilter,
        ToPercentageWithPrecision(
            *summary.partial_cluster_reads_passing_filter, summary.reads_passing_filter, kPercentagePrecision));
  }
}

// Writes the base count rows (pre- and post-filter).
template <typename Writer>
void WriteBaseCounts(Writer& writer, const ReadMetricsSummary& summary) {
  writer << std::make_tuple(kNameTotalBases, summary.total_bases, kNA, kNA);
  writer << std::make_tuple(
      kNameUnmappedBases,
      summary.unmapped_bases,
      kNameTotalBases,
      ToPercentageWithPrecision(summary.unmapped_bases, summary.total_bases, kPercentagePrecision));
  // aligned_bases includes bases from supplementary alignments and may double-count a read base
  // aligned in both the primary and a supplementary alignment, so a percentage of total_bases is
  // not meaningful (it can exceed 100%). Report unique_aligned_bases as the percentage instead.
  writer << std::make_tuple(kNameAlignedBases, summary.aligned_bases, kNameTotalBases, kNA);
  writer << std::make_tuple(
      kNameUniqueAlignedBases,
      summary.unique_aligned_bases,
      kNameTotalBases,
      ToPercentageWithPrecision(summary.unique_aligned_bases, summary.total_bases, kPercentagePrecision));
  writer << std::make_tuple(
      kNameSoftClippedBases,
      summary.soft_clipped_bases,
      kNameTotalBases,
      ToPercentageWithPrecision(summary.soft_clipped_bases, summary.total_bases, kPercentagePrecision));
  writer << std::make_tuple(
      kNameUnmappedSoftClippedBases,
      summary.unmapped_soft_clipped_bases,
      kNameTotalBases,
      ToPercentageWithPrecision(summary.unmapped_soft_clipped_bases, summary.total_bases, kPercentagePrecision));
  writer << std::make_tuple(kNameGCBases,
                            summary.gc_bases,
                            kNameTotalBases,
                            ToPercentageWithPrecision(summary.gc_bases, summary.total_bases, kPercentagePrecision));
  writer << std::make_tuple(
      kNameTotalBasesPassingFilter,
      summary.total_bases_passing_filter,
      kNameTotalBases,
      ToPercentageWithPrecision(summary.total_bases_passing_filter, summary.total_bases, kPercentagePrecision));
  // As with aligned_bases, aligned_bases_passing_filter may double-count, so its percentage is NA
  // and unique_aligned_bases_passing_filter carries the percentage of total_bases_passing_filter.
  writer << std::make_tuple(
      kNameAlignedBasesPassingFilter, summary.aligned_bases_passing_filter, kNameTotalBasesPassingFilter, kNA);
  writer << std::make_tuple(
      kNameUniqueAlignedBasesPassingFilter,
      summary.unique_aligned_bases_passing_filter,
      kNameTotalBasesPassingFilter,
      ToPercentageWithPrecision(
          summary.unique_aligned_bases_passing_filter, summary.total_bases_passing_filter, kPercentagePrecision));
  writer << std::make_tuple(
      kNameSoftClippedBasesPassingFilter,
      summary.soft_clipped_bases_passing_filter,
      kNameTotalBasesPassingFilter,
      ToPercentageWithPrecision(
          summary.soft_clipped_bases_passing_filter, summary.total_bases_passing_filter, kPercentagePrecision));
  writer << std::make_tuple(kNameUnmappedSoftClippedBasesPassingFilter,
                            summary.unmapped_soft_clipped_bases_passing_filter,
                            kNameTotalBasesPassingFilter,
                            ToPercentageWithPrecision(summary.unmapped_soft_clipped_bases_passing_filter,
                                                      summary.total_bases_passing_filter,
                                                      kPercentagePrecision));
  writer << std::make_tuple(
      kNameGCBasesPassingFilter,
      summary.gc_bases_passing_filter,
      kNameTotalBasesPassingFilter,
      ToPercentageWithPrecision(
          summary.gc_bases_passing_filter, summary.total_bases_passing_filter, kPercentagePrecision));
}

}  // namespace

ReadMetricsSummary::ReadMetricsSummary(const DatasetMetadata& dataset_metadata) {
  if (dataset_metadata.has_cluster_info) {
    mixed_strand_cluster_reads = 0;
    forward_cluster_reads = 0;
    reverse_cluster_reads = 0;
    full_cluster_reads = 0;
    mixed_full_and_partial_cluster_reads = 0;
    partial_cluster_reads = 0;
    mixed_strand_cluster_reads_passing_filter = 0;
    forward_cluster_reads_passing_filter = 0;
    reverse_cluster_reads_passing_filter = 0;
    full_cluster_reads_passing_filter = 0;
    mixed_full_and_partial_cluster_reads_passing_filter = 0;
    partial_cluster_reads_passing_filter = 0;
  }
  if (dataset_metadata.has_read_type_info) {
    full_length_reads = 0;
    partial_length_reads = 0;
    full_length_reads_passing_filter = 0;
    partial_length_reads_passing_filter = 0;
  }
}

ReadMetricsSummary& ReadMetricsSummary::operator+=(const ReadMetricsSummary& obj) {
  this->alignments_plus_unmapped_reads += obj.alignments_plus_unmapped_reads;
  this->alignment_records += obj.alignment_records;
  this->total_reads += obj.total_reads;
  this->unmapped_reads += obj.unmapped_reads;
  this->mapped_reads += obj.mapped_reads;
  this->forward_reads += obj.forward_reads;
  this->reverse_reads += obj.reverse_reads;
  this->full_length_reads += obj.full_length_reads;
  this->partial_length_reads += obj.partial_length_reads;
  this->mixed_strand_cluster_reads += obj.mixed_strand_cluster_reads;
  this->forward_cluster_reads += obj.forward_cluster_reads;
  this->reverse_cluster_reads += obj.reverse_cluster_reads;
  this->full_cluster_reads += obj.full_cluster_reads;
  this->mixed_full_and_partial_cluster_reads += obj.mixed_full_and_partial_cluster_reads;
  this->partial_cluster_reads += obj.partial_cluster_reads;
  this->secondary_alignments += obj.secondary_alignments;
  this->supplementary_alignments += obj.supplementary_alignments;
  this->duplicate_reads += obj.duplicate_reads;
  this->mapq_zero_reads += obj.mapq_zero_reads;
  this->reads_passing_filter += obj.reads_passing_filter;
  this->forward_reads_passing_filter += obj.forward_reads_passing_filter;
  this->reverse_reads_passing_filter += obj.reverse_reads_passing_filter;
  this->full_length_reads_passing_filter += obj.full_length_reads_passing_filter;
  this->partial_length_reads_passing_filter += obj.partial_length_reads_passing_filter;
  this->mixed_strand_cluster_reads_passing_filter += obj.mixed_strand_cluster_reads_passing_filter;
  this->forward_cluster_reads_passing_filter += obj.forward_cluster_reads_passing_filter;
  this->reverse_cluster_reads_passing_filter += obj.reverse_cluster_reads_passing_filter;
  this->full_cluster_reads_passing_filter += obj.full_cluster_reads_passing_filter;
  this->mixed_full_and_partial_cluster_reads_passing_filter += obj.mixed_full_and_partial_cluster_reads_passing_filter;
  this->partial_cluster_reads_passing_filter += obj.partial_cluster_reads_passing_filter;
  this->total_bases += obj.total_bases;
  this->unmapped_bases += obj.unmapped_bases;
  this->aligned_bases += obj.aligned_bases;
  this->unique_aligned_bases += obj.unique_aligned_bases;
  this->soft_clipped_bases += obj.soft_clipped_bases;
  this->unmapped_soft_clipped_bases += obj.unmapped_soft_clipped_bases;
  this->gc_bases += obj.gc_bases;
  this->total_bases_passing_filter += obj.total_bases_passing_filter;
  this->aligned_bases_passing_filter += obj.aligned_bases_passing_filter;
  this->unique_aligned_bases_passing_filter += obj.unique_aligned_bases_passing_filter;
  this->soft_clipped_bases_passing_filter += obj.soft_clipped_bases_passing_filter;
  this->unmapped_soft_clipped_bases_passing_filter += obj.unmapped_soft_clipped_bases_passing_filter;
  this->gc_bases_passing_filter += obj.gc_bases_passing_filter;
  return *this;
}

std::vector<std::string> ReadMetricsSummary::GetHeaders() {
  return {kNameMetricName, kNameCount, kNameDenominator, kNamePercentage};
}

void ReadMetricsSummary::WriteTsv(const fs::path& output_path, const io::TsvMetadata& tsv_metadata) const {
  std::ofstream ofstream(output_path);
  auto writer = csv::make_tsv_writer(ofstream).set_auto_flush(false);
  io::WriteTsvMetadata(ofstream, tsv_metadata);
  writer << GetHeaders();
  WriteReadCounts(writer, *this);
  WriteReadCountsPassingFilter(writer, *this);
  WriteBaseCounts(writer, *this);
}

}  // namespace xoos::alignment_metrics
