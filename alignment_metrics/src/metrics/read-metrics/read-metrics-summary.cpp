#include "metrics/read-metrics/read-metrics-summary.h"

#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>

#include <csv.hpp>

#include <xoos/log/logging.h>
#include <xoos/types/int.h>

#include "metrics/metrics-names.h"
#include "util/format-util.h"

namespace xoos::alignment_metrics {

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
  this->soft_clipped_bases += obj.soft_clipped_bases;
  this->gc_bases += obj.gc_bases;
  this->total_bases_passing_filter += obj.total_bases_passing_filter;
  this->aligned_bases_passing_filter += obj.aligned_bases_passing_filter;
  this->soft_clipped_bases_passing_filter += obj.soft_clipped_bases_passing_filter;
  this->gc_bases_passing_filter += obj.gc_bases_passing_filter;
  return *this;
}

std::vector<std::string> ReadMetricsSummary::GetHeaders() {
  return {kNameMetricName, kNameCount, kNameDenominator, kNamePercentage};
}

void ReadMetricsSummary::WriteTsv(const fs::path& output_path, const io::Comments& comments) const {
  std::ofstream ofstream(output_path);
  auto writer = csv::make_tsv_writer_buffered(ofstream);
  io::WriteTsvComments(writer, comments);
  writer << GetHeaders();
  writer << std::make_tuple(
      kNameTotalAlignmentsIncludingUnmapped, alignments_plus_unmapped_reads, kNotApplicable, kNotApplicable);
  writer << std::make_tuple(kNameSecondaryAlignments,
                            secondary_alignments,
                            kNameTotalAlignmentsIncludingUnmapped,
                            ToPercentageWithPrecision(secondary_alignments, alignments_plus_unmapped_reads, 2));
  writer << std::make_tuple(kNameSupplementaryAlignments,
                            supplementary_alignments,
                            kNameTotalAlignmentsIncludingUnmapped,
                            ToPercentageWithPrecision(supplementary_alignments, alignments_plus_unmapped_reads, 2));
  // Pre-filter read counts
  writer << std::make_tuple(kNameTotalReads,
                            total_reads,
                            kNameTotalAlignmentsIncludingUnmapped,
                            ToPercentageWithPrecision(total_reads, alignments_plus_unmapped_reads, 2));
  writer << std::make_tuple(
      kNameUnmappedReads, unmapped_reads, kNameTotalReads, ToPercentageWithPrecision(unmapped_reads, total_reads, 2));
  writer << std::make_tuple(
      kNameMappedReads, mapped_reads, kNameTotalReads, ToPercentageWithPrecision(mapped_reads, total_reads, 2));
  writer << std::make_tuple(
      kNameForwardReads, forward_reads, kNameMappedReads, ToPercentageWithPrecision(forward_reads, mapped_reads, 2));
  writer << std::make_tuple(
      kNameReverseReads, reverse_reads, kNameMappedReads, ToPercentageWithPrecision(reverse_reads, mapped_reads, 2));
  writer << std::make_tuple(kNameDuplicateReads,
                            duplicate_reads,
                            kNameMappedReads,
                            ToPercentageWithPrecision(duplicate_reads, mapped_reads, 2));
  writer << std::make_tuple(kNameMapQZeroReads,
                            mapq_zero_reads,
                            kNameMappedReads,
                            ToPercentageWithPrecision(mapq_zero_reads, mapped_reads, 2));
  if (full_length_reads.has_value()) {
    writer << std::make_tuple(kNameFullLengthReads,
                              *full_length_reads,
                              kNameTotalReads,
                              ToPercentageWithPrecision(*full_length_reads, total_reads, 2));
  }
  if (partial_length_reads.has_value()) {
    writer << std::make_tuple(kNamePartialLengthReads,
                              *partial_length_reads,
                              kNameTotalReads,
                              ToPercentageWithPrecision(*partial_length_reads, total_reads, 2));
  }
  if (mixed_strand_cluster_reads.has_value()) {
    writer << std::make_tuple(kNameMixedStrandClusterReads,
                              *mixed_strand_cluster_reads,
                              kNameTotalReads,
                              ToPercentageWithPrecision(*mixed_strand_cluster_reads, total_reads, 2));
  }
  if (forward_cluster_reads.has_value()) {
    writer << std::make_tuple(kNameForwardClusterReads,
                              *forward_cluster_reads,
                              kNameTotalReads,
                              ToPercentageWithPrecision(*forward_cluster_reads, total_reads, 2));
  }
  if (reverse_cluster_reads.has_value()) {
    writer << std::make_tuple(kNameReverseClusterReads,
                              *reverse_cluster_reads,
                              kNameTotalReads,
                              ToPercentageWithPrecision(*reverse_cluster_reads, total_reads, 2));
  }
  if (full_cluster_reads.has_value()) {
    writer << std::make_tuple(kNameFullClusterReads,
                              *full_cluster_reads,
                              kNameTotalReads,
                              ToPercentageWithPrecision(*full_cluster_reads, total_reads, 2));
  }
  if (mixed_full_and_partial_cluster_reads.has_value()) {
    writer << std::make_tuple(kNameMixedFullAndPartialClusterReads,
                              *mixed_full_and_partial_cluster_reads,
                              kNameTotalReads,
                              ToPercentageWithPrecision(*mixed_full_and_partial_cluster_reads, total_reads, 2));
  }
  if (partial_cluster_reads.has_value()) {
    writer << std::make_tuple(kNamePartialClusterReads,
                              *partial_cluster_reads,
                              kNameTotalReads,
                              ToPercentageWithPrecision(*partial_cluster_reads, total_reads, 2));
  }

  // Post-filter read counts
  writer << std::make_tuple(kNameReadsPassingFilter,
                            reads_passing_filter,
                            kNameTotalReads,
                            ToPercentageWithPrecision(reads_passing_filter, total_reads, 2));
  writer << std::make_tuple(kNameForwardReadsPassingFilter,
                            forward_reads_passing_filter,
                            kNameReadsPassingFilter,
                            ToPercentageWithPrecision(forward_reads_passing_filter, reads_passing_filter, 2));
  writer << std::make_tuple(kNameReverseReadsPassingFilter,
                            reverse_reads_passing_filter,
                            kNameReadsPassingFilter,
                            ToPercentageWithPrecision(reverse_reads_passing_filter, reads_passing_filter, 2));
  if (full_length_reads_passing_filter.has_value()) {
    writer << std::make_tuple(kNameFullLengthReadsPassingFilter,
                              *full_length_reads_passing_filter,
                              kNameReadsPassingFilter,
                              ToPercentageWithPrecision(*full_length_reads_passing_filter, reads_passing_filter, 2));
  }
  if (partial_length_reads_passing_filter.has_value()) {
    writer << std::make_tuple(kNamePartialLengthReadsPassingFilter,
                              *partial_length_reads_passing_filter,
                              kNameReadsPassingFilter,
                              ToPercentageWithPrecision(*partial_length_reads_passing_filter, reads_passing_filter, 2));
  }
  if (mixed_strand_cluster_reads_passing_filter.has_value()) {
    writer << std::make_tuple(
        kNameMixedStrandClusterReadsPassingFilter,
        *mixed_strand_cluster_reads_passing_filter,
        kNameReadsPassingFilter,
        ToPercentageWithPrecision(*mixed_strand_cluster_reads_passing_filter, reads_passing_filter, 2));
  }
  if (forward_cluster_reads_passing_filter.has_value()) {
    writer << std::make_tuple(
        kNameForwardClusterReadsPassingFilter,
        *forward_cluster_reads_passing_filter,
        kNameReadsPassingFilter,
        ToPercentageWithPrecision(*forward_cluster_reads_passing_filter, reads_passing_filter, 2));
  }
  if (reverse_cluster_reads_passing_filter.has_value()) {
    writer << std::make_tuple(
        kNameReverseClusterReadsPassingFilter,
        *reverse_cluster_reads_passing_filter,
        kNameReadsPassingFilter,
        ToPercentageWithPrecision(*reverse_cluster_reads_passing_filter, reads_passing_filter, 2));
  }
  if (full_cluster_reads_passing_filter.has_value()) {
    writer << std::make_tuple(kNameFullClusterReadsPassingFilter,
                              *full_cluster_reads_passing_filter,
                              kNameReadsPassingFilter,
                              ToPercentageWithPrecision(*full_cluster_reads_passing_filter, reads_passing_filter, 2));
  }
  if (mixed_full_and_partial_cluster_reads_passing_filter.has_value()) {
    writer << std::make_tuple(
        kNameMixedFullAndPartialClusterReadsPassingFilter,
        *mixed_full_and_partial_cluster_reads_passing_filter,
        kNameReadsPassingFilter,
        ToPercentageWithPrecision(*mixed_full_and_partial_cluster_reads_passing_filter, reads_passing_filter, 2));
  }
  if (partial_cluster_reads_passing_filter.has_value()) {
    writer << std::make_tuple(
        kNamePartialClusterReadsPassingFilter,
        *partial_cluster_reads_passing_filter,
        kNameReadsPassingFilter,
        ToPercentageWithPrecision(*partial_cluster_reads_passing_filter, reads_passing_filter, 2));
  }
  writer << std::make_tuple(kNameTotalBases, total_bases, kNotApplicable, kNotApplicable);
  writer << std::make_tuple(
      kNameUnmappedBases, unmapped_bases, kNameTotalBases, ToPercentageWithPrecision(unmapped_bases, total_bases, 2));
  writer << std::make_tuple(
      kNameAlignedBases, aligned_bases, kNameTotalBases, ToPercentageWithPrecision(aligned_bases, total_bases, 2));
  writer << std::make_tuple(kNameSoftClippedBases,
                            soft_clipped_bases,
                            kNameTotalBases,
                            ToPercentageWithPrecision(soft_clipped_bases, total_bases, 2));
  writer << std::make_tuple(
      kNameGCBases, gc_bases, kNameTotalBases, ToPercentageWithPrecision(gc_bases, total_bases, 2));
  writer << std::make_tuple(kNameTotalBasesPassingFilter,
                            total_bases_passing_filter,
                            kNameTotalBases,
                            ToPercentageWithPrecision(total_bases_passing_filter, total_bases, 2));
  writer << std::make_tuple(kNameAlignedBasesPassingFilter,
                            aligned_bases_passing_filter,
                            kNameTotalBasesPassingFilter,
                            ToPercentageWithPrecision(aligned_bases_passing_filter, total_bases_passing_filter, 2));
  writer << std::make_tuple(
      kNameSoftClippedBasesPassingFilter,
      soft_clipped_bases_passing_filter,
      kNameTotalBasesPassingFilter,
      ToPercentageWithPrecision(soft_clipped_bases_passing_filter, total_bases_passing_filter, 2));
  writer << std::make_tuple(kNameGCBasesPassingFilter,
                            gc_bases_passing_filter,
                            kNameTotalBasesPassingFilter,
                            ToPercentageWithPrecision(gc_bases_passing_filter, total_bases_passing_filter, 2));
}

}  // namespace xoos::alignment_metrics
