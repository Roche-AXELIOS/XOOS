#include "metrics/read-metrics/read-metrics-te.h"

#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>

#include <htslib/sam.h>

#include <csv.hpp>

#include "io/alignment-reader.h"
#include "metrics/metrics-names.h"
#include "util/format-util.h"

namespace xoos::alignment_metrics {

TargetEnrichmentReadMetrics::TargetEnrichmentReadMetrics(const fs::path& bam_path)
    : total_mapped_alignments_in_bam(GetTotalMappedAlignmentsInBam(bam_path)) {
}

TargetEnrichmentReadMetrics& TargetEnrichmentReadMetrics::operator+=(const TargetEnrichmentReadMetrics& obj) {
  // Note: total_mapped_alignments_in_bam is not summed because it is directly obtained from the BAM index
  // and is not calculated from per-thread metrics. If the other metrics object has already updated the TE or unmapped
  // metrics, we don't want to double count them.
  this->total_mapped_alignments_in_bam =
      std::max(this->total_mapped_alignments_in_bam, obj.total_mapped_alignments_in_bam);
  this->on_target_alignments += obj.on_target_alignments;
  this->on_target_reads_passing_filter += obj.on_target_reads_passing_filter;
  return *this;
}

std::vector<std::string> TargetEnrichmentReadMetrics::GetHeaders() {
  return {kNameMetricName, kNameCount, kNameDenominator, kNamePercentage};
}

void TargetEnrichmentReadMetrics::WriteTsv(const fs::path& output_path, const io::Comments& comments) const {
  std::ofstream ofstream(output_path);
  auto writer = csv::make_tsv_writer_buffered(ofstream);
  io::WriteTsvComments(writer, comments);
  writer << GetHeaders();
  writer << std::make_tuple(
      kNameTotalMappedAlignmentsInBam, total_mapped_alignments_in_bam, kNotApplicable, kNotApplicable);
  writer << std::make_tuple(kNameOnTargetAlignments,
                            on_target_alignments,
                            kNameTotalMappedAlignmentsInBam,
                            ToPercentageWithPrecision(on_target_alignments, total_mapped_alignments_in_bam, 2));
  writer << std::make_tuple(
      kNameOnTargetReadsPassingFilter,
      on_target_reads_passing_filter,
      kNameTotalMappedAlignmentsInBam,
      ToPercentageWithPrecision(on_target_reads_passing_filter, total_mapped_alignments_in_bam, 2));
}

u64 GetTotalMappedAlignmentsInBam(const fs::path& bam_path) {
  static u64 total_mapped_alignments_in_bam = 0;
  static bool total_mapped_alignments_in_bam_calculated = false;
  if (!total_mapped_alignments_in_bam_calculated) {
    const AlignmentReader alignment_reader = OpenAlignmentFile(bam_path);
    for (s32 i = 0; i < sam_hdr_nref(alignment_reader.header.get()); ++i) {
      u64 mapped = 0;
      u64 unmapped = 0;
      if (hts_idx_get_stat(alignment_reader.idx.get(), i, &mapped, &unmapped) == 0) {
        total_mapped_alignments_in_bam += mapped;
      }
    }
    total_mapped_alignments_in_bam_calculated = true;
  }
  return total_mapped_alignments_in_bam;
}

}  // namespace xoos::alignment_metrics
