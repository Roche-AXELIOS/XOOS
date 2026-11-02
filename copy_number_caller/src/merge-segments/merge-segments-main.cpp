#include "merge-segments/merge-segments-main.h"

#include <fstream>
#include <stdexcept>

#include <xoos/log/logging.h>
#include <xoos/types/vec.h>
#include <xoos/util/file-functions.h>

#include "cli/copy-number-caller-options.h"
#include "io/copy-number-caller-default-filenames.h"
#include "io/write-segments.h"
#include "likelihood/flag-copy-number-calls.h"
#include "merge-segments/merge-segments.h"
#include "observations/observations.h"
#include "seg-to-vcf/seg-to-vcf.h"
#include "segmentation/genomic-segments.h"
#include "segmentation/read-segments.h"
#include "segmentation/segment-type.h"
#include "segmentation/segments-header.h"
#include "util/mapq-util.h"

namespace xoos::cnc {

void MergeSegmentsMain(const CopyNumberCallerOptions& options) {
  const auto vcf_out_fname = options.output_dir / kDefaultGermlineCNCallsetVcfOutput;
  const auto merged_segments_out_fname = options.output_dir / kDefaultGermlineCNCallsetSegOutput;
  file::CheckOutputFilePathsForWriting({merged_segments_out_fname, vcf_out_fname}, options.overwrite);
  vec<GenomicSegment> segments = ReadSegments(options.segments_fname.value(), SegmentType::kGermlineLikelihood);
  const SegmentsHeader header = ReadHeaderFromSegments(options.segments_fname.value());
  Observations mapqs;
  if (options.merge_segments_options.recalculate_per_segment_data ||
      options.merge_segments_options.use_mapqs_observations) {
    if (!options.mapping_qualities_fname.has_value()) {
      Logging::Error("Must provide --mapqs option!");
      throw std::runtime_error("missing option");
    }
    std::ifstream mapqs_ifs(options.mapping_qualities_fname.value());
    mapqs = ReadObservations(mapqs_ifs);
    if (options.merge_segments_options.use_mapqs_observations) {
      AssignAvgMeanMapqPerSegment(segments, mapqs);
    }
  }
  auto merged_segments = MergeAdjacentEqualCopyNumberSegments(segments);
  merged_segments = MergeLowMapqSegments(merged_segments, options.merge_segments_options.min_mapq_threshold);
  if (options.merge_segments_options.recalculate_per_segment_data) {
    if (!options.logrs_fname.has_value()) {
      Logging::Error("Must provide --logrs with the --recalculate-per-segment-data option!");
      throw std::runtime_error("missing option");
    }
    std::ifstream logrs_ifs(options.logrs_fname.value());
    const Observations logrs = ReadObservations(logrs_ifs);
    AssignAvgMeanMapqPerSegment(merged_segments, mapqs);
    PopulateGenomicSegmentOptionalFields(merged_segments, logrs);
    FlagCallsByExpectedTotalCopyNumber(merged_segments,
                                       options.sample_metadata_options.sex.value(),
                                       options.likelihood_options.non_integer_tcn_filter_options);
  }
  FlagCallsByMeanMapq(merged_segments, options.merge_segments_options.min_mapq_threshold);
  for (auto& seg : merged_segments) {
    if (!seg.flags.has_value()) {
      seg.flags.emplace(vec<std::string>{});
    }
  }
  WriteSegments(
      merged_segments_out_fname, merged_segments, SegmentType::kGermlineLikelihood, options.command_line_info);
  WriteSegmentsToVcf(vcf_out_fname,
                     merged_segments,
                     {.seq_lengths = {},
                      .seq_order = {},
                      .sample_id = options.sample_metadata_options.sample_id,
                      .reference_file = header.reference_file,
                      .sex = options.sample_metadata_options.sex.value(),
                      .command_line_info = options.command_line_info});
}

}  // namespace xoos::cnc
