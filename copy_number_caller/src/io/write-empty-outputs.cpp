#include "io/write-empty-outputs.h"

#include <fstream>
#include <optional>

#include <xoos/io/metadata-util.h>
#include <xoos/sex_predict/sex.h>

#include "cli/copy-number-caller-options.h"
#include "coverage/coverage-check.h"
#include "coverage/coverage.h"
#include "germline-normal-wgs/germline-normal-wgs.h"
#include "intervals/augmented-intervals.h"
#include "io/column-names.h"
#include "io/fai.h"
#include "io/write-bigwig.h"
#include "io/write-igv-xml.h"
#include "io/write-metrics.h"
#include "io/write-segments.h"
#include "likelihood/likelihood.h"
#include "observations/observations.h"
#include "observations/vcf-check.h"
#include "purity-ploidy-search/purity-ploidy-search.h"
#include "seg-to-vcf/seg-to-vcf.h"
#include "seg-to-vcf/vcf-purity-source.h"
#include "segmentation/segment-type.h"
#include "somatic-tumor-normal-wgs/somatic-tumor-normal-wgs.h"
#include "tumor-only/tumor-only-wgs.h"
#include "util/mapq-util.h"

namespace xoos::cnc {

void WriteGermlineEmptyOutputs(const GermlineNormalWGSOutputPaths& paths,
                               const CoverageCheckResult& coverage_check,
                               const AugmentedIntervals& augmented_intervals,
                               const CopyNumberCallerOptions& options) {
  // SEG file — header only, no data rows
  WriteSegments(paths.likelihood_out, {}, SegmentType::kGermlineLikelihood, options.command_line_info);

  // VCF file — header only, no records
  WriteSegmentsToVcf(paths.vcf_out,
                     {},
                     {.seq_lengths = augmented_intervals.GetSeqLengths(),
                      .seq_order = GetContigOrder(options.reference_genome_fai_fname),
                      .sample_id = options.sample_metadata_options.sample_id,
                      .reference_file = augmented_intervals.GetReferenceFile(),
                      .sex = sex_predict::Sex::kUnknown,
                      .mode = options.mode,
                      .purity = kGermlinePurity,
                      .ploidy = kDiploidPloidy,
                      .command_line_info = options.command_line_info});

  // Metrics file
  WriteMetricsFile(paths.metrics_out,
                   augmented_intervals,
                   0,
                   sex_predict::Sex::kUnknown,
                   {},
                   std::nullopt,
                   std::nullopt,
                   options.mode,
                   options.command_line_info,
                   coverage_check);

  // Corrected coverage — header only
  {
    std::ofstream ofs(paths.normal_corrected_coverage_out);
    CoverageRecords{}.Write(ofs, options.command_line_info);
  }

  // Log ratios — header only
  {
    std::ofstream ofs(paths.logrs_out);
    Observations{}.Write(ofs, kColumnLogRatio, false, options.command_line_info);
  }

  // Log ratios bigwig
  WriteEmptyBigWig(paths.logrs_bw_out, options.reference_genome_fai_fname);

  // Mapping quality — header only
  {
    std::ofstream ofs(paths.mapping_qualities_out);
    Observations{}.Write(ofs, kMeanMapqColumnName, false, options.command_line_info);
  }

  // Log ratio segments (optional)
  if (paths.logr_segments_out.has_value()) {
    WriteSegments(paths.logr_segments_out.value(), {}, SegmentType::kLogROnly, options.command_line_info);
  }

  // BAF files (optional)
  if (paths.baf_out.has_value()) {
    std::ofstream ofs(paths.baf_out.value());
    Observations{}.Write(ofs, kColumnBAF, false, options.command_line_info);
  }
  if (paths.baf_bw_out.has_value()) {
    WriteEmptyBigWig(paths.baf_bw_out.value(), options.reference_genome_fai_fname);
  }

  // IGV XML
  WriteIGVXML(paths.igv_xml_out, std::make_optional(paths.logrs_bw_out), paths.baf_bw_out);
}

namespace {

/// Write the somatic empty-output set common to the tumor-normal and tumor-only workflows.
///
/// Duck-typed on the output-paths struct: both SomaticTumorNormalWGSOutputPaths and
/// SomaticTumorOnlyWGSOutputPaths expose the fields accessed here. The workflow-specific
/// likelihood SEG segment type is passed in; the tumor-normal path writes its extra
/// normal corrected-coverage file separately.
template <typename PathsT>
void WriteSomaticEmptyOutputsCommon(const PathsT& paths,
                                    const std::optional<CoverageCheckResult>& coverage_check,
                                    const AugmentedIntervals& augmented_intervals,
                                    const CopyNumberCallerOptions& options,
                                    const std::optional<VcfCheckResult>& vcf_check,
                                    SegmentType likelihood_segment_type) {
  // SEG file — header only
  WriteSegments(paths.likelihood_out, {}, likelihood_segment_type, options.command_line_info);

  // VCF file — header only
  WriteSegmentsToVcf(paths.vcf_out,
                     {},
                     {.seq_lengths = augmented_intervals.GetSeqLengths(),
                      .seq_order = GetContigOrder(options.reference_genome_fai_fname),
                      .sample_id = options.sample_metadata_options.sample_id,
                      .reference_file = augmented_intervals.GetReferenceFile(),
                      .sex = sex_predict::Sex::kUnknown,
                      .mode = options.mode,
                      .purity = 0.0,
                      .ploidy = kDiploidPloidy,
                      .vcf_purity_source = VcfPuritySource::kNone,
                      .total_copy_number_prior = std::nullopt,
                      .command_line_info = options.command_line_info,
                      .has_high_extreme_baf_proportion = false});

  // Metrics file
  WriteMetricsFile(paths.metrics_out,
                   augmented_intervals,
                   0,
                   sex_predict::Sex::kUnknown,
                   {},
                   std::nullopt,
                   std::nullopt,
                   options.mode,
                   options.command_line_info,
                   coverage_check,
                   vcf_check);

  // Tumor corrected coverage — header only
  {
    std::ofstream ofs(paths.tumor_corrected_coverage_out);
    CoverageRecords{}.Write(ofs, options.command_line_info);
  }

  // Log ratios — header only
  {
    std::ofstream ofs(paths.logrs_out);
    Observations{}.Write(ofs, kColumnLogRatio, false, options.command_line_info);
  }

  // Log ratios bigwig
  WriteEmptyBigWig(paths.logrs_bw_out, options.reference_genome_fai_fname);

  // Log ratio segments — header only
  WriteSegments(paths.logr_segments_out, {}, SegmentType::kLogROnly, options.command_line_info);

  // BAF segments (optional)
  if (paths.baf_segments_out.has_value()) {
    WriteSegments(paths.baf_segments_out.value(), {}, SegmentType::kBaf, options.command_line_info);
  }

  // Mapping quality — header only
  {
    std::ofstream ofs(paths.mapping_qualities_out);
    Observations{}.Write(ofs, kMeanMapqColumnName, false, options.command_line_info);
  }

  // Purity/ploidy grid — only when purity/ploidy are not user-supplied
  if (!options.sample_metadata_options.purity.has_value() || !options.sample_metadata_options.ploidy.has_value()) {
    std::ofstream ofs(paths.purity_ploidy_grid_out);
    if (!options.command_line_info.version.empty()) {
      io::WriteTsvMetadata(ofs, options.command_line_info);
    }
    ofs << kPurityPloidyGridOutPurity << "\t" << kPurityPloidyGridOutPloidy << "\t"
        << kPurityPloidyGridOutTotalJointLikelihood << "\t" << kPurityPloidyGridOutTotalLogRatioLikelihood << "\t"
        << kPurityPloidyGridOutTotalBAlleleFrequencyLikelihood << '\n';
  }

  // BAF files — header only
  {
    std::ofstream ofs(paths.baf_out);
    Observations{}.Write(ofs, kColumnBAF, false, options.command_line_info);
  }
  WriteEmptyBigWig(paths.baf_bw_out, options.reference_genome_fai_fname);

  // IGV XML
  WriteIGVXML(paths.igv_xml_out, std::make_optional(paths.logrs_bw_out), std::make_optional(paths.baf_bw_out));
}

}  // namespace

void WriteSomaticEmptyOutputs(const SomaticTumorNormalWGSOutputPaths& paths,
                              const std::optional<CoverageCheckResult>& coverage_check,
                              const AugmentedIntervals& augmented_intervals,
                              const CopyNumberCallerOptions& options,
                              const std::optional<VcfCheckResult>& vcf_check) {
  WriteSomaticEmptyOutputsCommon(
      paths, coverage_check, augmented_intervals, options, vcf_check, SegmentType::kSomaticWithBafLikelihood);

  // Normal corrected coverage — header only (tumor-normal only)
  std::ofstream ofs(paths.normal_corrected_coverage_out);
  CoverageRecords{}.Write(ofs, options.command_line_info);
}

void WriteSomaticTumorOnlyEmptyOutputs(const SomaticTumorOnlyWGSOutputPaths& paths,
                                       const std::optional<CoverageCheckResult>& coverage_check,
                                       const AugmentedIntervals& augmented_intervals,
                                       const CopyNumberCallerOptions& options,
                                       const std::optional<VcfCheckResult>& vcf_check) {
  WriteSomaticEmptyOutputsCommon(
      paths, coverage_check, augmented_intervals, options, vcf_check, SegmentType::kSomaticTumorOnlyWithBafLikelihood);
}

void WriteSomaticTumorOnlyEmptyOutputs(const SomaticTumorOnlyWGSOutputPaths& paths,
                                       const std::optional<CoverageCheckResult>& coverage_check,
                                       const AugmentedIntervals& augmented_intervals,
                                       const CopyNumberCallerOptions& options) {
  WriteSomaticTumorOnlyEmptyOutputs(paths, coverage_check, augmented_intervals, options, std::nullopt);
}

}  // namespace xoos::cnc
