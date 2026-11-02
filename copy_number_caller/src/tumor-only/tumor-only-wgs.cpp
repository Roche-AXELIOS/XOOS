#include "tumor-only/tumor-only-wgs.h"

#include <fstream>

#include <taskflow/taskflow.hpp>

#include <xoos/log/logging.h>
#include <xoos/sex_predict/sex.h>
#include <xoos/types/float.h>
#include <xoos/types/vec.h>
#include <xoos/util/file-functions.h>

#include "cli/copy-number-caller-modes.h"
#include "coverage/bam-to-corrected-coverage.h"
#include "coverage/coverage-check.h"
#include "coverage/coverage.h"
#include "intervals/check-chr-y-par.h"
#include "io/augmented-intervals-io.h"
#include "io/copy-number-caller-default-filenames.h"
#include "io/fai.h"
#include "io/write-empty-outputs.h"
#include "io/write-igv-xml.h"
#include "io/write-logr.h"
#include "io/write-metrics.h"
#include "io/write-segments.h"
#include "observations/observations.h"
#include "observations/vcf-check.h"
#include "one-sample-logr/self-normalize.h"
#include "seg-to-vcf/seg-to-vcf.h"
#include "seg-to-vcf/vcf-purity-source.h"
#include "segmentation/genomic-segments.h"
#include "segmentation/read-segments.h"
#include "segmentation/segment-and-call.h"
#include "segmentation/segment-type.h"
#include "tumor-only/tumor-only-bafs.h"
#include "util/mapq-util.h"

namespace xoos::cnc {

using xoos::file::PrefixedOutputPath;

struct SomaticTumorOnlyWGSData {
  Observations logrs;
  Observations mapqs;
  sex_predict::Sex sex = sex_predict::Sex::kUnknown;
  Observations ref_depths;
  Observations alt_depths;
  vec<GenomicSegment> segments;
  f64 tumor_purity = 0.99;
  f64 tumor_ploidy = 2;
  VcfPuritySource vcf_purity_source = VcfPuritySource::kNone;
  std::optional<TotalCopyNumberPrior> total_copy_number_prior;
  std::optional<CoverageCheckResult> coverage_check_result;
};

SomaticTumorOnlyWGSOutputPaths SetupDefaultSomaticTumorOnlyWGSOutputPaths(const CopyNumberCallerOptions& options) {
  const auto& dir = options.output_dir;
  const auto& prefix = options.output_prefix;
  SomaticTumorOnlyWGSOutputPaths paths{
      .likelihood_out = PrefixedOutputPath(dir, prefix, kDefaultSomaticCNCallsetSegOutput),
      .vcf_out = PrefixedOutputPath(dir, prefix, kDefaultSomaticCNCallsetVcfOutput),
      .metrics_out = PrefixedOutputPath(dir, prefix, kDefaultMetricsOutput),
      .igv_xml_out = PrefixedOutputPath(dir, prefix, kDefaultIgvXmlOutput),
      .logrs_out = PrefixedOutputPath(dir, prefix, kDefaultLogRsOutput),
      .logrs_bw_out = PrefixedOutputPath(dir, prefix, kDefaultLogRsBwOutput),
      .logr_segments_out = PrefixedOutputPath(dir, prefix, kDefaultLogRsSegOutput),
      .purity_ploidy_grid_out = PrefixedOutputPath(dir, prefix, kDefaultPurityPloidyGridOutput),
      .mapping_qualities_out = PrefixedOutputPath(dir, prefix, kDefaultMapQOutput),
      .augmented_intervals_out = PrefixedOutputPath(dir, prefix, kDefaultAugmentedIntervalsOutput),
      .tumor_coverage_out = PrefixedOutputPath(dir, prefix, kDefaultTumorCoverageOutput),
      .tumor_corrected_coverage_out = PrefixedOutputPath(dir, prefix, kDefaultTumorCorrectedCoverageOutput),
      .baf_out = PrefixedOutputPath(dir, prefix, kDefaultBafOutput),
      .baf_bw_out = PrefixedOutputPath(dir, prefix, kDefaultBafBwOutput),
      .taskflow_graph_out = std::nullopt,
  };
  if (options.save_baf_segments) {
    paths.baf_segments_out = PrefixedOutputPath(dir, prefix, kDefaultBafSegOutput);
  }
  if (options.save_taskflow_graph) {
    paths.taskflow_graph_out = PrefixedOutputPath(dir, prefix, kDefaultTaskflowGraphOutput);
  }
  vec<fs::path> output_paths = {
      paths.likelihood_out,
      paths.vcf_out,
      paths.metrics_out,
      paths.igv_xml_out,
      paths.logrs_out,
      paths.logrs_bw_out,
      paths.logr_segments_out,
      paths.purity_ploidy_grid_out,
      paths.mapping_qualities_out,
      paths.augmented_intervals_out,
      paths.tumor_coverage_out,
      paths.tumor_corrected_coverage_out,
      paths.baf_out,
      paths.baf_bw_out,
  };
  if (paths.baf_segments_out.has_value()) {
    output_paths.push_back(paths.baf_segments_out.value());
  }
  if (paths.taskflow_graph_out.has_value()) {
    output_paths.push_back(paths.taskflow_graph_out.value());
  }
  if (options.purity_ploidy_search_options.write_segment_likelihood) {
    output_paths.push_back(PrefixedOutputPath(dir, prefix, kDefaultSegmentLikelihoodDiploidBiasOutput));
    output_paths.push_back(PrefixedOutputPath(dir, prefix, kDefaultSegmentLikelihoodUniformOutput));
  }
  file::CheckOutputFilePathsForWriting(output_paths, options.overwrite);
  return paths;
}

namespace {
// tumor-only mode always runs the BAF-aware likelihood, so called segments carry this type
constexpr SegmentType kTumorOnlyOutputSegmentType = SegmentType::kSomaticTumorOnlyWithBafLikelihood;

/*
 * @brief Taskflow subroutine: single-sample tumor coverage, self-normalized logR, and sex.
 *
 * Single-sample analog of the tumor-normal CoverageAndLogRsTask: only the tumor BAM is
 * processed, and logR is produced by self-normalizing the tumor counts (as in
 * germline-normal-wgs) instead of dividing tumor by normal coverage.
 *
 * @param subflow Subflow
 * @param options Copy number caller options
 * @param paths Output file paths
 * @param augmented_intervals AugmentedIntervals for coverage calculation
 * @param data Shared working state to populate
 */
void CoverageAndLogRsTask(tf::Subflow& subflow,
                          const CopyNumberCallerOptions& options,
                          const SomaticTumorOnlyWGSOutputPaths& paths,
                          const AugmentedIntervals& augmented_intervals,
                          SomaticTumorOnlyWGSData& data) {
  Logging::Info("calculating and correcting coverage for {}", options.tumor_bam_fname.string());
  BamToCorrectedCoverageTaskFlowGraph tumor_corrected_coverage_tf_graph(
      options.tumor_bam_fname,
      augmented_intervals,
      paths.tumor_coverage_out,
      paths.tumor_corrected_coverage_out,
      options.calculate_coverage_options.exclude_flags,
      options.calculate_coverage_options.ignore_DN,
      options.gc_correct_options.first_span,
      options.command_line_info);
  subflow.composed_of(tumor_corrected_coverage_tf_graph).name("tumor coverage and correction");
  subflow.join();
  // low-coverage early exit
  if (tumor_corrected_coverage_tf_graph.IsLowCoverage()) {
    data.coverage_check_result = tumor_corrected_coverage_tf_graph.GetCoverageCheckResult();
    return;
  }
  const auto& tumor_coverage = tumor_corrected_coverage_tf_graph.GetResult();
  // no matched normal in tumor-only mode: predict sex from the tumor coverage unless overridden
  data.sex = sex_predict::ResolveSex(options.sample_metadata_options.sex_selection,
                                     [&tumor_coverage] { return tumor_coverage.PredictSex(); });
  // in tumor-only mode the tumor counts are normalized against themselves (median autosomal)
  Logging::Info("self-normalizing counts for {}", options.tumor_bam_fname.string());
  data.logrs = SelfNormalizeCounts(tumor_coverage);
  WriteLogRFiles(
      data.logrs, paths.logrs_out, paths.logrs_bw_out, options.reference_genome_fai_fname, options.command_line_info);
  data.mapqs = GetAvgMapqsFromCoverage(tumor_coverage);
  std::ofstream mapq_ofs(paths.mapping_qualities_out);
  data.mapqs.Write(mapq_ofs, kMeanMapqColumnName, true, options.command_line_info);
}

/*
 * @brief Taskflow subroutine: BAF consumption, segmentation, purity/ploidy search, and likelihood.
 *
 * Single-sample analog of the tumor-normal SegmentationAndCallingTask. The BAFs are extracted from the
 * tumor VCF against the population panel, so there is no matched normal and no
 * somatic-VAF purity fallback: purity/ploidy comes solely from the SCNA grid search. Both the search and
 * the likelihood run in tumor-only mode, so germline vs. somatic segments are classified from copy number.
 */
void SegmentationAndCallingTask(tf::Subflow& subflow,
                                const CopyNumberCallerOptions& options,
                                const SomaticTumorOnlyWGSOutputPaths& paths,
                                const vec<GenomicSegment>& seed_segments,
                                const RefAltObservations& ref_alt_depths,
                                SomaticTumorOnlyWGSData& data) {
  // Tumor-only low-purity handling: log a warning only. Unlike the tumor-normal path there is no matched
  // normal and thus no somatic-VAF purity fallback; purity/ploidy comes solely from the SCNA grid search.
  const auto on_low_purity = [&data] {
    Logging::Info(
        "Estimated tumor purity is quite low ({:.2f}). Copy number calls may be unreliable at low purity levels.",
        data.tumor_purity);
  };
  SegmentAndCall(subflow, options, paths, seed_segments, ref_alt_depths, data, on_low_purity);
}
}  // namespace

CopyNumberCallerRet SomaticTumorOnlyWGS(const AugmentedIntervals& augmented_intervals,
                                        const vec<GenomicSegment>& seed_segments,
                                        const RefAltObservations& ref_alt_depths,
                                        const CopyNumberCallerOptions& options,
                                        const SomaticTumorOnlyWGSOutputPaths& paths) {
  SomaticTumorOnlyWGSData data;
  tf::Taskflow taskflow;
  tf::Task coverage_and_logrs_task =
      taskflow
          .emplace([&augmented_intervals, &options, &paths, &data](tf::Subflow& subflow) {
            CoverageAndLogRsTask(subflow, options, paths, augmented_intervals, data);
          })
          .name("coverage and logrs");
  const tf::Task segmentation_and_likelihood_task =
      taskflow
          .emplace([&options, &paths, &seed_segments, &ref_alt_depths, &data](tf::Subflow& subflow) {
            if (IsLowCoverage(data.coverage_check_result)) {
              return;
            }
            SegmentationAndCallingTask(subflow, options, paths, seed_segments, ref_alt_depths, data);
          })
          .name("segmentation and likelihood");
  coverage_and_logrs_task.precede(segmentation_and_likelihood_task);
  if (paths.taskflow_graph_out.has_value()) {
    std::ofstream ofs(paths.taskflow_graph_out.value());
    taskflow.dump(ofs);
  }
  Logging::Info("using {} threads", options.threads);
  tf::Executor executor(options.threads);
  executor.run(taskflow).get();

  // low-coverage early exit
  if (IsLowCoverage(data.coverage_check_result)) {
    CopyNumberCallerRet ret{};
    ret.coverage_check_result = data.coverage_check_result;
    return ret;
  }

  // add metadata to all the segments
  for (auto& seg : data.segments) {
    seg.id = options.sample_metadata_options.sample_id;
    seg.sex = data.sex;
  }
  return {.ref_depths = data.ref_depths,
          .alt_depths = data.alt_depths,
          .segments = data.segments,
          .sex = data.sex,
          .tumor_purity = data.tumor_purity,
          .tumor_ploidy = data.tumor_ploidy,
          .vcf_purity_source = data.vcf_purity_source,
          .total_copy_number_prior = data.total_copy_number_prior,
          .coverage_check_result = data.coverage_check_result};
}

void SomaticTumorOnlyWGSMain(const CopyNumberCallerOptions& options) {
  if (options.mode != CopyNumberCallerModes::kSomaticTumorOnlyWGS) {
    Logging::Error("invalid CopyNumberCaller mode for SomaticTumorOnlyWGSMain");
    throw std::runtime_error("invalid CopyNumberCaller mode for SomaticTumorOnlyWGSMain");
  }
  const auto paths = SetupDefaultSomaticTumorOnlyWGSOutputPaths(options);
  AugmentedIntervals augmented_intervals = LoadOrGenerateWholeGenomeIntervals(options, paths.augmented_intervals_out);

  // Extract BAF observations from the single-sample tumor VCF, confirming het SNP sites against the
  // population panel. There is no matched normal, so somatic VAFs are not produced.
  const auto ref_alt_depths = ExtractTumorOnlyBafs(options);

  // Graceful exit: no variants remain after depth/BAF filtering
  if (ref_alt_depths.ref_obvs.contigs.empty() || ref_alt_depths.alt_obvs.contigs.empty()) {
    Logging::Warn("No variants found in VCF after applying filters. Writing empty outputs.");
    WriteSomaticTumorOnlyEmptyOutputs(paths, std::nullopt, augmented_intervals, options, ref_alt_depths.vcf_check);
    return;
  }

  // Graceful exit: low het-SNP fraction
  if (IsVcfInsufficient(ref_alt_depths.vcf_check)) {
    const auto& vc = *ref_alt_depths.vcf_check;
    Logging::Warn(
        "Only {:.1f}% of the evaluated variants are heterozygous germline variants, which is below the expected "
        "threshold of {:.0f}%. Writing empty outputs.",
        vc.het_snp_fraction.value_or(0.0) * 100.0,
        kDefaultMinHetSnpFraction * 100.0);
    WriteSomaticTumorOnlyEmptyOutputs(paths, std::nullopt, augmented_intervals, options, ref_alt_depths.vcf_check);
    return;
  }

  const vec<GenomicSegment> seed_segments = ReadSegments(options.seed_segments_fname, SegmentType::kSeed);
  // PAR masking is driven by the tumor BAM (no matched normal in tumor-only mode)
  augmented_intervals = RemovePARIfChrYPARUnmasked(options.tumor_bam_fname, augmented_intervals, seed_segments);

  const auto ret = SomaticTumorOnlyWGS(augmented_intervals, seed_segments, ref_alt_depths, options, paths);

  // Low-coverage early exit — write header-only outputs and return
  if (IsLowCoverage(ret.coverage_check_result)) {
    Logging::Warn("Coverage too low to continue; writing empty outputs.");
    WriteSomaticTumorOnlyEmptyOutputs(paths, ret.coverage_check_result, augmented_intervals, options);
    return;
  }

  WriteSegments(paths.likelihood_out, ret.segments, kTumorOnlyOutputSegmentType, options.command_line_info);

  bool has_high_extreme_baf_proportion = false;
  if (ret.tumor_purity.has_value()) {
    has_high_extreme_baf_proportion =
        CheckExtremeBafProportionAndWarn(ret.ref_depths, ret.alt_depths, ret.tumor_purity.value());
  }
  WriteSegmentsToVcf(paths.vcf_out,
                     ret.segments,
                     {.seq_lengths = augmented_intervals.GetSeqLengths(),
                      .seq_order = GetContigOrder(options.reference_genome_fai_fname),
                      .sample_id = options.sample_metadata_options.sample_id,
                      .reference_file = augmented_intervals.GetReferenceFile(),
                      .sex = ret.sex,
                      .mode = options.mode,
                      .purity = ret.tumor_purity,
                      .ploidy = ret.tumor_ploidy,
                      .vcf_purity_source = ret.vcf_purity_source,
                      .total_copy_number_prior = ret.total_copy_number_prior,
                      .command_line_info = options.command_line_info,
                      .has_high_extreme_baf_proportion = has_high_extreme_baf_proportion});
  WriteMetricsFile(paths.metrics_out,
                   augmented_intervals,
                   ret.ref_depths.contigs.size(),
                   ret.sex,
                   ret.segments,
                   ret.tumor_purity,
                   ret.tumor_ploidy,
                   options.mode,
                   options.command_line_info,
                   std::nullopt);
  WriteIGVXML(paths.igv_xml_out, std::make_optional(paths.logrs_bw_out), std::make_optional(paths.baf_bw_out));
}

}  // namespace xoos::cnc
