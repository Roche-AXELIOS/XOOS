#include "germline-normal-wgs/germline-normal-wgs.h"

#include <fstream>

#include <taskflow/algorithm/for_each.hpp>

#include <xoos/io/metadata-util.h>
#include <xoos/log/logging.h>
#include <xoos/types/vec.h>
#include <xoos/util/file-functions.h>

#include "coverage/bam-to-corrected-coverage.h"
#include "coverage/coverage-check.h"
#include "intervals/check-chr-y-par.h"
#include "io/augmented-intervals-io.h"
#include "io/copy-number-caller-default-filenames.h"
#include "io/write-baf-values.h"
#include "io/write-empty-outputs.h"
#include "io/write-igv-xml.h"
#include "io/write-logr.h"
#include "io/write-metrics.h"
#include "io/write-segments.h"
#include "likelihood/flag-copy-number-calls.h"
#include "likelihood/likelihood.h"
#include "merge-segments/merge-segments.h"
#include "observations/filter-observations.h"
#include "observations/observations.h"
#include "one-sample-logr/self-normalize.h"
#include "seg-to-vcf/seg-to-vcf.h"
#include "segmentation/genomic-segments.h"
#include "segmentation/parent-specific-binary-segmentation-taskflow-graph.h"
#include "segmentation/read-segments.h"
#include "util/copy-number-caller-ret.h"
#include "util/mapq-util.h"

namespace xoos::cnc {

using xoos::file::PrefixedOutputPath;

GermlineNormalWGSOutputPaths SetupDefaultGermlineNormalWGSOutputPaths(const CopyNumberCallerOptions& options) {
  const auto& dir = options.output_dir;
  const auto& prefix = options.output_prefix;
  GermlineNormalWGSOutputPaths paths{
      .likelihood_out = PrefixedOutputPath(dir, prefix, kDefaultGermlineCNCallsetSegOutput),
      .vcf_out = PrefixedOutputPath(dir, prefix, kDefaultGermlineCNCallsetVcfOutput),
      .metrics_out = PrefixedOutputPath(dir, prefix, kDefaultMetricsOutput),
      .igv_xml_out = PrefixedOutputPath(dir, prefix, kDefaultIgvXmlOutput),
      .logrs_out = PrefixedOutputPath(dir, prefix, kDefaultLogRsOutput),
      .logrs_bw_out = PrefixedOutputPath(dir, prefix, kDefaultLogRsBwOutput),
      .mapping_qualities_out = PrefixedOutputPath(dir, prefix, kDefaultMapQOutput),
      .augmented_intervals_out = PrefixedOutputPath(dir, prefix, kDefaultAugmentedIntervalsOutput),
      .normal_coverage_out = PrefixedOutputPath(dir, prefix, kDefaultCoverageOutput),
      .normal_corrected_coverage_out = PrefixedOutputPath(dir, prefix, kDefaultCorrectedCoverageOutput),
  };
  if (options.vcf_fname.has_value()) {
    paths.baf_out = PrefixedOutputPath(dir, prefix, kDefaultBafOutput);
    paths.baf_bw_out = PrefixedOutputPath(dir, prefix, kDefaultBafBwOutput);
  }
  if (options.save_log_ratio_segments) {
    paths.logr_segments_out = PrefixedOutputPath(dir, prefix, kDefaultLogRsSegOutput);
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
      paths.mapping_qualities_out,
      paths.augmented_intervals_out,
      paths.normal_coverage_out,
      paths.normal_corrected_coverage_out,
  };
  if (paths.logr_segments_out.has_value()) {
    output_paths.push_back(paths.logr_segments_out.value());
  }
  if (paths.taskflow_graph_out.has_value()) {
    output_paths.push_back(paths.taskflow_graph_out.value());
  }
  if (paths.baf_out.has_value()) {
    output_paths.push_back(paths.baf_out.value());
  }
  if (paths.baf_bw_out.has_value()) {
    output_paths.push_back(paths.baf_bw_out.value());
  }
  file::CheckOutputFilePathsForWriting(output_paths, options.overwrite);
  return paths;
}

void ApplyRescuedSecondaryMapqOverride(LikelihoodOptions& likelihood_options, const bool has_rescued_secondaries) {
  if (!has_rescued_secondaries) {
    return;
  }
  if (likelihood_options.mapq_cutoff_for_calls_is_user_set) {
    Logging::Info(
        "Rescued secondary alignments (YF:i:1) detected during coverage calculation, "
        "but --min-mapq-for-calls was explicitly set by the user to {}. Not overriding.",
        likelihood_options.mapq_cutoff_for_calls);
    return;
  }
  if (likelihood_options.mapq_cutoff_for_calls > 0) {
    Logging::Info(
        "Rescued secondary alignments (YF:i:1) detected during coverage calculation. "
        "Overriding --min-mapq-for-calls to 0 to allow common duplications to be predicted.");
    likelihood_options.mapq_cutoff_for_calls = 0;
  }
}

void GermlineNormalWGSMain(const CopyNumberCallerOptions& options) {
  const auto paths = SetupDefaultGermlineNormalWGSOutputPaths(options);
  AugmentedIntervals augmented_intervals = LoadOrGenerateWholeGenomeIntervals(options, paths.augmented_intervals_out);
  std::optional<Observations> ref_obvs{};
  std::optional<Observations> alt_obvs{};
  // load ref/alt observation if user desires to print out BAFs. note that BAFs are not currently used in the actual
  // calling algorithm for germline
  if (options.vcf_fname.has_value()) {
    io::VcfReader vcf_reader_for_baf(options.vcf_fname.value());
    auto ref_alt_depths =
        GetDepthsFromVcf(vcf_reader_for_baf, true, false, options.baf_filter_options, options.sample_metadata_options);
    ref_obvs = std::move(ref_alt_depths.ref_obvs);
    alt_obvs = std::move(ref_alt_depths.alt_obvs);
    if (ref_obvs->contigs.empty() || alt_obvs->contigs.empty()) {
      Logging::Warn("No variants found in VCF after applying filters");
    }
  }
  const vec<GenomicSegment> seed_segments = ReadSegments(options.seed_segments_fname, SegmentType::kSeed);
  augmented_intervals = RemovePARIfChrYPARUnmasked(options.normal_bam_fname, augmented_intervals, seed_segments);
  auto ret = GermlineNormalWGS(augmented_intervals, seed_segments, ref_obvs, alt_obvs, options, paths);

  // Handle low coverage early exit — write header-only outputs and return
  if (IsLowCoverage(ret.coverage_check_result)) {
    WriteGermlineEmptyOutputs(paths, *ret.coverage_check_result, augmented_intervals, options);
    return;
  }

  WriteSegments(paths.likelihood_out, ret.segments, SegmentType::kGermlineLikelihood, options.command_line_info);
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
                      .command_line_info = options.command_line_info});
  WriteMetricsFile(paths.metrics_out,
                   augmented_intervals,
                   0,
                   ret.sex,
                   ret.segments,
                   ret.tumor_purity,
                   ret.tumor_ploidy,
                   options.mode,
                   options.command_line_info,
                   std::nullopt);
  WriteIGVXML(paths.igv_xml_out, std::make_optional(paths.logrs_bw_out), paths.baf_bw_out);
}

struct GermlineNormalWGSData {
  Observations logrs;
  Observations mapqs;
  sex_predict::Sex sex = sex_predict::Sex::kUnknown;
  vec<GenomicSegment> segments;
  std::optional<Observations> filtered_ref_obvs;
  std::optional<Observations> filtered_alt_obvs;
  std::optional<CoverageCheckResult> coverage_check_result;
  bool has_rescued_secondaries = false;
};

/*
 * @brief Taskflow subroutine to calculate coverage and logR values
 * @param subflow Subflow
 * @param options Copy number caller options
 * @param augmented_intervals AugmentedIntervals for coverage calculation
 * @param data Data structure to hold results
 */
namespace {
void CoverageAndLogRsTask(tf::Subflow& subflow,
                          const CopyNumberCallerOptions& options,
                          const GermlineNormalWGSOutputPaths& paths,
                          const AugmentedIntervals& augmented_intervals,
                          GermlineNormalWGSData& data) {
  Logging::Info("calculating and correcting coverage for {}", options.normal_bam_fname.string());
  BamToCorrectedCoverageTaskFlowGraph normal_corrected_coverage_tf_graph(
      options.normal_bam_fname,
      augmented_intervals,
      paths.normal_coverage_out,
      paths.normal_corrected_coverage_out,
      options.calculate_coverage_options.exclude_flags,
      options.calculate_coverage_options.ignore_DN,
      options.gc_correct_options.first_span,
      options.command_line_info);
  subflow.composed_of(normal_corrected_coverage_tf_graph).name("normal coverage and correction");
  subflow.join();
  // Check for low coverage early exit
  if (normal_corrected_coverage_tf_graph.IsLowCoverage()) {
    data.coverage_check_result = normal_corrected_coverage_tf_graph.GetCoverageCheckResult();
    return;
  }
  const auto& normal_coverage = normal_corrected_coverage_tf_graph.GetResult();
  data.has_rescued_secondaries = normal_coverage.has_rescued_secondaries.load();
  // assign sex
  data.sex = sex_predict::ResolveSex(options.sample_metadata_options.sex_selection,
                                     [&normal_coverage] { return normal_coverage.PredictSex(); });
  // in Germline Mode, the counts are normalized against themselves
  Logging::Info("self-normalizing counts for {}", options.normal_bam_fname.string());
  data.logrs = SelfNormalizeCounts(normal_coverage);
  // populate other return values and write to disk if necessary
  WriteLogRFiles(
      data.logrs, paths.logrs_out, paths.logrs_bw_out, options.reference_genome_fai_fname, options.command_line_info);
  data.mapqs = GetAvgMapqsFromCoverage(normal_coverage);
  std::ofstream mapq_ofs(paths.mapping_qualities_out);
  data.mapqs.Write(mapq_ofs, kMeanMapqColumnName, true, options.command_line_info);
}
}  // namespace

/*
 * @brief Taskflow subroutine to segment and call copy number
 * @param subflow Subflow
 * @param options Copy number caller options
 * @param augmented_intervals AugmentedIntervals for coverage calculation
 * @param seed_segments Seed segments for segmentation
 * @param data Data structure to hold results
 */
namespace {
void SegmentationAndCallingTask(tf::Subflow& subflow,
                                const CopyNumberCallerOptions& options,
                                const GermlineNormalWGSOutputPaths& paths,
                                const AugmentedIntervals& /*augmented_intervals*/,
                                const vec<GenomicSegment>& seed_segments,
                                GermlineNormalWGSData& data) {
  // segmentation
  const SegmentationOptions segmentation_options = options.segmentation_options;
  SampleMetadataOptions sample_metadata_options = options.sample_metadata_options;
  sample_metadata_options.sex = data.sex;
  ParentSpecificBinarySegmentationTaskflowGraph segment_graph(
      data.logrs, {}, seed_segments, segmentation_options, sample_metadata_options);
  subflow.composed_of(segment_graph);
  subflow.join();
  if (paths.logr_segments_out.has_value()) {
    WriteLogRSegments(segment_graph.GetLogRSegments(),
                      paths.logr_segments_out.value(),
                      data.logrs,
                      options.sample_metadata_options.sample_id,
                      data.sex,
                      options.command_line_info);
  }
  data.segments = segment_graph.GetResult();
  for (auto& seg : data.segments) {
    seg.id = options.sample_metadata_options.sample_id;
    seg.sex = data.sex;
  }
  Logging::Info("found {} segments", data.segments.size());
  // call copy numbers for each segment
  const LikelihoodModel likelihood_model = LikelihoodModel::kLogrSummarizedOnly;
  // make sure that alternative model is only used when in germline mode and when we have BAF observations
  LikelihoodOptions likelihood_options = options.likelihood_options;
  // If rescued secondary alignments were detected during coverage calculation, disable MAPQ-based
  // call flagging so that common duplications represented in the pangenome graph are not falsely
  // flagged as deletions. Only overrides when the user did not explicitly set --min-mapq-for-calls.
  ApplyRescuedSecondaryMapqOverride(likelihood_options, data.has_rescued_secondaries);
  // force purity/ploidy to be 0.99/2 when in germline mode
  sample_metadata_options.purity = kGermlinePurity;
  sample_metadata_options.ploidy = kDiploidPloidy;
  data.segments = CalculateLikelihoods(
      data.segments, data.logrs, {}, {}, data.mapqs, likelihood_model, likelihood_options, sample_metadata_options);
  // merge segments if necessary after calling copy numbers
  data.segments = MergeAdjacentEqualCopyNumberSegments(data.segments);
  // merge segments with low mean mapping quality
  data.segments = MergeLowMapqSegments(data.segments, likelihood_options.mapq_cutoff_for_calls);
  // flag potentially bad calls
  for (auto& seg : data.segments) {
    if (!seg.flags.has_value()) {
      seg.flags.emplace(vec<std::string>{});
    }
  }
  AssignAvgMeanMapqPerSegment(data.segments, data.mapqs);
  PopulateGenomicSegmentOptionalFields(data.segments, data.logrs);
  FlagCallsByMeanMapq(data.segments, likelihood_options.mapq_cutoff_for_calls);
  FlagCallsByExpectedTotalCopyNumber(data.segments, data.sex, likelihood_options.non_integer_tcn_filter_options);
}
}  // namespace

CopyNumberCallerRet GermlineNormalWGS(const AugmentedIntervals& augmented_intervals,
                                      const vec<GenomicSegment>& seed_segments,
                                      const std::optional<Observations>& ref_obvs,
                                      const std::optional<Observations>& alt_obvs,
                                      const CopyNumberCallerOptions& options,
                                      const GermlineNormalWGSOutputPaths& paths) {
  // define tasks
  GermlineNormalWGSData data;
  tf::Taskflow taskflow;
  tf::Task coverage_and_logrs_task =
      taskflow
          .emplace([&augmented_intervals, &options, &paths, &data](tf::Subflow& subflow) {
            CoverageAndLogRsTask(subflow, options, paths, augmented_intervals, data);
          })
          .name("coverage and logrs");
  // segment — skip if low coverage
  const tf::Task segmentation_and_likelihood_task =
      taskflow.emplace([&options, &paths, &augmented_intervals, &seed_segments, &data](tf::Subflow& subflow) {
        if (IsLowCoverage(data.coverage_check_result)) {
          return;
        }
        SegmentationAndCallingTask(subflow, options, paths, augmented_intervals, seed_segments, data);
      });
  // BAFs — skip if low coverage
  const tf::Task write_bafs_task = taskflow.emplace([&options, &paths, &ref_obvs, &alt_obvs, &data] {
    if (IsLowCoverage(data.coverage_check_result)) {
      return;
    }
    if (ref_obvs.has_value() && alt_obvs.has_value()) {
      data.filtered_ref_obvs = GetOverlappingObservations(data.logrs, ref_obvs.value());
      data.filtered_alt_obvs = GetOverlappingObservations(data.logrs, alt_obvs.value());
      WriteBafFiles(data.filtered_ref_obvs.value(),
                    data.filtered_alt_obvs.value(),
                    paths.baf_out,
                    paths.baf_bw_out,
                    options.reference_genome_fai_fname,
                    options.command_line_info);
    }
  });
  // define taskflow DAG
  coverage_and_logrs_task.precede(segmentation_and_likelihood_task);
  coverage_and_logrs_task.precede(write_bafs_task);
  if (paths.taskflow_graph_out.has_value()) {
    std::ofstream ofs(paths.taskflow_graph_out.value());
    taskflow.dump(ofs);
  }
  Logging::Info("using {} threads", options.threads);
  tf::Executor executor(options.threads);
  executor.run(taskflow).get();

  // Handle low coverage early exit
  if (IsLowCoverage(data.coverage_check_result)) {
    CopyNumberCallerRet ret{};
    ret.coverage_check_result = data.coverage_check_result;
    return ret;
  }

  // add some metadata to all the segments
  for (auto& seg : data.segments) {
    seg.id = options.sample_metadata_options.sample_id;
    seg.sex = data.sex;
  }
  return {.ref_depths = {},
          .alt_depths = {},
          .segments = data.segments,
          .sex = data.sex,
          .tumor_purity = kGermlinePurity,
          .tumor_ploidy = kDiploidPloidy};
}
}  // namespace xoos::cnc
