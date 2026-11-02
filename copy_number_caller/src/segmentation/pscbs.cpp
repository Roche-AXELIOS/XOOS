#include "segmentation/pscbs.h"

#include <fstream>

#include <taskflow/algorithm/for_each.hpp>

#include <xoos/io/metadata-util.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/log/logging.h>
#include <xoos/sex_predict/sex.h>
#include <xoos/types/vec.h>
#include <xoos/util/file-functions.h>

#include "cli/copy-number-caller-options.h"
#include "cli/sample-metadata-options.h"
#include "io/copy-number-caller-default-filenames.h"
#include "io/write-baf-values.h"
#include "io/write-bigwig.h"
#include "io/write-segments.h"
#include "observations/filter-observations.h"
#include "observations/observations.h"
#include "segmentation/parent-specific-binary-segmentation-taskflow-graph.h"
#include "segmentation/read-segments.h"
#include "segmentation/segment-type.h"
#include "segmentation/segmentation-options.h"

namespace xoos::cnc {

namespace {
SegmentType GetSegmentType(const CopyNumberCallerOptions& options) {
  if (options.vcf_fname.has_value() || options.bafs_fname.has_value() || options.ad_fname.has_value()) {
    return SegmentType::kBaf;
  }
  return SegmentType::kLogROnly;
}

/**
 * @brief Write the BAF value and bigwig outputs shared by the VCF and BAFs input paths.
 */
void WriteBafOutputs(const Observations& ref_depths,
                     const Observations& alt_depths,
                     const CopyNumberCallerOptions& options) {
  const bool has_reference = !options.reference_genome_fai_fname.empty();
  WriteBafFiles(ref_depths,
                alt_depths,
                options.output_dir / kDefaultBafOutput,
                has_reference ? std::make_optional(options.output_dir / kDefaultBafBwOutput) : std::nullopt,
                has_reference ? std::make_optional(options.reference_genome_fai_fname) : std::nullopt,
                options.command_line_info);
}

/**
 * @brief Extract overlapping ref/alt depths from a VCF and write the BAF outputs.
 */
Observations ComputeDhValuesFromVcf(const Observations& logrs, const CopyNumberCallerOptions& options) {
  io::VcfReader vcf_reader(options.vcf_fname.value().string());
  const bool has_normal = options.sample_metadata_options.normal_sample_name.has_value();
  const bool has_tumor = options.sample_metadata_options.tumor_sample_name.has_value();
  Observations ref_depths;
  Observations alt_depths;
  if (has_normal && has_tumor) {
    auto ref_alts =
        GetDepthsFromVcf(vcf_reader, false, true, options.baf_filter_options, options.sample_metadata_options);
    ref_depths = std::move(ref_alts.ref_obvs);
    alt_depths = std::move(ref_alts.alt_obvs);
  } else if (has_normal) {
    auto ref_alts =
        GetDepthsFromVcf(vcf_reader, false, false, options.baf_filter_options, options.sample_metadata_options);
    ref_depths = std::move(ref_alts.ref_obvs);
    alt_depths = std::move(ref_alts.alt_obvs);
  } else if (has_tumor) {
    SampleMetadataOptions tumor_as_normal = options.sample_metadata_options;
    tumor_as_normal.normal_sample_name = options.sample_metadata_options.tumor_sample_name;
    BafFilterOptions tumor_filter{options.baf_filter_options};
    tumor_filter.normal_sample_min_depth = options.baf_filter_options.tumor_sample_min_depth;
    tumor_filter.normal_sample_min_baf = options.baf_filter_options.tumor_sample_min_baf;
    tumor_filter.normal_sample_max_baf = options.baf_filter_options.tumor_sample_max_baf;
    auto ref_alts = GetDepthsFromVcf(vcf_reader, false, false, tumor_filter, tumor_as_normal);
    ref_depths = std::move(ref_alts.ref_obvs);
    alt_depths = std::move(ref_alts.alt_obvs);
  }
  ref_depths = GetOverlappingObservations(logrs, ref_depths);
  alt_depths = GetOverlappingObservations(logrs, alt_depths);
  Observations dh_values = GetDhFromDepths(ref_depths, alt_depths);
  dh_values.SetBAFSegObvsStatus(true);
  WriteBafOutputs(ref_depths, alt_depths, options);
  return dh_values;
}

/**
 * @brief Compute DH (inverted BAF) values from whichever BAF source is configured.
 *
 * Returns empty Observations when no BAF source (VCF, BAFs, or allele-depth TSV) is provided.
 */
Observations ComputeDhValues(const Observations& logrs, const CopyNumberCallerOptions& options) {
  if (options.vcf_fname.has_value()) {
    return ComputeDhValuesFromVcf(logrs, options);
  }
  if (options.bafs_fname.has_value()) {
    const auto ref_alts = GetDepthsFromFile(options.bafs_fname.value());
    auto ref_depths = GetOverlappingObservations(logrs, ref_alts.ref_obvs);
    auto alt_depths = GetOverlappingObservations(logrs, ref_alts.alt_obvs);
    Observations dh_values = GetDhFromDepths(ref_depths, alt_depths);
    dh_values.SetBAFSegObvsStatus(true);
    WriteBafOutputs(ref_depths, alt_depths, options);
    return dh_values;
  }
  if (options.ad_fname.has_value()) {
    const auto ref_alts = GetDepthsFromFile(options.ad_fname.value());
    Observations dh_values = GetDhFromDepths(ref_alts.ref_obvs, ref_alts.alt_obvs);
    dh_values = GetOverlappingObservations(logrs, dh_values);
    dh_values.SetBAFSegObvsStatus(true);
    return dh_values;
  }
  return {};
}
}  // namespace

/**
 * @brief entry-function for the Parent-Specific-Binary-Segmentation algorithm
 * @param cli_opts contains file names etc to be parsed for input to ParentSpecificBinarySegmentation function
 */
void ParentSpecificBinarySegmentationMain(const CopyNumberCallerOptions& options) {
  const auto logr_segments_out = options.output_dir / kDefaultLogRsSegOutput;
  const auto seg_out = options.output_dir / kDefaultSegmentationSegOutput;
  const bool has_baf_source =
      options.vcf_fname.has_value() || options.bafs_fname.has_value() || options.ad_fname.has_value();
  if (!has_baf_source && !options.reference_genome_fai_fname.empty()) {
    Logging::Warn("--reference is ignored when no BAF source (--vcf, --bafs, or --allele-depth-tsv) is provided");
  }
  // BAF segments (before hierarchical pruning) are written whenever a BAF source is provided.
  const auto baf_segments_out = has_baf_source ? options.output_dir / kDefaultBafSegOutput : fs::path{};
  vec<fs::path> output_paths = {logr_segments_out, seg_out};
  if (has_baf_source) {
    output_paths.push_back(baf_segments_out);
  }
  if (options.vcf_fname.has_value() || options.bafs_fname.has_value()) {
    output_paths.push_back(options.output_dir / kDefaultBafOutput);
    if (!options.reference_genome_fai_fname.empty()) {
      output_paths.push_back(options.output_dir / kDefaultBafBwOutput);
    }
  }
  file::CheckOutputFilePathsForWriting(output_paths, options.overwrite);
  std::ifstream logrs_stream(options.logrs_fname.value());
  const Observations logrs = ReadObservations(logrs_stream);
  const Observations dh_values = ComputeDhValues(logrs, options);
  vec<GenomicSegment> seed_segments;
  const std::ifstream seed_segments_stream;
  const auto& command_line_info = options.command_line_info;
  if (!options.seed_segments_fname.empty()) {
    seed_segments = ReadSegments(options.seed_segments_fname, SegmentType::kSeed);
  } else {
    Logging::Error("Seed segments are required input!");
    throw std::runtime_error("missing seed segments");
  }
  vec<GenomicSegment> segments = ParentSpecificBinarySegmentation(logrs,
                                                                  dh_values,
                                                                  seed_segments,
                                                                  options.segmentation_options,
                                                                  options.sample_metadata_options,
                                                                  options.threads,
                                                                  logr_segments_out,
                                                                  baf_segments_out,
                                                                  command_line_info);
  for (auto& seg : segments) {
    seg.id = options.sample_metadata_options.sample_id;
    seg.sex = options.sample_metadata_options.sex;
  }
  const SegmentType mode = GetSegmentType(options);
  WriteSegments(seg_out, segments, mode, command_line_info);
}

/**
 * @brief implementation of the Parent-Specific-Binary-Segmentation algorithm
 * described in Olshen et. al 2011
 * (https://doi.org/10.1093/bioinformatics/btr329). Briefly, this algorithm
 * performs Circular Binary Segmentation (CBS) first on the tumor-normal logr
 * values, then feeds the resulting segments into a second CBS step based on
 * B-allele fractions at germline heterozygous variants
 * @param logrs tumor-normal(panel) log-ratios for the sample
 * @param dh_vals (optional) the dh (inverted BAF) values for the target
 * @param seed_segments  seed segments (if empty, will not be used). Common seeds are partitions dividing chromosomes
 by centromeres/telomeres
 * @param max_p maximum p-value for calling a segment significant (default: 0.01)
 * @param n_permutations  number of permutations for significance testing
 * @param min_obs_per_segment minimum number of observations for a segment
 * @param n_threads number of threads to use
 */
vec<GenomicSegment> ParentSpecificBinarySegmentation(const Observations& logrs,
                                                     const Observations& dh_vals,
                                                     const vec<GenomicSegment>& seed_segments,
                                                     const SegmentationOptions& segmentation_options,
                                                     const SampleMetadataOptions& sample_metadata_options,
                                                     const size_t n_threads,
                                                     const fs::path& logr_segments_out,
                                                     const fs::path& baf_segments_out,
                                                     const io::CommandLineInfo& command_line_info) {
  tf::Taskflow taskflow;
  ParentSpecificBinarySegmentationTaskflowGraph pscbs_graph(
      logrs, dh_vals, seed_segments, segmentation_options, sample_metadata_options);
  tf::Task pscbs_task = taskflow.composed_of(pscbs_graph);
  vec<GenomicSegment> res;
  const tf::Task pscbs_res_task = taskflow.emplace([&res,
                                                    &logrs,
                                                    &dh_vals,
                                                    &pscbs_graph,
                                                    &logr_segments_out,
                                                    &baf_segments_out,
                                                    &command_line_info,
                                                    &sample_metadata_options]() {
    res = pscbs_graph.GetResult();
    // repopulate segments with original array data
    PopulateGenomicSegmentOptionalFields(res, logrs);
    if (dh_vals.starts.size() != 0) {
      PopulateGenomicSegmentOptionalFields(res, dh_vals);
    }
    if (!logr_segments_out.empty()) {
      WriteLogRSegments(pscbs_graph.GetLogRSegments(),
                        logr_segments_out,
                        logrs,
                        sample_metadata_options.sample_id,
                        sample_metadata_options.sex.value_or(sex_predict::Sex::kUnknown),
                        command_line_info);
    }
    // Write the BAF segments produced before hierarchical pruning. Only relevant when a BAF
    // source was provided; otherwise the BAF segment list is empty.
    if (!baf_segments_out.empty()) {
      auto& baf_segments = pscbs_graph.GetBafSegments();
      for (auto& seg : baf_segments) {
        seg.id = sample_metadata_options.sample_id;
        seg.sex = sample_metadata_options.sex;
      }
      WriteBafSegments(baf_segments, baf_segments_out, command_line_info);
    }
  });
  pscbs_task.precede(pscbs_res_task);
  tf::Executor executor(n_threads);
  executor.run(taskflow).get();
  return res;
}

}  // namespace xoos::cnc
