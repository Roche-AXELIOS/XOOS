#pragma once

#include <functional>

#include <taskflow/taskflow.hpp>

#include <xoos/log/logging.h>
#include <xoos/types/float.h>
#include <xoos/types/vec.h>

#include "cli/copy-number-caller-options.h"
#include "cli/sample-metadata-options.h"
#include "io/write-baf-values.h"
#include "io/write-segments.h"
#include "likelihood/likelihood-model.h"
#include "likelihood/likelihood-options.h"
#include "likelihood/likelihood.h"
#include "observations/filter-observations.h"
#include "observations/observations.h"
#include "purity-ploidy-search/purity-ploidy-search.h"
#include "seg-to-vcf/vcf-purity-source.h"
#include "segmentation/genomic-segments.h"
#include "segmentation/parent-specific-binary-segmentation-taskflow-graph.h"
#include "segmentation/segmentation-options.h"

namespace xoos::cnc {

/// Predicted tumor purity below which copy-number calls are flagged as potentially unreliable.
constexpr f64 kLowPredictedPurity = 0.3;

/// Callback invoked after the SCNA grid search when the estimated purity is below kLowPredictedPurity.
///
/// The tumor-normal workflow uses this hook to attempt a somatic-VAF-based purity re-estimation;
/// the tumor-only workflow typically uses this hook for a warning-only behavior (it has no matched-normal somatic
/// VAFs).

/// Shared segmentation, purity/ploidy search, and copy-number calling body for the somatic workflows.
///
/// Duck-typed on the paths (SomaticTumorNormalWGSOutputPaths / SomaticTumorOnlyWGSOutputPaths) and
/// data (SomaticTumorNormalWGSData / SomaticTumorOnlyWGSData) structs, which expose the same fields
/// accessed here. The single divergence between the tumor-normal and tumor-only paths — what to do
/// when the searched purity is low — is injected via @p on_low_purity.
template <typename PathsT, typename DataT>
void SegmentAndCall(tf::Subflow& subflow,
                    const CopyNumberCallerOptions& options,
                    const PathsT& paths,
                    const vec<GenomicSegment>& seed_segments,
                    const RefAltObservations& ref_alt_depths,
                    DataT& data,
                    const std::function<void()>& on_low_purity) {
  const SegmentationOptions segmentation_options = options.segmentation_options;
  SampleMetadataOptions sample_metadata_options = options.sample_metadata_options;
  sample_metadata_options.sex = data.sex;
  // filter the SNPs to only those that overlap LogR intervals. Non-intersecting SNPs are irrelevant to the algorithm
  data.ref_depths = GetOverlappingObservations(data.logrs, ref_alt_depths.ref_obvs);
  data.alt_depths = GetOverlappingObservations(data.logrs, ref_alt_depths.alt_obvs);
  WriteBafFiles(data.ref_depths,
                data.alt_depths,
                paths.baf_out,
                paths.baf_bw_out,
                options.reference_genome_fai_fname,
                options.command_line_info);
  // calculate DH values (decrease in heterozygosity) from depths. These will be used for 2nd round of segmentation
  const auto filtered_dh_vals = GetDhFromDepths(data.ref_depths, data.alt_depths);
  ParentSpecificBinarySegmentationTaskflowGraph segment_graph(
      data.logrs, filtered_dh_vals, seed_segments, segmentation_options, sample_metadata_options);
  subflow.composed_of(segment_graph);
  subflow.join();
  WriteLogRSegments(segment_graph.GetLogRSegments(),
                    paths.logr_segments_out,
                    data.logrs,
                    options.sample_metadata_options.sample_id,
                    data.sex,
                    options.command_line_info);
  data.segments = segment_graph.GetResult();
  for (auto& seg : data.segments) {
    seg.id = options.sample_metadata_options.sample_id;
    seg.sex = data.sex;
  }
  if (paths.baf_segments_out.has_value()) {
    WriteBafSegments(data.segments, paths.baf_segments_out.value(), options.command_line_info);
  }
  Logging::Info("found {} segments", data.segments.size());
  // Recompute mean_logr and num_obs from raw per-locus logRs before the purity/ploidy search. Pruning
  // replaces mean_logr on BAF-split segments with the cluster-weighted value.
  PopulateGenomicSegmentOptionalFields(data.segments, data.logrs);
  // determine tumor purity/ploidy: user-provided values take precedence, otherwise search the SCNA grid
  if (sample_metadata_options.purity.has_value() && sample_metadata_options.ploidy.has_value()) {
    data.tumor_purity = sample_metadata_options.purity.value();
    data.tumor_ploidy = sample_metadata_options.ploidy.value();
    data.vcf_purity_source = VcfPuritySource::kUserInput;
  } else {
    const auto purity_ploidy_out = PurityPloidySearchWithWGDQC(data.segments,
                                                               data.logrs,
                                                               data.ref_depths,
                                                               data.alt_depths,
                                                               options.purity_ploidy_search_options,
                                                               paths.purity_ploidy_grid_out,
                                                               options.output_prefix,
                                                               data.sex,
                                                               options.command_line_info);
    data.tumor_purity = purity_ploidy_out.purity;
    data.tumor_ploidy = purity_ploidy_out.ploidy;
    data.vcf_purity_source = VcfPuritySource::kScna;
    data.total_copy_number_prior = purity_ploidy_out.total_copy_number_prior;
    if (data.tumor_purity < kLowPredictedPurity && on_low_purity != nullptr) {
      on_low_purity();
    }
  }
  // call copy numbers
  const LikelihoodModel likelihood_model = LikelihoodModel::kSerialSummarized;
  const LikelihoodOptions likelihood_options = options.likelihood_options;
  sample_metadata_options.purity = data.tumor_purity;
  sample_metadata_options.ploidy = data.tumor_ploidy;
  Logging::Info("Calculating likelihood and calling copy numbers using tumor purity {:.2f} and ploidy {:.2f}",
                data.tumor_purity,
                data.tumor_ploidy);
  data.segments = CalculateLikelihoods(data.segments,
                                       data.logrs,
                                       data.ref_depths,
                                       data.alt_depths,
                                       data.mapqs,
                                       likelihood_model,
                                       likelihood_options,
                                       sample_metadata_options);
}

}  // namespace xoos::cnc
