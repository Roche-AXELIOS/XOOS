#pragma once
#include <optional>

#include <xoos/types/fs.h>
#include <xoos/types/vec.h>

#include "cli/copy-number-caller-options.h"
#include "intervals/augmented-intervals.h"
#include "observations/observations.h"
#include "segmentation/genomic-segments.h"
#include "util/copy-number-caller-ret.h"

namespace xoos::cnc {

struct SomaticTumorOnlyWGSOutputPaths {
  fs::path likelihood_out;
  fs::path vcf_out;
  fs::path metrics_out;
  fs::path igv_xml_out;
  fs::path logrs_out;
  fs::path logrs_bw_out;
  fs::path logr_segments_out;
  fs::path purity_ploidy_grid_out;
  fs::path mapping_qualities_out;
  fs::path augmented_intervals_out;
  fs::path tumor_coverage_out;
  fs::path tumor_corrected_coverage_out;
  fs::path baf_out;
  fs::path baf_bw_out;
  std::optional<fs::path> baf_segments_out;
  std::optional<fs::path> taskflow_graph_out;
};

SomaticTumorOnlyWGSOutputPaths SetupDefaultSomaticTumorOnlyWGSOutputPaths(const CopyNumberCallerOptions& options);

/*
 * @brief Somatic tumor-only WGS copy number calling function.
 *
 * Runs the tumor-only workflow and returns the shared CopyNumberCallerRet. On the
 * low-coverage early exit, only CopyNumberCallerRet::coverage_check_result is populated.
 *
 * @param augmented_intervals AugmentedIntervals for coverage calculation
 * @param seed_segments Seed segments for segmentation
 * @param ref_alt_depths Ref/alt allele depths extracted from the tumor VCF (BAF observations)
 * @param options Copy number caller options
 * @param paths Output file paths
 * @return CopyNumberCallerRet structure containing results
 */
CopyNumberCallerRet SomaticTumorOnlyWGS(const AugmentedIntervals& augmented_intervals,
                                        const vec<GenomicSegment>& seed_segments,
                                        const RefAltObservations& ref_alt_depths,
                                        const CopyNumberCallerOptions& options,
                                        const SomaticTumorOnlyWGSOutputPaths& paths);

/*
 * @brief Main function for somatic tumor-only WGS copy number calling. Meant to be called via command line.
 *
 * @param options Copy number caller options
 */
void SomaticTumorOnlyWGSMain(const CopyNumberCallerOptions& options);

}  // namespace xoos::cnc
