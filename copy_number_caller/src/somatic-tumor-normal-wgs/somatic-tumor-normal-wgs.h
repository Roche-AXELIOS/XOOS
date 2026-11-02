#pragma once
#include <optional>

#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "cli/copy-number-caller-options.h"
#include "intervals/augmented-intervals.h"
#include "observations/observations.h"
#include "segmentation/genomic-segments.h"
#include "util/copy-number-caller-ret.h"

namespace xoos::cnc {

struct SomaticTumorNormalWGSOutputPaths {
  fs::path likelihood_out;
  fs::path vcf_out;
  fs::path metrics_out;
  fs::path igv_xml_out;
  fs::path logrs_out;
  fs::path logrs_bw_out;
  fs::path logr_segments_out;
  std::optional<fs::path> baf_segments_out;
  fs::path purity_ploidy_grid_out;
  fs::path mapping_qualities_out;
  std::optional<fs::path> taskflow_graph_out;
  fs::path augmented_intervals_out;
  fs::path tumor_coverage_out;
  fs::path tumor_corrected_coverage_out;
  fs::path normal_coverage_out;
  fs::path normal_corrected_coverage_out;
  fs::path baf_out;
  fs::path baf_bw_out;
};

SomaticTumorNormalWGSOutputPaths SetupDefaultSomaticTumorNormalWGSOutputPaths(const CopyNumberCallerOptions& options);

/*
 * @brief Main function for somatic tumor-normal WGS copy number calling. Meant to be called via command line
 * @param options Copy number caller options
 * @return CopyNumberCallerRet structure containing results
 */
void SomaticTumorNormalWGSMain(const CopyNumberCallerOptions& options);

/*
 * @brief Somatic tumor-normal WGS copy number calling function
 * @param augmented_intervals AugmentedIntervals for coverage calculation
 * @param seed_segments Seed segments for segmentation
 * @param ref_depths Reference allele depths from VCF
 * @param alt_depths Alternate allele depths from VCF
 * @param options Copy number caller options
 * @param paths Output file paths
 * @return CopyNumberCallerRet structure containing results
 */
CopyNumberCallerRet SomaticTumorNormalWGS(const AugmentedIntervals& augmented_intervals,
                                          const vec<GenomicSegment>& seed_segments,
                                          const RefAltObservations& ref_alt_depths,
                                          const CopyNumberCallerOptions& options,
                                          const SomaticTumorNormalWGSOutputPaths& paths);
}  // namespace xoos::cnc
