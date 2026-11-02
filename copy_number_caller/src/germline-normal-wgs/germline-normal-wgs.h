#pragma once
#include <optional>

#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "cli/copy-number-caller-options.h"
#include "intervals/augmented-intervals.h"
#include "likelihood/likelihood-options.h"
#include "observations/observations.h"
#include "segmentation/genomic-segments.h"
#include "util/copy-number-caller-ret.h"

namespace xoos::cnc {

struct GermlineNormalWGSOutputPaths {
  fs::path likelihood_out;
  fs::path vcf_out;
  fs::path metrics_out;
  fs::path igv_xml_out;
  fs::path logrs_out;
  fs::path logrs_bw_out;
  std::optional<fs::path> logr_segments_out;
  fs::path mapping_qualities_out;
  std::optional<fs::path> taskflow_graph_out;
  fs::path augmented_intervals_out;
  fs::path normal_coverage_out;
  fs::path normal_corrected_coverage_out;
  std::optional<fs::path> baf_out;
  std::optional<fs::path> baf_bw_out;
};

GermlineNormalWGSOutputPaths SetupDefaultGermlineNormalWGSOutputPaths(const CopyNumberCallerOptions& options);

/*
 * @brief Override mapq_cutoff_for_calls to 0 when rescued secondary alignments (YF:i:1) were
 *        detected during coverage calculation. The override is skipped when the user explicitly
 *        set --min-mapq-for-calls on the command line, or when the cutoff is already <= 0.
 * @param likelihood_options Likelihood options to (possibly) modify in-place
 * @param has_rescued_secondaries Whether rescued secondary alignments were detected in the BAM
 */
void ApplyRescuedSecondaryMapqOverride(LikelihoodOptions& likelihood_options, bool has_rescued_secondaries);

/*
 * @brief Main function for germline normal WGS copy number calling. Meant to be called via command line
 * @param options Copy number caller options
 * @return CopyNumberCallerRet structure containing results
 */
void GermlineNormalWGSMain(const CopyNumberCallerOptions& options);
/*
 * @brief Germline normal WGS copy number calling function
 * @param augmented_intervals AugmentedIntervals for coverage calculation
 * @param seed_segments Seed segments for segmentation
 * @param ref_obvs Optional reference allele observations from VCF
 * @param alt_obvs Optional alternate allele observations from VCF
 * @param options Copy number caller options
 * @param paths Output file paths
 * @return CopyNumberCallerRet structure containing results
 */
CopyNumberCallerRet GermlineNormalWGS(const AugmentedIntervals& augmented_intervals,
                                      const vec<GenomicSegment>& seed_segments,
                                      const std::optional<Observations>& ref_obvs,
                                      const std::optional<Observations>& alt_obvs,
                                      const CopyNumberCallerOptions& options,
                                      const GermlineNormalWGSOutputPaths& paths);
}  // namespace xoos::cnc
