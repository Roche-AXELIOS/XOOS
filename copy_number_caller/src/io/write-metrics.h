#pragma once
#include <optional>

#include <xoos/io/metadata-util.h>
#include <xoos/sex_predict/sex.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/vec.h>

#include "cli/copy-number-caller-modes.h"
#include "coverage/coverage-check.h"
#include "intervals/augmented-intervals.h"
#include "observations/vcf-check.h"
#include "segmentation/genomic-segments.h"

namespace xoos::cnc {
/**
 * @brief Writes summary metrics about the copy number run
 *
 * @param metrics_fname fs::path to the output metrics file.
 * @param augmented_intervals AugmentedIntervals object containing target regions.
 * @param n_snps number of SNPs used for analysis
 * @param sex Sample sex
 * @param segments Vector of GenomicSegment objects to analyze.
 * @param purity Sample purity. Only used for somatic modes.
 * @param ploidy Sample ploidy. Only used for somatic modes.
 * @param mode The mode of the CopyNumberCaller. Both somatic modes (tumor-normal and tumor-only) share the
 *        somatic metric layout; see IsSomaticLikeMode.
 * @param command_line_info Metadata to write at the top of the file.
 * @param coverage_check Optional coverage check result. When provided, appends
 *        "Median coverage" and "Percentage of low coverage windows" rows.
 * @param vcf_check Optional VCF check result. In tumor-normal mode, appends
 *        "Number of heterozygous germline SNPs", "Number of somatic variants",
 *        and "Het SNP fraction" rows when provided. Ignored in every other mode:
 *        these rows require a matched normal, and the tumor-only metrics
 *        contract excludes them.
 */
void WriteMetricsFile(const fs::path& metrics_fname,
                      const AugmentedIntervals& augmented_intervals,
                      size_t n_snps,
                      sex_predict::Sex sex,
                      const vec<GenomicSegment>& segments,
                      const std::optional<f64>& purity,
                      const std::optional<f64>& ploidy,
                      CopyNumberCallerModes mode,
                      const io::CommandLineInfo& command_line_info,
                      const std::optional<CoverageCheckResult>& coverage_check = std::nullopt,
                      const std::optional<VcfCheckResult>& vcf_check = std::nullopt);
}  // namespace xoos::cnc
