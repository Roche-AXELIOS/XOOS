#pragma once
#include <optional>

#include <xoos/io/metadata-util.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>

#include "coverage/coverage-check.h"
#include "germline-normal-wgs/germline-normal-wgs.h"
#include "observations/vcf-check.h"
#include "somatic-tumor-normal-wgs/somatic-tumor-normal-wgs.h"
#include "tumor-only/tumor-only-wgs.h"

namespace xoos::cnc {

/**
 * @brief Write all output files for a germline workflow early exit.
 *
 * Writes header-only versions of all expected output files so downstream
 * consumers don't fail on missing files. The coverage file is already written
 * before this is called.
 */
void WriteGermlineEmptyOutputs(const GermlineNormalWGSOutputPaths& paths,
                               const CoverageCheckResult& coverage_check,
                               const AugmentedIntervals& augmented_intervals,
                               const CopyNumberCallerOptions& options);

/**
 * @brief Write all output files for a somatic workflow early exit.
 *
 * Handles both low-coverage and VCF data-quality graceful exits.
 * When coverage_check is provided, coverage metrics are included in the metrics file.
 * When vcf_check is provided, VCF-related metrics are included in the metrics file.
 */
void WriteSomaticEmptyOutputs(const SomaticTumorNormalWGSOutputPaths& paths,
                              const std::optional<CoverageCheckResult>& coverage_check,
                              const AugmentedIntervals& augmented_intervals,
                              const CopyNumberCallerOptions& options,
                              const std::optional<VcfCheckResult>& vcf_check = std::nullopt);

/**
 * @brief Write all output files for a somatic tumor-only workflow early exit.
 *
 * Single-sample analog of WriteSomaticEmptyOutputs: writes header-only versions of every
 * tumor-only output (no matched-normal corrected coverage). Handles both the low-coverage
 * and VCF data-quality graceful exits.
 */
void WriteSomaticTumorOnlyEmptyOutputs(const SomaticTumorOnlyWGSOutputPaths& paths,
                                       const std::optional<CoverageCheckResult>& coverage_check,
                                       const AugmentedIntervals& augmented_intervals,
                                       const CopyNumberCallerOptions& options,
                                       const std::optional<VcfCheckResult>& vcf_check);

/**
 * @brief Overload for the low-coverage exit, where no VCF data-quality result is available.
 */
void WriteSomaticTumorOnlyEmptyOutputs(const SomaticTumorOnlyWGSOutputPaths& paths,
                                       const std::optional<CoverageCheckResult>& coverage_check,
                                       const AugmentedIntervals& augmented_intervals,
                                       const CopyNumberCallerOptions& options);

}  // namespace xoos::cnc
