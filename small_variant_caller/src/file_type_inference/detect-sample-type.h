#pragma once

#include <cstddef>
#include <optional>

#include <xoos/types/float.h>
#include <xoos/types/fs.h>

#include "core/sample-type.h"
#include "file_type_inference/sample-type-metrics.h"

namespace xoos::svc {

/**
 * @brief Runtime options for sample-type detection.
 * @details Shared by the standalone detector and filter_variants integration. The AF and phred thresholds
 * control metric accumulation and final classification, while @p threads controls whether an indexed VCF
 * can be scanned by contig in parallel.
 */
struct DetectSampleTypeOptions {
  /// Maximum allele frequency for alternate allele observations counted as low-frequency errors.
  f64 max_af{kDefaultMaxAfForSampleTypeDetection};
  /// Maximum phred value classified as FFPE; larger values classify as cell-line.
  f64 max_ffpe_phred{kDefaultMaxFfpePhredForSampleTypeDetection};
  /// Requested number of VCF scanning threads.
  size_t threads{1};
};

/**
 * @brief Compute sample-type detection metrics from a VCF.
 * @details Scans the selected sample in the VCF, counts low-allele-frequency SNV observations as an
 * error estimate, and combines that count with the median DP and callable length to produce the phred
 * metric used for FFPE versus cell-line classification. If @p bed_file is provided, only records
 * overlapping the BED intervals contribute to the metric and the BED interval length is used as the
 * denominator length.
 * @param vcf_file Path to the input VCF or VCF.GZ file.
 * @param bed_file Optional BED file limiting the records and denominator length used for detection.
 * @param options Detection thresholds and threading options.
 * @return Metrics used to classify the sample type.
 * @throws xoos::error::Error if no usable observations are found or the target length is zero.
 */
SampleTypeMetrics ComputeSampleTypeMetrics(const fs::path& vcf_file,
                                           const std::optional<fs::path>& bed_file,
                                           const DetectSampleTypeOptions& options);

/**
 * @brief Classify a sample-type Phred quality score.
 * @details Scores at or below @p max_ffpe_phred are treated as FFPE; higher scores are treated as
 * cell-line. Callers use this function after computing metrics in either the standalone detector or the
 * filter_variants VCF pre-scan path.
 * @param phred Computed sample-type Phred quality score.
 * @param max_ffpe_phred Maximum score still classified as FFPE.
 * @return SampleType::kFfpe when @p phred is below or equal to the threshold, otherwise SampleType::kCellLine.
 */
SampleType ClassifySampleType(f64 phred, f64 max_ffpe_phred);

}  // namespace xoos::svc
