#pragma once

#include <xoos/io/vcf/vcf-header.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "util/region-util.h"

namespace xoos::svc {

constexpr f64 kDefaultMaxAfForSampleTypeDetection = 0.10;
constexpr f64 kDefaultMaxFfpePhredForSampleTypeDetection = 42.5;
constexpr u64 kDefaultWholeGenomeLengthForSampleTypeDetection = 3'000'000'000ULL;

/**
 * @brief Final metrics used to classify FFPE versus cell-line sample type.
 * @details These values are produced after VCF observation accumulation. The Phred quality score is derived from
 * error_count, total_length, and median_dp and is compared against the configured FFPE threshold.
 */
struct SampleTypeMetrics {
  /// Sum of low-frequency SNV alternate allele depths used as the error count.
  u64 error_count{};
  /// Denominator length from BED intervals, VCF contig headers, or fallback genome length.
  u64 total_length{};
  /// Median DP among records that contributed at least one error observation.
  u32 median_dp{};
  /// Phred-scaled error-rate estimate used for sample-type classification.
  f64 phred{};
};

/**
 * @brief Mutable accumulator for sample-type observations during VCF scanning.
 * @details Each sequential or per-thread scan writes into this lightweight struct. Parallel scans merge
 * per-task accumulators before final metrics are computed.
 */
struct SampleTypeMetricAccumulator {
  /// Running sum of low-frequency SNV alternate allele depths.
  u64 error_count{};
  /// DP values for records that contributed to error_count; used to compute median DP.
  vec<u32> dps;
};

/**
 * @brief Sum contig lengths declared in a VCF header.
 * @details Used as the whole-VCF denominator for sample-type detection. Missing contigs contribute no
 * length; callers decide whether to fall back when the returned total is zero.
 * @param header VCF header containing optional ##contig length metadata.
 * @return Total length across all contigs with length metadata.
 */
u64 SumVcfHeaderContigLengths(const io::VcfHeaderPtr& header);

/**
 * @brief Sum interval lengths in a chromosome interval map.
 * @details Used when sample-type detection is limited to BED targets, so the phred denominator matches
 * the regions that can contribute observations.
 * @param intervals_by_chromosome BED-style intervals grouped by chromosome.
 * @return Sum of interval end minus start across all chromosomes.
 */
u64 SumChromIntervalLengths(const ChromIntervalsMap& intervals_by_chromosome);

/**
 * @brief Add one VCF record to sample-type metric accumulation.
 * @details Reads DP and AD for @p sample_idx, counts SNV alternate allele depths below @p max_af as
 * low-frequency error observations, and stores the record DP when it contributes at least one observation.
 * Non-SNV records, missing fields, zero DP, and alleles above the AF threshold are ignored.
 * @param record VCF record to evaluate.
 * @param sample_idx Zero-based sample index used for FORMAT/DP and FORMAT/AD extraction.
 * @param max_af Maximum allele frequency for observations counted as errors.
 * @param accumulator Accumulator updated in place.
 */
void TallySampleTypeRecord(const io::VcfRecordPtr& record,
                           s32 sample_idx,
                           f64 max_af,
                           SampleTypeMetricAccumulator& accumulator);

/**
 * @brief Add one VCF record to sample-type metrics using an already extracted sample DP.
 * @details This avoids a second FORMAT/DP lookup when the caller has already read DP values for other
 * pre-scan work. FORMAT/AD is still read only after @p dp is known to be positive and the record is SNV-like.
 * @param record VCF record to evaluate.
 * @param sample_idx Zero-based sample index used for FORMAT/AD extraction.
 * @param dp Non-negative sample DP value already extracted from FORMAT/DP.
 * @param max_af Maximum allele frequency for observations counted as errors.
 * @param accumulator Accumulator updated in place.
 */
void TallySampleTypeRecordWithDp(
    const io::VcfRecordPtr& record, s32 sample_idx, u32 dp, f64 max_af, SampleTypeMetricAccumulator& accumulator);

/**
 * @brief Merge a per-task sample-type accumulator into a combined accumulator.
 * @details Used by parallel VCF scans after each worker has processed an independent chromosome region.
 * @param combined Destination accumulator updated in place.
 * @param current Source accumulator passed by value so callers can move per-task state into this helper.
 */
void MergeSampleTypeAccumulators(SampleTypeMetricAccumulator& combined, SampleTypeMetricAccumulator current);

/**
 * @brief Convert accumulated sample-type observations into final metrics.
 * @details Computes median DP from contributing records and calculates phred as
 * -10 * log10(error_count / (median_dp * total_length)). This function centralizes validation so both
 * standalone and filter_variants pre-scan paths fail consistently when input evidence is insufficient.
 * @param accumulator Error counts and DP values collected from VCF records.
 * @param total_length Denominator length from BED intervals, VCF header contigs, or fallback whole-genome length.
 * @param vcf_file VCF path used only for error messages.
 * @param max_af Maximum AF threshold used during accumulation, included in error messages.
 * @return Final metrics used by ClassifySampleType().
 * @throws xoos::error::Error if @p total_length is zero, no observations were collected, or median DP is zero.
 */
SampleTypeMetrics FinalizeSampleTypeMetrics(const SampleTypeMetricAccumulator& accumulator,
                                            u64 total_length,
                                            const fs::path& vcf_file,
                                            f64 max_af);

}  // namespace xoos::svc
