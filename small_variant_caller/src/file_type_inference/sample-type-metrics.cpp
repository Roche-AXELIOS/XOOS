#include "sample-type-metrics.h"

#include <cmath>
#include <iterator>
#include <numeric>
#include <ranges>
#include <utility>

#include <xoos/error/error.h>
#include <xoos/util/math.h>

#include "core/vcf-fields.h"

namespace xoos::svc {

namespace {
/**
 * @brief Coordinates for extracting one allele depth from flattened FORMAT/AD values.
 * @details FORMAT/AD is laid out as sample-major allele values. Carrying the sample index, allele count,
 * and allele index together makes the flattened index calculation explicit at the call site.
 */
struct AltAdIndex {
  /// Zero-based VCF sample index.
  s32 sample_idx{};
  /// Number of alleles in the record, including REF.
  s32 num_alleles{};
  /// Allele index within the record; alternate alleles start at 1.
  s32 allele_index{};
};

/**
 * @brief Extract non-negative DP for one VCF sample.
 * @details Missing FORMAT/DP values, out-of-range sample indexes, and negative values are normalized to
 * zero so callers can skip records that do not provide usable depth.
 * @param record VCF record containing FORMAT/DP values.
 * @param sample_idx Zero-based sample index to read.
 * @return DP value for the sample, or zero when unavailable or invalid.
 */
u32 GetSampleDp(const io::VcfRecordPtr& record, const s32 sample_idx) {
  const auto values = record->GetFormatFieldNoCheck<s32>(kFieldDp);
  if (sample_idx < 0 || !std::cmp_less(sample_idx, values.size())) {
    return 0;
  }
  const s32 value = values.at(static_cast<size_t>(sample_idx));
  return value > 0 ? static_cast<u32>(value) : 0;
}

/**
 * @brief Extract non-negative alternate allele depth for one sample allele.
 * @details FORMAT/AD is stored as a flattened sample-by-allele array. This helper converts sample and
 * allele indexes into the flattened AD index and normalizes invalid or negative values to zero.
 * @param ad_values Flattened FORMAT/AD values from the VCF record.
 * @param index Sample, allele-count, and allele-index metadata for the lookup.
 * @return Alternate allele depth, or zero when unavailable or invalid.
 */
u32 GetAltAd(const vec<s32>& ad_values, const AltAdIndex& index) {
  const auto alt_idx = (index.sample_idx * index.num_alleles) + index.allele_index;
  if (alt_idx < 0 || !std::cmp_less(alt_idx, ad_values.size())) {
    return 0;
  }
  const s32 value = ad_values.at(static_cast<size_t>(alt_idx));
  return value > 0 ? static_cast<u32>(value) : 0;
}

/**
 * @brief Count low-allele-frequency SNV observations in one VCF record.
 * @details Only single-nucleotide alternate alleles are considered. An alternate allele contributes its AD
 * count when AD is positive and at or below @p max_af multiplied by the sample DP. This approximates the error
 * evidence used by the sample-type phred calculation.
 * @param record VCF record to evaluate.
 * @param sample_idx Zero-based sample index used for AD extraction.
 * @param dp Sample DP for the record.
 * @param max_af Maximum allele frequency for an observation to be counted.
 * @return Sum of contributing alternate allele depths for the record.
 */
u64 CountLowAlleleFrequencySnvErrors(const io::VcfRecordPtr& record,
                                     const s32 sample_idx,
                                     const u32 dp,
                                     const f64 max_af) {
  const s32 num_alleles = record->NumAlleles();
  if (num_alleles <= 1 || record->Allele(0).size() != 1) {
    return 0;
  }

  const auto ad_values = record->GetFormatFieldNoCheck<s32>(kFieldAd);
  if (ad_values.empty()) {
    return 0;
  }

  u64 error_count = 0;
  const f64 max_alt_count = max_af * static_cast<f64>(dp);
  for (s32 allele_index = 1; allele_index < num_alleles; ++allele_index) {
    if (record->Allele(allele_index).size() != 1) {
      continue;
    }
    const u32 alt_ad = GetAltAd(ad_values,
                                {
                                    .sample_idx = sample_idx,
                                    .num_alleles = num_alleles,
                                    .allele_index = allele_index,
                                });
    if (alt_ad == 0) {
      continue;
    }
    if (static_cast<f64>(alt_ad) <= max_alt_count) {
      error_count += alt_ad;
    }
  }
  return error_count;
}
}  // namespace

u64 SumVcfHeaderContigLengths(const io::VcfHeaderPtr& header) {
  const auto contig_lengths = header->GetContigLengths();
  return std::accumulate(contig_lengths.begin(), contig_lengths.end(), 0ULL, [](const u64 total, const auto& item) {
    return total + item.second;
  });
}

u64 SumChromIntervalLengths(const ChromIntervalsMap& intervals_by_chromosome) {
  u64 total_length = 0;
  for (const auto& intervals : intervals_by_chromosome | std::views::values) {
    for (const auto& interval : intervals) {
      total_length += interval.end - interval.start;
    }
  }
  return total_length;
}

void TallySampleTypeRecord(const io::VcfRecordPtr& record,
                           const s32 sample_idx,
                           const f64 max_af,
                           SampleTypeMetricAccumulator& accumulator) {
  const u32 dp = GetSampleDp(record, sample_idx);
  TallySampleTypeRecordWithDp(record, sample_idx, dp, max_af, accumulator);
}

void TallySampleTypeRecordWithDp(const io::VcfRecordPtr& record,
                                 const s32 sample_idx,
                                 const u32 dp,
                                 const f64 max_af,
                                 SampleTypeMetricAccumulator& accumulator) {
  if (dp == 0) {
    return;
  }

  const u64 record_error_count = CountLowAlleleFrequencySnvErrors(record, sample_idx, dp, max_af);
  if (record_error_count == 0) {
    return;
  }

  accumulator.error_count += record_error_count;
  accumulator.dps.emplace_back(dp);
}

void MergeSampleTypeAccumulators(SampleTypeMetricAccumulator& combined, SampleTypeMetricAccumulator current) {
  combined.error_count += current.error_count;
  combined.dps.insert(
      combined.dps.end(), std::make_move_iterator(current.dps.begin()), std::make_move_iterator(current.dps.end()));
}

SampleTypeMetrics FinalizeSampleTypeMetrics(const SampleTypeMetricAccumulator& accumulator,
                                            const u64 total_length,
                                            const fs::path& vcf_file,
                                            const f64 max_af) {
  if (total_length == 0) {
    throw error::Error("Cannot infer sample type because target region length is zero");
  }
  if (accumulator.error_count == 0 || accumulator.dps.empty()) {
    throw error::Error(
        "Cannot infer sample type from '{}' because no SNVs with AD > 0 and AF <= {} were found", vcf_file, max_af);
  }

  auto dps = accumulator.dps;
  const u32 median_dp = math::Median(dps);
  if (median_dp == 0) {
    throw error::Error("Cannot infer sample type from '{}' because median DP is zero", vcf_file);
  }

  const f64 denominator = static_cast<f64>(median_dp) * static_cast<f64>(total_length);
  const f64 phred = -10.0 * std::log10(static_cast<f64>(accumulator.error_count) / denominator);
  return {.error_count = accumulator.error_count, .total_length = total_length, .median_dp = median_dp, .phred = phred};
}

}  // namespace xoos::svc
