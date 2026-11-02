#include "metrics/coverage-metrics/gc-bias-metrics.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include <csv.hpp>

#include <xoos/io/metadata-util.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/util/format-functions.h>

#include "metrics/metrics-names.h"

namespace xoos::alignment_metrics {

constexpr u8 kGcBiasNormalizedCoveragePrecision = 6;
constexpr u8 kGcBiasDropoutPrecision = 6;
constexpr u8 kGcBiasQuintileNormalizedCoveragePrecision = 6;

GcBiasMetrics& GcBiasMetrics::operator+=(const GcBiasMetrics& other) {
  for (size_t i = 0; i < windows_per_gc_pct.size(); ++i) {
    windows_per_gc_pct[i] += other.windows_per_gc_pct[i];
    reads_per_gc_pct[i] += other.reads_per_gc_pct[i];
    post_filter_reads_per_gc_pct[i] += other.post_filter_reads_per_gc_pct[i];
  }
  total_windows += other.total_windows;
  total_reads += other.total_reads;
  total_filtered_reads += other.total_filtered_reads;
  return *this;
}

void GcBiasMetrics::AddRead(size_t gc_bin_index, bool passed_read_filter) {
  if (gc_bin_index <= kMaxGcPercentage) {
    ++reads_per_gc_pct[gc_bin_index];
    ++total_reads;
    if (passed_read_filter) {
      ++post_filter_reads_per_gc_pct[gc_bin_index];
      ++total_filtered_reads;
    }
  }
}

f64 GcBiasMetrics::CalculateMeanReadsPerWindow(const bool post_filter) const {
  u64 reads = 0;
  u64 windows = 0;
  for (size_t i = 0; i < reads_per_gc_pct.size(); ++i) {
    reads += post_filter ? post_filter_reads_per_gc_pct[i] : reads_per_gc_pct[i];
    windows += windows_per_gc_pct[i];
  }
  if (windows > 0) {
    return static_cast<f64>(reads) / static_cast<f64>(windows);
  }
  return 0.0;
}

vec<f64> GcBiasMetrics::CalculateNormalizedCoverage(const bool post_filter) const {
  vec<f64> normalized_coverage_per_gc(kMaxGcPercentage + 1, 0.0);
  const f64 mean_reads_per_window = CalculateMeanReadsPerWindow(post_filter);
  // Check if the mean reads per window is effectively zero to avoid division by zero.
  if (std::abs(mean_reads_per_window) < std::numeric_limits<f64>::epsilon()) {
    return normalized_coverage_per_gc;
  }
  for (size_t i = 0; i < windows_per_gc_pct.size(); ++i) {
    if (windows_per_gc_pct[i] > 0) {
      normalized_coverage_per_gc[i] =
          (static_cast<f64>(post_filter ? post_filter_reads_per_gc_pct[i] : reads_per_gc_pct[i]) /
           windows_per_gc_pct[i]) /
          mean_reads_per_window;
    }
  }
  return normalized_coverage_per_gc;
}

vec<f64> GcBiasMetrics::CalculateErrorBar(const bool post_filter) const {
  vec<f64> error_bar_per_gc(kMaxGcPercentage + 1, 0.0);
  const f64 mean_reads_per_window = CalculateMeanReadsPerWindow(post_filter);
  for (size_t i = 0; i < windows_per_gc_pct.size(); ++i) {
    if (windows_per_gc_pct[i] > 0 && mean_reads_per_window > 0.0) {
      const u32 reads = post_filter ? post_filter_reads_per_gc_pct[i] : reads_per_gc_pct[i];
      error_bar_per_gc[i] = (std::sqrt(static_cast<f64>(reads)) / windows_per_gc_pct[i]) / mean_reads_per_window;
    }
  }
  return error_bar_per_gc;
}

vec<f64> GcBiasMetrics::CalculateBinDropout(const bool post_filter) const {
  vec<f64> dropout_per_gc(kMaxGcPercentage + 1, 0.0);
  for (size_t i = 0; i < windows_per_gc_pct.size(); ++i) {
    const u32 reads = post_filter ? post_filter_reads_per_gc_pct[i] : reads_per_gc_pct[i];
    const u64 reads_total = post_filter ? total_filtered_reads : total_reads;
    if (total_windows != 0 && reads_total != 0) {
      const f64 ref_percent = windows_per_gc_pct[i] / static_cast<f64>(total_windows);
      const f64 read_percent = reads / static_cast<f64>(reads_total);
      dropout_per_gc[i] = (ref_percent - read_percent) * 100.0;
    }
  }
  return dropout_per_gc;
}

std::array<f64, kNumQuintiles> GcBiasMetrics::CalculateQuintileNormalizedCoverage(const bool post_filter) const {
  std::array<f64, kNumQuintiles> quintile_normalized_coverage = {0.0, 0.0, 0.0, 0.0, 0.0};
  constexpr u32 kQuintileSize = 20;
  const f64 mean_reads_per_window = CalculateMeanReadsPerWindow(post_filter);
  // Check if the mean reads per window is effectively zero to avoid division by zero.
  if (std::abs(mean_reads_per_window) < std::numeric_limits<f64>::epsilon()) {
    return quintile_normalized_coverage;
  }
  vec<u32> quintile_reads(kNumQuintiles, 0);
  vec<u32> quintile_windows(kNumQuintiles, 0);
  for (size_t i = 0; i < windows_per_gc_pct.size(); ++i) {
    const size_t quintile_index = std::min(i / kQuintileSize, static_cast<size_t>(kNumQuintiles - 1));
    if (windows_per_gc_pct[i] > 0) {
      quintile_reads[quintile_index] += post_filter ? post_filter_reads_per_gc_pct[i] : reads_per_gc_pct[i];
      quintile_windows[quintile_index] += windows_per_gc_pct[i];
    }
  }
  for (size_t i = 0; i < quintile_normalized_coverage.size(); ++i) {
    if (quintile_windows[i] > 0) {
      quintile_normalized_coverage.at(i) =
          (static_cast<f64>(quintile_reads[i]) / quintile_windows[i]) / mean_reads_per_window;
    }
  }
  return quintile_normalized_coverage;
}

GcBiasMetrics::DropoutMetrics GcBiasMetrics::CalculateATGCDropout(const bool post_filter) const {
  DropoutMetrics at_gc_dropout;  // at_gc_dropout.at_dropout = AT dropout, at_gc_dropout.gc_dropout = GC dropout
  vec<f64> dropout_per_gc = CalculateBinDropout(post_filter);
  // sum the AT dropout values (bins 0 through 50, i.e. GC% <= 50)
  for (size_t i = 0; i <= dropout_per_gc.size() / 2; ++i) {  // NOSONAR cpp:S109
    at_gc_dropout.at_dropout += std::max(0.0, dropout_per_gc[i]);
  }
  // sum the GC dropout values (bins 50 through 100, i.e. GC% > 50)
  for (size_t i = (dropout_per_gc.size() / 2) + 1; i < dropout_per_gc.size(); ++i) {  // NOSONAR cpp:S109
    at_gc_dropout.gc_dropout += std::max(0.0, dropout_per_gc[i]);
  }
  return at_gc_dropout;
}

void GcBiasMetrics::WriteTsv(const fs::path& output_file, const io::TsvMetadata& tsv_metadata) const {
  std::ofstream ofstream(output_file);
  auto writer = csv::make_tsv_writer(ofstream).set_auto_flush(false);
  io::WriteTsvMetadata(ofstream, tsv_metadata);
  // Write header
  writer << vec<std::string>{kNameGcPercentage,
                             kNameWindowCount,
                             kNameReadCount,
                             kNamePostFilterReadCount,
                             kNameNormalizedCoverage,
                             kNameNormalizedPostFilterCoverage,
                             kNameErrorBar,
                             kNamePostFilterErrorBar};
  vec<f64> normalized_coverage = CalculateNormalizedCoverage(false);
  vec<f64> normalized_post_filter_coverage = CalculateNormalizedCoverage(true);
  vec<f64> error_bar = CalculateErrorBar(false);
  vec<f64> post_filter_error_bar = CalculateErrorBar(true);
  for (size_t i = 0; i < kMaxGcPercentage + 1; ++i) {
    writer << std::make_tuple(
        i,
        windows_per_gc_pct[i],
        reads_per_gc_pct[i],
        post_filter_reads_per_gc_pct[i],
        util::ToStringWithPrecision(normalized_coverage[i], kGcBiasNormalizedCoveragePrecision),
        util::ToStringWithPrecision(normalized_post_filter_coverage[i], kGcBiasNormalizedCoveragePrecision),
        util::ToStringWithPrecision(error_bar[i], kGcBiasNormalizedCoveragePrecision),
        util::ToStringWithPrecision(post_filter_error_bar[i], kGcBiasNormalizedCoveragePrecision));
  }
}

void GcBiasMetrics::WriteSummaryTsv(const fs::path& output_file, const io::TsvMetadata& tsv_metadata) const {
  std::ofstream ofstream(output_file);
  auto writer = csv::make_tsv_writer(ofstream).set_auto_flush(false);
  io::WriteTsvMetadata(ofstream, tsv_metadata);
  // Write header
  writer << vec<std::string>{kNameMetricName, kNameValue};
  // Write summary metrics (total_reads, total_windows, AT dropout, GC dropout, normalized coverage for each quintile)
  auto quintile_normalized_coverage = CalculateQuintileNormalizedCoverage(false);
  DropoutMetrics dropout_metrics = CalculateATGCDropout(false);
  writer << std::make_tuple(kNameTotalReads, total_reads);
  writer << std::make_tuple(kNameTotalWindows, total_windows);
  writer << std::make_tuple(kNameATDropout,
                            util::ToStringWithPrecision(dropout_metrics.at_dropout, kGcBiasDropoutPrecision));
  writer << std::make_tuple(kNameGCDropout,
                            util::ToStringWithPrecision(dropout_metrics.gc_dropout, kGcBiasDropoutPrecision));
  for (size_t i = 0; i < quintile_normalized_coverage.size(); ++i) {
    writer << std::make_tuple(
        kNamesQuintilesNormalizedCoverage.at(i),
        util::ToStringWithPrecision(quintile_normalized_coverage.at(i), kGcBiasQuintileNormalizedCoveragePrecision));
  }
  // Write post-filter summary metrics
  auto post_filter_quintile_normalized_coverage = CalculateQuintileNormalizedCoverage(true);
  DropoutMetrics post_filter_dropout_metrics = CalculateATGCDropout(true);
  writer << std::make_tuple(kNamePostFilterTotalReads, total_filtered_reads);
  writer << std::make_tuple(
      kNamePostFilterATDropout,
      util::ToStringWithPrecision(post_filter_dropout_metrics.at_dropout, kGcBiasDropoutPrecision));
  writer << std::make_tuple(
      kNamePostFilterGCDropout,
      util::ToStringWithPrecision(post_filter_dropout_metrics.gc_dropout, kGcBiasDropoutPrecision));
  for (size_t i = 0; i < post_filter_quintile_normalized_coverage.size(); ++i) {
    writer << std::make_tuple(kNamesPostFilterQuintilesNormalizedCoverage.at(i),
                              util::ToStringWithPrecision(post_filter_quintile_normalized_coverage.at(i),
                                                          kGcBiasQuintileNormalizedCoveragePrecision));
  }
}
}  // namespace xoos::alignment_metrics
