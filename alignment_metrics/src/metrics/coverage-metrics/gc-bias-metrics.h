#pragma once

#include <array>

#include <xoos/io/metadata-util.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::alignment_metrics {

// Maximum GC percentage value (inclusive)
constexpr u32 kMaxGcPercentage = 100;
constexpr u32 kNumQuintiles = 5;

/**
 * @brief Container for GC bias metrics.
 *
 * Tracks counts of windows mapped to GC percentages by individual percents,
 * and the number of reads whose start positions map to windows with that GC percent value.
 */
class GcBiasMetrics {
 public:
  vec<u32> windows_per_gc_pct = vec<u32>(kMaxGcPercentage + 1, 0);
  vec<u32> reads_per_gc_pct = vec<u32>(kMaxGcPercentage + 1, 0);
  vec<u32> post_filter_reads_per_gc_pct = vec<u32>(kMaxGcPercentage + 1, 0);
  u64 total_windows = 0;
  u64 total_reads = 0;
  u64 total_filtered_reads = 0;

  struct DropoutMetrics {
    f64 at_dropout = 0.0;
    f64 gc_dropout = 0.0;
  };

  /**
   * @brief Calculate the mean number of reads per window across all GC percentage bins.
   * @param post_filter If true, use post-filter read counts; otherwise, use pre-filter read counts.
   * @return The mean number of reads per window.
   */
  f64 CalculateMeanReadsPerWindow(bool post_filter) const;
  /**
   * @brief Calculate the normalized coverage for each GC percentage bin.
   * @param post_filter If true, use post-filter read counts; otherwise, use pre-filter read counts.
   * @return A vector of normalized coverage values for each GC percentage bin.
   */
  vec<f64> CalculateNormalizedCoverage(bool post_filter) const;
  /**
   * @brief Calculate the error bar width for each GC percentage bin.
   * @param post_filter If true, use post-filter read counts; otherwise, use pre-filter read counts.
   * @return A vector of error bar widths for each GC percentage bin.
   */
  vec<f64> CalculateErrorBar(bool post_filter) const;
  /**
   * @brief Calculate the dropout for each GC percentage bin.
   * @param post_filter If true, use post-filter read counts; otherwise, use pre-filter read counts.
   * @return A vector of dropout values for each GC percentage bin.
   */
  vec<f64> CalculateBinDropout(bool post_filter) const;
  /**
   * @brief Calculate the normalized coverage for each quintile of GC percentage bins.
   * @param post_filter If true, use post-filter read counts; otherwise, use pre-filter read counts.
   * @return A vector of normalized coverage values for each quintile of GC percentage bins.
   */
  std::array<f64, kNumQuintiles> CalculateQuintileNormalizedCoverage(bool post_filter) const;
  /**
   * @brief Calculate the AT Dropout (bins 0 through 50) and GC Dropout (bins 50 through 100) across GC percentage bins.
   * @param post_filter If true, use post-filter read counts; otherwise, use pre-filter read counts.
   * @return A DropoutMetrics struct containing AT and GC dropout values.
   */
  DropoutMetrics CalculateATGCDropout(bool post_filter) const;
  /**
   * @brief Write the GC bias metrics to a TSV file.
   * @param output_file The path to the output TSV file.
   * @param tsv_metadata Metadata to include in the TSV header.
   */
  void WriteTsv(const fs::path& output_file, const io::TsvMetadata& tsv_metadata) const;
  /**
   * @brief Write the GC bias summary metrics to a TSV file.
   * @param output_file The path to the output TSV file.
   * @param tsv_metadata Metadata to include in the TSV header.
   */
  void WriteSummaryTsv(const fs::path& output_file, const io::TsvMetadata& tsv_metadata) const;
  /**
   * @brief Increment the read count for a specific GC percentage bin.
   * @param gc_bin_index The index of the GC percentage bin (0-100).
   * @param passed_read_filter Whether the read passed read-level filtering.
   */
  void AddRead(size_t gc_bin_index, bool passed_read_filter);

  GcBiasMetrics& operator+=(const GcBiasMetrics& other);
};

}  // namespace xoos::alignment_metrics
