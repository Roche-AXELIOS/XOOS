#include "cascade-heterozygosity.h"

#include <algorithm>

#include <xoos/log/logging.h>

#include "cascade-contig-span.h"

namespace xoos::sex_predict {

vec<ScanWindow> SampleWindows(const s32 tid, const s64 start, const s64 end, const u32 count, const s64 size) {
  vec<ScanWindow> windows;
  const s64 span = end - start;
  if (span <= 0 || count == 0) {
    return windows;
  }
  if (span <= static_cast<s64>(count) * size) {
    windows.push_back({.tid = tid, .start = start, .end = end});
    return windows;
  }
  const s64 stride = span / count;
  windows.reserve(count);
  for (u32 i = 0; i < count; ++i) {
    const s64 window_start = start + (static_cast<s64>(i) * stride);
    windows.push_back({.tid = tid, .start = window_start, .end = std::min(window_start + size, end)});
  }
  return windows;
}

std::variant<SexMeasurements, SexPrediction> MeasureHeterozygosity(const fs::path& bam_file_path,
                                                                   const DepthProfile& depth_profile,
                                                                   const RegionInfo& region_x,
                                                                   SexMeasurements measurements,
                                                                   const u8 min_mapq,
                                                                   const u32 num_threads,
                                                                   const Sex ratio_sex) {
  using enum SexEvidence;
  using enum SexUncertaintyReason;

  const auto x_span = TrimmedSpan(depth_profile.contigs[depth_profile.x_index].len);
  if (!x_span.has_value()) {
    return SexPrediction{.sex = ratio_sex,
                         .evidence = kChrXYRatioFallback,
                         .confidence = SexConfidence::kLow,
                         .reason = kInsufficientHetSites,
                         .measurements = measurements};
  }
  const vec<ScanWindow> x_windows =
      SampleWindows(region_x.tid, x_span->first, x_span->second, kNumHetWindowsPerGroup, kHetWindowSize);

  // Spread the same total window budget across the sampled autosomes so both
  // groups are measured with comparable statistical weight.
  const u32 windows_per_autosome =
      std::max(1U, kNumHetWindowsPerGroup / static_cast<u32>(depth_profile.autosomes.size()));
  vec<ScanWindow> autosome_windows;
  for (const AutosomeContig& autosome : depth_profile.autosomes) {
    const auto span = TrimmedSpan(autosome.len);
    if (!span.has_value()) {
      continue;
    }
    const vec<ScanWindow> sampled =
        SampleWindows(autosome.tid, span->first, span->second, windows_per_autosome, kHetWindowSize);
    autosome_windows.insert(autosome_windows.end(), sampled.begin(), sampled.end());
  }

  const HetCounts chr_x_het =
      CountHeterozygousSites(bam_file_path, x_windows, min_mapq, kDefaultMinBaseQuality, num_threads);
  const HetCounts autosome_het =
      CountHeterozygousSites(bam_file_path, autosome_windows, min_mapq, kDefaultMinBaseQuality, num_threads);

  measurements.chr_x_het_rate = chr_x_het.Rate();
  measurements.autosome_het_rate = autosome_het.Rate();
  measurements.chr_x_evaluated_sites = chr_x_het.evaluated_sites;
  measurements.autosome_evaluated_sites = autosome_het.evaluated_sites;
  Logging::Debug(
      "sex_predict: chrX het {}/{} sites (rate={:.6f}, median depth={}), autosome het {}/{} sites (rate={:.6f}, "
      "median depth={})",
      chr_x_het.het_sites,
      chr_x_het.evaluated_sites,
      chr_x_het.Rate(),
      chr_x_het.median_site_depth,
      autosome_het.het_sites,
      autosome_het.evaluated_sites,
      autosome_het.Rate(),
      autosome_het.median_site_depth);

  if (chr_x_het.evaluated_sites < kMinEvaluatedHetSites || autosome_het.evaluated_sites < kMinEvaluatedHetSites ||
      autosome_het.Rate() <= 0.0) {
    // Not enough evaluable sites to trust a rate -- typically low coverage.
    return SexPrediction{.sex = ratio_sex,
                         .evidence = kChrXYRatioFallback,
                         .confidence = SexConfidence::kLow,
                         .reason = kInsufficientHetSites,
                         .measurements = measurements};
  }

  measurements.het_ratio = chr_x_het.Rate() / autosome_het.Rate();
  Logging::Debug("sex_predict: het_ratio={:.4f} (chrX depth copies={:.3f})",
                 measurements.het_ratio,
                 measurements.chr_x_depth_copies);
  return measurements;
}

}  // namespace xoos::sex_predict
