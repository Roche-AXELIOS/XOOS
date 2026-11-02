#include "cascade-depth-profile.h"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <limits>
#include <locale>
#include <numeric>
#include <string_view>
#include <utility>

#include <htslib/sam.h>

#include <taskflow/taskflow.hpp>

#include <xoos/log/logging.h>

#include "cascade-contig-span.h"
#include "scan-chunk.h"

namespace xoos::sex_predict {

namespace {

// Maximum coefficient of variation of per-autosome windowed depth for the
// baseline to be trusted. Above this the sampled autosomes disagree (e.g. a
// copy-number-noisy tumor with arm-level CNAs across them), so the copy-state
// estimate would be unreliable and PredictSexDetailed() falls back to the
// chrX/chrY ratio. A diploid genome's windowed autosome depths are tightly
// clustered (CV well below this).
constexpr f64 kMaxAutosomeDepthCv = 0.30;

// Number of largest autosomes sampled to build the diploid depth baseline.
// Sampling several -- not just chr1 -- lets their depth dispersion reveal a
// copy-number-noisy (e.g. aneuploid tumor) genome, and avoids trusting a
// single chromosome that may itself be aneuploid.
constexpr u32 kDefaultNumAutosomesToSample = 5;

// A coefficient of variation needs at least two observations to say anything
// about consistency. This also permits chromosome-chunked BAMs when two or
// more autosomes carry coverage.
constexpr size_t kMinCoveredAutosomes = 2;

// A recognized autosome contig: "chr1".."chr22" or "1".."22". Anything with a
// suffix (decoys like "chr1_KI270706v1_random"), the allosomes, chrM, and
// unplaced contigs are rejected.
bool IsAutosomeName(std::string name) {
  std::ranges::transform(name, name.begin(), [](const char c) { return std::tolower(c, std::locale::classic()); });
  std::string_view digits = name;
  constexpr size_t kChrPrefixLength = 3;
  if (digits.starts_with("chr")) {
    digits.remove_prefix(kChrPrefixLength);
  }
  // "1".."22" -- at most two digits.
  constexpr size_t kMaxAutosomeDigits = 2;
  if (digits.empty() || digits.size() > kMaxAutosomeDigits) {
    return false;
  }
  s32 value = 0;
  for (const char c : digits) {
    // std::isdigit's locale overload returns bool, so compare as one.
    if (!std::isdigit(c, std::locale::classic())) {
      return false;
    }
    constexpr s32 kDecimalBase = 10;
    value = (value * kDecimalBase) + (c - '0');
  }
  constexpr s32 kMaxAutosome = 22;
  return value >= 1 && value <= kMaxAutosome;
}

// The `num_autosomes` largest autosomes in the header, longest first.
vec<AutosomeContig> LargestAutosomes(const sam_hdr_t* const header, const u32 num_autosomes) {
  vec<AutosomeContig> autosomes;
  for (s32 i = 0; i < header->n_targets; ++i) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic) — htslib C array access
    if (IsAutosomeName(header->target_name[i])) {
      autosomes.push_back({.tid = i, .len = TargetLen(header, i)});
    }
  }
  std::ranges::sort(autosomes, [](const AutosomeContig& lhs, const AutosomeContig& rhs) { return lhs.len > rhs.len; });
  if (autosomes.size() > num_autosomes) {
    autosomes.resize(num_autosomes);
  }
  return autosomes;
}

// Per-contig depth is the median over covered windows of this width, not
// total_reads/contig_length. The median ignores centromeres and assembly gaps
// (which otherwise dilute the denominator -- worst on chrY, where ~30 Mb of its
// 57 Mb is heterochromatin stored as N in hg38) and is robust to focal CNAs.
constexpr s64 kDepthWindowSize = 1'000'000;

// Median of the non-zero entries of `values`; 0.0 if there are none. Windows
// with no coverage are dropped rather than averaged in: they are usually
// reference-assembly gaps, not evidence of low copy number.
f64 MedianOfCovered(vec<f64> values) {
  std::erase_if(values, [](const f64 value) { return value <= 0.0; });
  if (values.empty()) {
    return 0.0;
  }
  std::ranges::sort(values);
  const size_t upper_middle = values.size() / 2;
  if (values.size() % 2 != 0) {  // NOSONAR cpp:S109 -- two is the stable modulo base for parity.
    return values[upper_middle];
  }
  // The literal is the divisor for the arithmetic mean of the two middle values.
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers,readability-magic-numbers)
  return (values[upper_middle - 1] + values[upper_middle]) / 2.0;
}

// Fixed-width chunks tiling `span` ([span.first, span.second)), tagged with
// `region_index`. Takes the span as a pair (rather than separate start/end
// params) so the two bounds can't be transposed at a call site.
vec<ScanChunk> ChunksForSpan(const s32 tid, const std::pair<s64, s64>& span, const s64 window, const u32 region_index) {
  vec<ScanChunk> chunks;
  s64 window_start = span.first;
  while (window_start < span.second) {
    const s64 window_end = std::min(window_start + window, span.second);
    chunks.push_back({.tid = tid, .start = window_start, .end = window_end, .region_index = region_index});
    window_start = window_end;
  }
  return chunks;
}

// Scan every contig in `contigs` in fixed windows and return each one's median
// window depth, in reads per base. All windows across all contigs run through a
// single taskflow pass sized to `num_threads`; each window writes to its own
// slot, so no synchronization is needed beyond the barrier.
vec<f64> MedianWindowDepths(const fs::path& bam_file_path,
                            const vec<ScanContig>& contigs,
                            const u8 min_mapq,
                            const u32 num_threads) {
  vec<ScanChunk> chunks;
  for (u32 i = 0; i < contigs.size(); ++i) {
    const auto span = TrimmedSpan(contigs[i].len);
    if (!span.has_value()) {
      continue;
    }
    vec<ScanChunk> contig_chunks = ChunksForSpan(contigs[i].tid, *span, kDepthWindowSize, i);
    chunks.insert(chunks.end(), contig_chunks.begin(), contig_chunks.end());
  }
  vec<f64> medians(contigs.size(), 0.0);
  if (chunks.empty()) {
    return medians;
  }

  const u32 thread_count = std::min(std::max(num_threads, 1U), static_cast<u32>(chunks.size()));
  vec<io::AlignmentReader> readers = io::OpenAlignmentReaders(bam_file_path, thread_count);
  vec<f64> window_depths(chunks.size(), 0.0);
  std::atomic<u64> unused{0};

  tf::Executor executor(thread_count);
  tf::Taskflow taskflow;
  for (size_t i = 0; i < chunks.size(); ++i) {
    taskflow.emplace([&readers, &executor, &chunks, &window_depths, i, min_mapq, &unused]() {
      const auto& reader = readers[ToUnsigned(executor.this_worker_id())];
      const ScanChunk& chunk = chunks[i];
      const u64 count =
          CountFilteredReads(reader.bam_file.get(), reader.idx.get(), chunk, min_mapq, /*min_reads=*/0, unused);
      window_depths[i] = static_cast<f64>(count) / static_cast<f64>(chunk.end - chunk.start);
    });
  }
  executor.run(taskflow).get();

  vec<vec<f64>> per_contig(contigs.size());
  for (size_t i = 0; i < chunks.size(); ++i) {
    per_contig[chunks[i].region_index].push_back(window_depths[i]);
  }
  for (size_t i = 0; i < contigs.size(); ++i) {
    medians[i] = MedianOfCovered(std::move(per_contig[i]));
  }
  return medians;
}

// Coefficient of variation of the non-zero entries of `values`; infinity when
// there is nothing to measure. Used as the autosome dispersion self-check.
f64 CoefficientOfVariation(vec<f64> values) {
  std::erase_if(values, [](const f64 value) { return value <= 0.0; });
  if (values.empty()) {
    return std::numeric_limits<f64>::infinity();
  }
  const f64 mean = std::accumulate(values.begin(), values.end(), 0.0) / static_cast<f64>(values.size());
  if (mean <= 0.0) {
    return std::numeric_limits<f64>::infinity();
  }
  f64 sum_sq = 0.0;
  for (const f64 value : values) {
    sum_sq += (value - mean) * (value - mean);
  }
  return std::sqrt(sum_sq / static_cast<f64>(values.size())) / mean;
}

// Copies carried by the autosome baseline, which is diploid by construction.
// Scaling a depth ratio by this converts it to absolute copy number.
constexpr f64 kDiploidCopyNumber = 2.0;

// chrX and chrY: always appended to the sampled autosomes as a pair.
constexpr size_t kNumAllosomes = 2;

// Continuous copy number of a chromosome relative to the diploid baseline;
// 0.0 when there is no baseline to compare against.
f64 CopiesFromDepth(const f64 chrom_depth, const f64 diploid_depth) {
  if (diploid_depth <= 0.0) {
    return 0.0;
  }
  return kDiploidCopyNumber * chrom_depth / diploid_depth;
}

}  // namespace

std::variant<DepthProfile, SexPrediction> BuildDepthProfile(const fs::path& bam_file_path,
                                                            const io::AlignmentReader& reader,
                                                            const RegionInfo& region_x,
                                                            const RegionInfo& region_y,
                                                            const u8 min_mapq,
                                                            const u32 num_threads,
                                                            const Sex ratio_sex) {
  using enum SexEvidence;
  using enum SexUncertaintyReason;

  const vec<AutosomeContig> autosomes = LargestAutosomes(reader.hdr.get(), kDefaultNumAutosomesToSample);
  if (autosomes.empty()) {
    return SexPrediction{.sex = ratio_sex,
                         .evidence = kChrXYRatioFallback,
                         .confidence = SexConfidence::kLow,
                         .reason = kNoAutosomeBaseline,
                         .measurements = std::nullopt};
  }

  vec<ScanContig> contigs;
  contigs.reserve(autosomes.size() + kNumAllosomes);
  for (const AutosomeContig& autosome : autosomes) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic) — htslib C array access
    contigs.push_back({.tid = autosome.tid, .len = autosome.len, .name = reader.hdr->target_name[autosome.tid]});
  }
  const size_t x_index = contigs.size();
  contigs.push_back({.tid = region_x.tid, .len = TargetLen(reader.hdr.get(), region_x.tid), .name = region_x.name});
  const size_t y_index = contigs.size();
  contigs.push_back({.tid = region_y.tid, .len = TargetLen(reader.hdr.get(), region_y.tid), .name = region_y.name});

  const vec<f64> depths = MedianWindowDepths(bam_file_path, contigs, min_mapq, num_threads);
  const vec<f64> autosome_depths(depths.begin(), depths.begin() + static_cast<s64>(autosomes.size()));
  const auto covered_autosome_count =
      std::ranges::count_if(autosome_depths, [](const f64 depth) { return depth > 0.0; });
  if (std::cmp_less(covered_autosome_count, kMinCoveredAutosomes)) {
    return SexPrediction{.sex = ratio_sex,
                         .evidence = kChrXYRatioFallback,
                         .confidence = SexConfidence::kLow,
                         .reason = kNoAutosomeBaseline,
                         .measurements = std::nullopt};
  }
  const f64 autosome_median = MedianOfCovered(autosome_depths);
  const f64 autosome_cv = CoefficientOfVariation(autosome_depths);

  for (size_t i = 0; i < contigs.size(); ++i) {
    const f64 copies = CopiesFromDepth(depths[i], autosome_median);
    Logging::Debug(
        "sex_predict: {} median window depth={:.6f} reads/bp ({:.3f} copies)", contigs[i].name, depths[i], copies);
  }
  Logging::Debug("sex_predict: autosome baseline median={:.6f} reads/bp, cv={:.4f}", autosome_median, autosome_cv);

  SexMeasurements measurements;
  measurements.autosome_depth_cv = autosome_cv;
  if (autosome_median <= 0.0) {
    return SexPrediction{.sex = ratio_sex,
                         .evidence = kChrXYRatioFallback,
                         .confidence = SexConfidence::kLow,
                         .reason = kNoAutosomeBaseline,
                         .measurements = std::nullopt};
  }
  measurements.chr_x_depth_copies = CopiesFromDepth(depths[x_index], autosome_median);
  measurements.chr_y_depth_copies = CopiesFromDepth(depths[y_index], autosome_median);

  if (autosome_cv > kMaxAutosomeDepthCv) {
    // The autosomes disagree too much to normalize against -- typically an
    // aneuploid tumor with arm-level changes across the sampled chromosomes.
    return SexPrediction{.sex = ratio_sex,
                         .evidence = kChrXYRatioFallback,
                         .confidence = SexConfidence::kLow,
                         .reason = kAutosomeDepthDispersed,
                         .measurements = measurements};
  }

  return DepthProfile{
      .autosomes = autosomes, .contigs = contigs, .x_index = x_index, .y_index = y_index, .measurements = measurements};
}

}  // namespace xoos::sex_predict
