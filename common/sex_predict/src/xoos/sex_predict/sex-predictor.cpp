#include "sex-predictor.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <filesystem>
#include <locale>
#include <optional>
#include <random>
#include <string>
#include <vector>

#include <taskflow/taskflow.hpp>

#include <xoos/io/alignment-reader.h>
#include <xoos/log/logging.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "scan-chunk.h"
#include "sex.h"

namespace xoos::sex_predict {

// Find the first contig in the BAM header matching one of `chromosomes` (case-insensitive).
std::optional<RegionInfo> GetChromosomeForStrings(const io::SamHdrPtr& bam_header,
                                                  const std::vector<std::string>& chromosomes) {
  if (chromosomes.empty()) {
    return std::nullopt;
  }
  for (s32 i = 0; i < bam_header->n_targets; ++i) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic) — htslib C array access
    std::string target_name = bam_header->target_name[i];
    std::ranges::transform(
        target_name, target_name.begin(), [](const char c) { return std::tolower(c, std::locale::classic()); });
    if (std::ranges::find(chromosomes, target_name) != chromosomes.end()) {
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic) — htslib C array access
      return RegionInfo{.tid = i, .name = bam_header->target_name[i]};
    }
  }
  return std::nullopt;
}

std::optional<RegionInfo> GetXChromosome(const io::SamHdrPtr& bam_header) {
  return GetChromosomeForStrings(bam_header, std::vector<std::string>{"chrx", "x"});
}

std::optional<RegionInfo> GetYChromosome(const io::SamHdrPtr& bam_header) {
  return GetChromosomeForStrings(bam_header, std::vector<std::string>{"chry", "y"});
}

namespace {

// Chunk size for the parallel scan; 1 Mb measured fastest for this module's
// chrX+chrY scan size.
constexpr s64 kScanChunkSize = 1'000'000;

// Split one chromosome into fixed-width scan chunks over 0-based, half-open
// [start, end) coordinates (the XOOS convention). The final chunk is
// truncated to chrom_len, not padded past it.
vec<ScanChunk> ChunksForChromosome(const s32 tid, const s64 chrom_len, const u32 region_index) {
  vec<ScanChunk> chunks;
  s64 window_start = 0;
  while (window_start < chrom_len) {
    const s64 window_end = std::min(window_start + kScanChunkSize, chrom_len);
    chunks.push_back({.tid = tid, .start = window_start, .end = window_end, .region_index = region_index});
    window_start = window_end;
  }
  return chunks;
}

// Region indices used by the chrX/chrY ratio scan (BuildChunks) and its
// ScanWorker routing below.
constexpr u32 kRegionChrX = 0;
constexpr u32 kRegionChrY = 1;

// Build the full chunk list: shuffle each chromosome's chunks independently
// (fixed seed, for reproducibility), then interleave the two shuffled lists
// one-for-one. Only matters for `min_reads > 0` -- shuffling avoids
// concentrating an early-stopped scan on a chromosome's PAR-heavy start;
// interleaving prevents chrX chunks alone from starving chrY. See the
// module README for the empirical rationale.
vec<ScanChunk> BuildChunks(const s32 x_tid, const s64 x_len, const s32 y_tid, const s64 y_len) {
  vec<ScanChunk> x_chunks = ChunksForChromosome(x_tid, x_len, kRegionChrX);
  vec<ScanChunk> y_chunks = ChunksForChromosome(y_tid, y_len, kRegionChrY);
  constexpr u32 kChunkShuffleSeed = 20'260'812U;
  // Fixed seed for reproducible chunk order; not security-sensitive.
  // NOLINTNEXTLINE(cert-msc32-c,cert-msc51-cpp)
  std::mt19937 rng(kChunkShuffleSeed);  // NOSONAR(cpp:S2245)
  std::ranges::shuffle(x_chunks, rng);
  std::ranges::shuffle(y_chunks, rng);

  vec<ScanChunk> chunks;
  chunks.reserve(x_chunks.size() + y_chunks.size());
  const size_t max_chunk_count = std::max(x_chunks.size(), y_chunks.size());
  for (size_t i = 0; i < max_chunk_count; ++i) {
    if (i < x_chunks.size()) {
      chunks.push_back(x_chunks[i]);
    }
    if (i < y_chunks.size()) {
      chunks.push_back(y_chunks[i]);
    }
  }
  return chunks;
}

// Per-thread scan state: each worker needs its own reader (a single
// htsFile isn't safe for concurrent use across threads) plus running totals.
struct ScanWorker {
  io::AlignmentReader reader;
  u64 x_count = 0;
  u64 y_count = 0;
};

}  // namespace

// Classify a chrX/chrY ratio against the two ratio gates, shared by
// PredictSex() and step 1 of PredictSexDetailed() (see sex-predictor-cascade.cpp).
Sex ClassifyByRatio(const f64 ratio, const f32 female_threshold, const f32 male_threshold_or_na) {
  using enum Sex;
  if (std::isnan(ratio)) {
    return kUnknown;
  }
  if (ratio > female_threshold) {
    return kFemale;
  }
  if (ratio > male_threshold_or_na) {
    return kUnknown;
  }
  return kMale;
}

std::optional<ChrXYReadCounts> ScanChrXYReadCounts(const fs::path& bam_file_path,
                                                   const u8 min_mapq,
                                                   const u64 min_reads,
                                                   const u32 num_threads) {
  io::AlignmentReader header_reader = io::OpenAlignmentReader(bam_file_path);
  const auto region_x = GetXChromosome(header_reader.hdr);
  const auto region_y = GetYChromosome(header_reader.hdr);
  if (!region_x.has_value() || !region_y.has_value()) {
    Logging::Warn("Could not find chrX and/or chrY contig in BAM header; cannot predict sex.");
    return std::nullopt;
  }

  const s64 chr_x_len = TargetLen(header_reader.hdr.get(), region_x->tid);
  const s64 chr_y_len = TargetLen(header_reader.hdr.get(), region_y->tid);
  const vec<ScanChunk> chunks = BuildChunks(region_x->tid, chr_x_len, region_y->tid, chr_y_len);

  // One reader per thread; worker 0 reuses the reader already opened above
  // instead of opening the file again.
  const u32 thread_count = std::max(1U, std::min(std::max(num_threads, 1U), static_cast<u32>(chunks.size())));
  vec<ScanWorker> workers;
  workers.reserve(thread_count);
  workers.push_back(ScanWorker{.reader = std::move(header_reader)});
  for (u32 i = 1; i < thread_count; ++i) {
    workers.push_back(ScanWorker{.reader = io::OpenAlignmentReader(bam_file_path)});
  }

  // Shared total of filtered reads; only used when min_reads > 0.
  std::atomic<u64> total_counted{0};

  tf::Executor executor(thread_count);
  tf::Taskflow taskflow;
  for (const auto& chunk : chunks) {
    taskflow.emplace([&workers, &executor, chunk, min_mapq, min_reads, &total_counted]() {
      // Skip chunks not yet started once another has already met the threshold.
      if (min_reads > 0 && total_counted.load(std::memory_order_relaxed) >= min_reads) {
        return;
      }
      auto& worker = workers[ToUnsigned(executor.this_worker_id())];
      const u64 count = CountFilteredReads(
          worker.reader.bam_file.get(), worker.reader.idx.get(), chunk, min_mapq, min_reads, total_counted);
      if (chunk.region_index == kRegionChrX) {
        worker.x_count += count;
      } else {
        worker.y_count += count;
      }
    });
  }
  executor.run(taskflow).get();

  ChrXYReadCounts totals;
  for (const auto& worker : workers) {
    totals.chr_x += worker.x_count;
    totals.chr_y += worker.y_count;
  }
  return totals;
}

Sex PredictSex(const fs::path& bam_file_path, const u8 min_mapq, const u64 min_reads, const u32 num_threads) {
  const auto counts = ScanChrXYReadCounts(bam_file_path, min_mapq, min_reads, num_threads);
  if (!counts.has_value()) {
    return Sex::kUnknown;
  }
  return ClassifyByRatio(
      static_cast<f64>(counts->chr_x) / static_cast<f64>(counts->chr_y), kMinRatioForSexDet, kMinRatioForSexNA);
}

Sex PredictSex(const fs::path& bam_file_path, const u8 min_mapq, const u64 min_reads) {
  return PredictSex(bam_file_path, min_mapq, min_reads, kDefaultNumThreads);
}

Sex PredictSex(const fs::path& bam_file_path, const u8 min_mapq) {
  return PredictSex(bam_file_path, min_mapq, kDefaultMinReads);
}

Sex PredictSex(const fs::path& bam_file_path) {
  return PredictSex(bam_file_path, kDefaultMinMapQ);
}

}  // namespace xoos::sex_predict
