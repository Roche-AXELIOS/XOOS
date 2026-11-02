#include "allele-counter.h"

#include <algorithm>
#include <array>
#include <filesystem>
#include <functional>
#include <limits>
#include <span>
#include <utility>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <taskflow/taskflow.hpp>

#include <xoos/error/error.h>
#include <xoos/io/alignment-reader.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/io/htslib-util/read-util.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::sex_predict {

f64 HetCounts::Rate() const {
  if (evaluated_sites == 0) {
    return 0.0;
  }
  return static_cast<f64>(het_sites) / static_cast<f64>(evaluated_sites);
}

namespace {

// SAM spec sentinel for "mapping quality not available"; excluded explicitly
// since it is never `< min_mapq`. Mirrors sex-predictor.cpp's read filter.
constexpr u8 kMapQNotAvailable = 255;

// BAM uses 0xff for every base when base qualities are unavailable. It must
// not be interpreted as a quality score above the configured threshold.
constexpr u8 kBaseQualityNotAvailable = 255;

// Unsigned to satisfy clang-tidy's hicpp-signed-bitwise (core.flag is unsigned).
constexpr u32 kExcludeFlags = ToUnsigned(BAM_FUNMAP) | ToUnsigned(BAM_FSECONDARY) | ToUnsigned(BAM_FQCFAIL) |
                              ToUnsigned(BAM_FDUP) | ToUnsigned(BAM_FSUPPLEMENTARY);

// A/C/G/T. Indices into the per-site allele count array; named so the switch
// below reads as bases rather than as bare offsets.
constexpr size_t kNumBases = 4;
constexpr size_t kBaseA = 0;
constexpr size_t kBaseC = 1;
constexpr size_t kBaseG = 2;
constexpr size_t kBaseT = 3;

// Indices into the *sorted-descending* allele counts: the major allele and the
// candidate minor allele that decides heterozygosity.
constexpr size_t kMajorAllele = 0;
constexpr size_t kMinorAllele = 1;

/// State handed to the htslib pileup engine's read-fetch callback.
struct PileupSource {
  htsFile* bam_file = nullptr;
  hts_itr_t* itr = nullptr;
  u8 min_mapq = 0;
  s32 read_error = 0;
};

/// `bam_plp_init` callback: yield the next read passing the same flag/MAPQ
/// filter the depth scan uses.
///
/// This MUST return a raw `s32` (htslib's `bam_plp_auto_f` contract: >= 0 on
/// success, < 0 at end of input). Do NOT substitute `io::SamItrNext` -- it
/// returns `bool`, and `false == 0` reads as success to the pileup engine,
/// which then loops forever on the same record.
///
/// `data`/`read` cannot be const-qualified: the signature is fixed by
/// htslib's `bam_plp_auto_f` typedef, and `read` is genuinely written by
/// `sam_itr_next` below.
s32 NextFilteredRead(void* data, bam1_t* read) noexcept {  // NOSONAR cpp:S994 -- fixed by bam_plp_auto_f's typedef
  // Recovering typed context from htslib's void* callback userdata; no alternative in the C API.
  auto* const source = static_cast<PileupSource*>(data);  // NOSONAR cpp:S5356,S5357,M23_093
  while (true) {
    const s32 ret = sam_itr_next(source->bam_file, source->itr, read);
    if (ret < 0) {
      if (ret < -1) {
        source->read_error = ret;
      }
      return ret;
    }
    if ((read->core.flag & kExcludeFlags) != 0) {
      continue;
    }
    if (read->core.qual < source->min_mapq || read->core.qual == kMapQNotAvailable) {
      continue;
    }
    return ret;
  }
}

/// Per-thread accumulator. Each worker owns a reader (a single htsFile is not
/// safe for concurrent use), a depth histogram over every covered site, and the
/// depths of sites meeting the heterozygous-allele criterion.
///
/// Heterozygosity is decided during the scan, but the depth band cannot be --
/// it depends on the median across all windows, which is not known until the
/// scan finishes. Keeping a histogram of all sites plus the depths of the few
/// het candidates lets both be resolved afterwards in one pass, without holding
/// a per-site record for tens of millions of positions.
struct HetScanWorker {
  io::AlignmentReader reader;
  vec<u64> depth_histogram = vec<u64>(kDepthHistogramSize, 0);
  vec<u32> het_candidate_depths;
};

/// Tally A/C/G/T at one pileup column, skipping deletions, reference skips, and
/// bases below `min_base_quality`.
std::array<u32, kNumBases> CountBasesAtColumn(const std::span<const bam_pileup1_t> column, const u8 min_base_quality) {
  std::array<u32, kNumBases> counts{};
  for (const bam_pileup1_t& entry : column) {
    // is_del/is_refskip are single-bit bitfields, not bool.
    if (entry.is_del != 0U || entry.is_refskip != 0U) {
      continue;
    }
    const u32 qpos = ToUnsigned(entry.qpos);
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic) — htslib C array access
    const u8 base_quality = bam_get_qual(entry.b)[qpos];
    if (base_quality == kBaseQualityNotAvailable || base_quality < min_base_quality) {
      continue;
    }
    switch (io::GetBase(bam_get_seq(entry.b), qpos)) {
      case 'A':
        ++counts[kBaseA];
        break;
      case 'C':
        ++counts[kBaseC];
        break;
      case 'G':
        ++counts[kBaseG];
        break;
      case 'T':
        ++counts[kBaseT];
        break;
      default:
        // 'N' and any other IUPAC ambiguity code carry no allele information.
        break;
    }
  }
  return counts;
}

/// Whether the top two alleles at a site look heterozygous. `total` is the sum
/// of `counts`, passed in to avoid recomputing it.
bool IsHeterozygous(std::array<u32, kNumBases> counts, const u32 total) {
  if (total == 0) {
    return false;
  }
  std::ranges::sort(counts, std::greater{});
  const u32 second = counts[kMinorAllele];
  if (counts[kMajorAllele] < kMinAlleleSupportingReads || second < kMinAlleleSupportingReads) {
    return false;
  }
  const f64 fraction = static_cast<f64>(second) / static_cast<f64>(total);
  return fraction >= kMinHetAlleleFraction && fraction <= kMaxHetAlleleFraction;
}

/// Pileup one window, recording every covered site's depth and the depths of
/// heterozygous candidates.
void ScanWindowForHeterozygosity(HetScanWorker& worker,
                                 const ScanWindow& window,
                                 const u8 min_mapq,
                                 const u8 min_base_quality) {
  const io::HtsItrPtr itr = io::SamItrQueryI(worker.reader.idx.get(), window.tid, window.start, window.end);
  PileupSource source{.bam_file = worker.reader.bam_file.get(), .itr = itr.get(), .min_mapq = min_mapq};
  // BamPlpPtr has a custom deleter (BamPlpDeleter); std::make_unique cannot
  // construct a unique_ptr with a non-default deleter.
  const io::BamPlpPtr pileup(bam_plp_init(NextFilteredRead, &source));  // NOSONAR cpp:M23_286,S5356
  // Default max depth would silently truncate deep columns, distorting both the
  // histogram and the allele fractions derived from it.
  bam_plp_set_maxcnt(pileup.get(), std::numeric_limits<s32>::max());

  s32 column_tid = 0;
  s32 column_pos = 0;
  s32 column_depth = 0;
  while (true) {
    const bam_pileup1_t* const column = bam_plp_auto(pileup.get(), &column_tid, &column_pos, &column_depth);
    if (column == nullptr) {
      break;
    }
    // The iterator yields every read overlapping the window, so columns outside
    // it are only partially covered and must not be counted.
    if (column_pos < window.start || column_pos >= window.end) {
      continue;
    }
    const std::array<u32, kNumBases> counts =
        CountBasesAtColumn({column, static_cast<size_t>(column_depth)}, min_base_quality);
    const u32 total = counts[kBaseA] + counts[kBaseC] + counts[kBaseG] + counts[kBaseT];
    if (total == 0) {
      continue;
    }
    const u32 bucket = std::min(total, kDepthHistogramSize - 1);
    ++worker.depth_histogram[bucket];
    if (IsHeterozygous(counts, total)) {
      worker.het_candidate_depths.push_back(total);
    }
  }
  if (source.read_error < -1) {
    throw error::Error("Failed to read BAM '{}' while measuring heterozygosity (htslib return code {})",
                       worker.reader.bam_file->fn,
                       source.read_error);
  }
}

/// Median of a depth histogram, ignoring sites below `kMinSiteDepth` so that
/// shallow edge coverage does not drag the reference depth down.
u32 MedianDepth(const vec<u64>& histogram) {
  u64 total = 0;
  for (u32 depth = kMinSiteDepth; depth < kDepthHistogramSize; ++depth) {
    total += histogram[depth];
  }
  if (total == 0) {
    return 0;
  }
  const u64 midpoint = total / 2;
  u64 cumulative = 0;
  for (u32 depth = kMinSiteDepth; depth < kDepthHistogramSize; ++depth) {
    cumulative += histogram[depth];
    if (cumulative > midpoint) {
      return depth;
    }
  }
  return kDepthHistogramSize - 1;
}

}  // namespace

HetCounts CountHeterozygousSites(const fs::path& bam_file_path,
                                 const vec<ScanWindow>& windows,
                                 const u8 min_mapq,
                                 const u8 min_base_quality,
                                 const u32 num_threads) {
  if (windows.empty()) {
    return HetCounts{};
  }

  const u32 thread_count = std::min(std::max(num_threads, 1U), static_cast<u32>(windows.size()));
  vec<io::AlignmentReader> readers = io::OpenAlignmentReaders(bam_file_path, thread_count);
  vec<HetScanWorker> workers;
  workers.reserve(thread_count);
  for (auto& reader : readers) {
    workers.push_back(HetScanWorker{
        .reader = std::move(reader), .depth_histogram = vec<u64>(kDepthHistogramSize, 0), .het_candidate_depths = {}});
  }

  tf::Executor executor(thread_count);
  tf::Taskflow taskflow;
  for (const auto& window : windows) {
    taskflow.emplace([&workers, &executor, window, min_mapq, min_base_quality]() {
      ScanWindowForHeterozygosity(workers[ToUnsigned(executor.this_worker_id())], window, min_mapq, min_base_quality);
    });
  }
  executor.run(taskflow).get();

  vec<u64> histogram(kDepthHistogramSize, 0);
  vec<u32> het_candidate_depths;
  for (const auto& worker : workers) {
    for (size_t i = 0; i < kDepthHistogramSize; ++i) {
      histogram[i] += worker.depth_histogram[i];
    }
    het_candidate_depths.insert(
        het_candidate_depths.end(), worker.het_candidate_depths.begin(), worker.het_candidate_depths.end());
  }

  const u32 median = MedianDepth(histogram);
  if (median == 0) {
    return HetCounts{};
  }

  // Resolve the relative depth band now that the median is known, then apply it
  // to both the denominator (all covered sites) and the numerator (candidates).
  const auto lower = std::max(static_cast<u32>(kMinRelativeSiteDepth * median), kMinSiteDepth);
  const auto upper = static_cast<u32>(kMaxRelativeSiteDepth * median);

  HetCounts result;
  result.median_site_depth = median;
  for (u32 depth = lower; depth <= std::min(upper, kDepthHistogramSize - 1); ++depth) {
    result.evaluated_sites += histogram[depth];
  }
  for (const u32 depth : het_candidate_depths) {
    if (depth >= lower && depth <= upper) {
      ++result.het_sites;
    }
  }
  return result;
}

}  // namespace xoos::sex_predict
