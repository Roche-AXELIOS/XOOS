#pragma once

#include <filesystem>
#include <optional>
#include <thread>
#include <vector>

#include <xoos/io/alignment-reader.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>

#include "sex-predictor-cascade.h"  // IWYU pragma: keep
#include "sex.h"

namespace xoos::sex_predict {

// Default minimum MAPQ for a read to count towards chrX/chrY coverage --
// excludes low-MAPQ (likely multi-mapped) reads that could otherwise bias
// the ratio towards Male.
constexpr u8 kDefaultMinMapQ = 20;

// Default number of taskflow worker threads for the full scan: all available
// hardware threads, matching xoos::cli::AddThreadCountOption's own
// 0-means-hardware_concurrency() convention. Every scan already clamps its
// effective thread count to at least 1, so a platform where this reports 0
// (rare, but documented as possible) is handled downstream.
inline const u32 kDefaultNumThreads = std::thread::hardware_concurrency();

// Default `min_reads`: 0 disables early stopping (every read is counted).
// See PredictSex()'s `min_reads` documentation below.
constexpr u64 kDefaultMinReads = 0;

// Ratio thresholds for the whole-chromosome, MAPQ-filtered chrX/chrY ratio.
// Sourced from PureCN's getSexFromCoverage() defaults; see the README's
// "Threshold derivation" section for the citation and full derivation.
constexpr f32 kMinRatioForSexDet = 25.0F;
constexpr f32 kMinRatioForSexNA = 20.0F;

struct RegionInfo {
  s32 tid;
  std::string name;
};

std::optional<RegionInfo> GetXChromosome(const io::SamHdrPtr& bam_header);
std::optional<RegionInfo> GetYChromosome(const io::SamHdrPtr& bam_header);
std::optional<RegionInfo> GetChromosomeForStrings(const io::SamHdrPtr& bam_header,
                                                  const std::vector<std::string>& chromosomes);

/// Raw, post-filter chrX/chrY read totals from a scan (see PredictSex()'s
/// `min_reads` documentation for what "post-filter" and early stopping mean).
struct ChrXYReadCounts {
  u64 chr_x = 0;
  u64 chr_y = 0;
};

/// Scan a BAM's chrX/chrY and return the raw, post-filter read totals
/// PredictSex() classifies from. `std::nullopt` if chrX and/or chrY aren't in
/// the BAM header. Exposed separately so tests can verify `min_reads` early
/// stopping directly; see PredictSex() for parameter documentation.
std::optional<ChrXYReadCounts> ScanChrXYReadCounts(const fs::path& bam_file_path,
                                                   u8 min_mapq,
                                                   u64 min_reads,
                                                   u32 num_threads);

// Classify a chrX/chrY ratio against the two ratio gates, shared by
// PredictSex() and step 1 of PredictSexDetailed() (see sex-predictor-cascade.h).
Sex ClassifyByRatio(f64 ratio, f32 female_threshold, f32 male_threshold_or_na);

/// Predict sex from a BAM's chrX/chrY read coverage.
/// @param min_mapq Minimum MAPQ (inclusive) for a primary, mapped,
/// non-duplicate, non-QC-fail read to count.
/// @param min_reads Minimum total (chrX+chrY, post-filter) number of reads at
/// which the scan may stop early and classify from whatever's been counted
/// so far. 0 (the default) disables early stopping. Only 1,000,000 has been
/// verified safe against real data -- see the README's "min_reads: opt-in
/// early stopping" section before using a smaller value.
/// @note An exact, parallel scan of the BAM, not an O(1) index lookup --
/// do not call this in a tight loop.
Sex PredictSex(const fs::path& bam_file_path, u8 min_mapq, u64 min_reads, u32 num_threads);

/// Overload using `kDefaultNumThreads` worker threads.
Sex PredictSex(const fs::path& bam_file_path, u8 min_mapq, u64 min_reads);

/// Overload using `kDefaultMinReads` (0, no early stopping).
Sex PredictSex(const fs::path& bam_file_path, u8 min_mapq);

/// Overload using `kDefaultMinMapQ`.
Sex PredictSex(const fs::path& bam_file_path);

}  // namespace xoos::sex_predict
