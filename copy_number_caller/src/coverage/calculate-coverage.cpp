#include "coverage/calculate-coverage.h"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <gsl/gsl>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <utility>

#include <htslib/sam.h>

#include <csv.hpp>

#include <taskflow/taskflow.hpp>

#include <xoos/error/error.h>
#include <xoos/io/alignment-reader.h>
#include <xoos/io/htslib-util/cigar-util.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/io/malloc-ptr.h>
#include <xoos/io/metadata-util.h>
#include <xoos/log/logging.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/util/file-functions.h>

#include "cli/copy-number-caller-options.h"
#include "coverage/calculate-coverage-options.h"
#include "coverage/calculate-coverage-taskflow-graph.h"
#include "coverage/coverage.h"
#include "intervals/augmented-intervals.h"
#include "io/column-names.h"
#include "io/copy-number-caller-default-filenames.h"
#include "util/parse-util.h"

namespace xoos::cnc {

namespace {

bool HasRescuedReadTag(const bam1_t* const record) {
  const auto* const tag_ptr = bam_aux_get(record, "YF");
  return tag_ptr != nullptr && bam_aux2i(tag_ptr) == 1;
}

// Outcome of applying the read filter to a single record.
struct FilterResult {
  // whether the read contributes to coverage
  bool keep;
  // whether a secondary read was kept only via its YF:i:1 rescue tag
  bool rescued;
};

// Apply the flag-based read filtering policy, including the pangenome
// secondary-read rescue: a read excluded solely because it is secondary is kept
// when it carries YF:i:1.
FilterResult EvaluateReadFilter(const bam1_t* const record, const u32 exclude_flags) {
  const u32 flags = record->core.flag;
  const bool excluded = (flags & exclude_flags) != 0U;
  if (!excluded) {
    return {.keep = true, .rescued = false};
  }
  const bool excludes_secondary = (exclude_flags & static_cast<u32>(BAM_FSECONDARY)) != 0U;
  const bool is_secondary = (flags & static_cast<u32>(BAM_FSECONDARY)) != 0U;
  if (excludes_secondary && is_secondary && HasRescuedReadTag(record)) {
    // Re-evaluate ignoring the secondary bit: keep (and mark rescued) only if
    // no other excluded flag is set.
    const u32 excl_no_sec = exclude_flags & ~static_cast<u32>(BAM_FSECONDARY);
    if ((flags & excl_no_sec) == 0U) {
      return {.keep = true, .rescued = true};
    }
  }
  return {.keep = false, .rescued = false};
}

void WriteMapqs(const CoverageRecords& cov, const io::CommandLineInfo& command_line_info, const fs::path& out_file) {
  std::ofstream ofs(out_file);
  if (!ofs.is_open()) {
    Logging::Error("Failed to open mapping quality output file: {}", out_file.string());
    throw std::runtime_error("Failed to open mapping quality output file");
  }
  auto writer = csv::make_tsv_writer(ofs);
  if (!command_line_info.version.empty()) {
    io::WriteTsvMetadata(ofs, command_line_info);
  }
  const vec<std::string> column_names{kColumnContigAsComment, kColumnStart, kColumnEnd, kColumnMeanMapq};
  writer << column_names;
  for (size_t i = 0; i < cov.region.size(); ++i) {
    const auto& [contig, start, end] = ParseRegionString(cov.region[i]);
    // Fixed six-decimal-place output matches the other per-interval tracks and preserves decimals for whole-number
    // means.
    writer << std::make_tuple(
        contig, start, end, FloatAsString(cov.mean_mapping_quality[static_cast<Eigen::Index>(i)]));
  }
}

// A parsed interval within a chunk, with its index into CoverageRecords vectors.
struct ChunkInterval {
  s32 tid;
  // 0-based half-open
  s64 start;
  s64 end;
  Eigen::Index eidx;
};

struct RegListResult {
  io::HtsRegListPtr reglist;
  u32 num_chroms;
};

// Build an hts_reglist_t covering all intervals in the chunk. `intervals` MUST be sorted by
// (tid, start); CalculateCoverageRegion establishes this locally so this routine does not depend
// on any upstream ordering. htslib takes ownership of the calloc'd arrays and frees them with
// free() via HtsRegListDeleter, so they are allocated with calloc into MallocPtr rather than a
// std::vector/new.
RegListResult BuildRegList(const vec<ChunkInterval>& intervals, const io::AlignmentReader& reader) {
  // Count distinct chromosomes; valid because intervals are tid-sorted (see contract above).
  u32 num_chroms = 0;
  for (size_t i = 0; i < intervals.size(); ++i) {
    if (i == 0 || intervals[i].tid != intervals[i - 1].tid) {
      ++num_chroms;
    }
  }

  // clang-format off
  io::HtsRegListPtr reglist{static_cast<hts_reglist_t*>(std::calloc(num_chroms, sizeof(hts_reglist_t))),  // NOLINT(hicpp-no-malloc,cppcoreguidelines-no-malloc,cppcoreguidelines-owning-memory) NOSONAR(S1231,S984,S5356,S5357,S5488,S103)
                            io::HtsRegListDeleter{static_cast<int>(num_chroms)}};  // NOLINT(cppcoreguidelines-owning-memory) NOSONAR(S103)
  // clang-format on
  const gsl::span<hts_reglist_t> reglist_span{reglist.get(), num_chroms};

  size_t iv_pos = 0;
  for (u32 c = 0; c < num_chroms; ++c) {
    const s32 tid = intervals[iv_pos].tid;

    size_t count = 0;
    while (iv_pos + count < intervals.size() && intervals[iv_pos + count].tid == tid) {
      ++count;
    }

    // Own the interval array until it is handed to the reglist, so an exception in between cannot leak it.
    // clang-format off
    io::MallocPtr<hts_pair_pos_t[]> pairs{static_cast<hts_pair_pos_t*>(std::calloc(count, sizeof(hts_pair_pos_t)))};  // NOLINT(hicpp-no-malloc,cppcoreguidelines-no-malloc,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,cppcoreguidelines-avoid-c-arrays,cppcoreguidelines-owning-memory) NOSONAR(S1231,S984,S5356,S5357,S5945,S5488,S103)
    // clang-format on

    const gsl::span<hts_pair_pos_t> pairs_span{pairs.get(), count};
    for (size_t j = 0; j < count; ++j) {
      pairs_span[j] = hts_pair_pos_t{.beg = intervals[iv_pos + j].start, .end = intervals[iv_pos + j].end};
    }

    reglist_span[c].tid = tid;
    reglist_span[c].reg = sam_hdr_tid2name(reader.hdr.get(), tid);
    // The reglist now owns the interval array.
    reglist_span[c].intervals = pairs.release();
    reglist_span[c].count = static_cast<u32>(count);
    iv_pos += count;
  }

  return {.reglist = std::move(reglist), .num_chroms = num_chroms};
}

// Compute reference bases of `record` overlapping [iv_start, iv_end) via CIGAR.
// If ignore_dn=false (default), D/N ops count as covered positions (matches samtools bedcov).
// If ignore_dn=true, only M/X/= ops are counted.
s32 CigarOverlapBases(const bam1_t* const record, const s64 iv_start, const s64 iv_end, const bool ignore_dn) {
  s32 overlap = 0;
  s64 ref_pos = record->core.pos;
  const gsl::span<const u32> cigar{bam_get_cigar(record), record->core.n_cigar};
  for (const u32 cigar_op : cigar) {
    const u32 op = bam_cigar_op(cigar_op);
    if (io::ConsumesReference(op)) {
      const auto len = static_cast<s32>(bam_cigar_oplen(cigar_op));
      const bool is_gap = (op == BAM_CDEL || op == BAM_CREF_SKIP);
      if (!ignore_dn || !is_gap) {
        const s64 ovlp = std::min(ref_pos + len, iv_end) - std::max(ref_pos, iv_start);
        if (ovlp > 0) {
          overlap += static_cast<s32>(ovlp);
        }
      }
      ref_pos += len;
    }
  }
  return overlap;
}

}  // namespace

void CalculateCoverageMain(const CopyNumberCallerOptions& options) {
  Logging::Info("Reading augmented intervals file");
  std::ifstream augmented_intervals_istream(options.augmented_intervals_fname.value());
  const AugmentedIntervals augmented_intervals(augmented_intervals_istream);

  const auto coverage_out = options.output_dir / kDefaultCoverageOutput;
  const auto mapping_qualities_out = options.output_dir / kDefaultMapQOutput;
  file::CheckOutputFilePathsForWriting({coverage_out, mapping_qualities_out}, options.overwrite);

  tf::Taskflow taskflow;
  CalculateCoverageTaskflowGraph calculate_coverage_tf_graph(options.normal_bam_fname,
                                                             augmented_intervals,
                                                             options.calculate_coverage_options.exclude_flags,
                                                             options.calculate_coverage_options.ignore_DN);
  taskflow.composed_of(calculate_coverage_tf_graph);
  tf::Executor executor(options.threads);
  executor.run(taskflow).get();
  Logging::Info("Writing coverages to {}", coverage_out.string());
  const CoverageRecords cov = calculate_coverage_tf_graph.GetResult();
  if (!cov.has_duplicates) {
    Logging::Info("The input BAM file does not have any reads marked as duplicate");
  }
  std::ofstream ofs(coverage_out);
  if (!ofs.is_open()) {
    throw error::Error("Failed to open coverage output file: {}", coverage_out.string());
  }
  WriteMapqs(cov, options.command_line_info, mapping_qualities_out);
  cov.Write(ofs, options.command_line_info);
}

namespace {

// Parse and validate every region in [start_idx, end_idx), returning the intervals sorted by BAM
// (tid, start). Sorting here means the reglist grouping and the streaming interval pointer in
// CalculateCoverageRegion are correct by construction, independent of the upstream FAI order.
vec<ChunkInterval> ParseChunkIntervals(const CoverageRecords& res,
                                       const io::AlignmentReader& reader,
                                       const s32 start_idx,
                                       const s32 end_idx) {
  vec<ChunkInterval> chunk_intervals;
  chunk_intervals.reserve(static_cast<size_t>(end_idx - start_idx));
  for (s32 i = start_idx; i < end_idx; ++i) {
    const auto [chr, start_raw, end_raw] = ParseRegionString(res.region[static_cast<u64>(i)]);
    const auto start = static_cast<s64>(start_raw);
    const auto end = static_cast<s64>(end_raw);
    if (end <= start) {
      throw std::runtime_error("Error: record end < start in regions file");
    }
    const s32 tid = io::SamHdrName2Tid(reader.hdr.get(), chr);
    if (tid < 0) {
      throw error::Error("Contig '{}' not found in BAM header", chr);
    }
    chunk_intervals.push_back({.tid = tid, .start = start, .end = end, .eidx = static_cast<Eigen::Index>(i)});
  }
  std::ranges::sort(chunk_intervals, [](const ChunkInterval& lhs, const ChunkInterval& rhs) {
    return std::tie(lhs.tid, lhs.start) < std::tie(rhs.tid, rhs.start);
  });
  return chunk_intervals;
}

// Accumulate a single read's coverage into every chunk interval it overlaps. `iv_ptr` is advanced
// past intervals that ended before this read; it must persist across calls within one chunk since
// both reads and intervals arrive in (tid, pos) order, making the scan monotonic.
void AccumulateReadCoverage(const bam1_t* const record,
                            const vec<ChunkInterval>& chunk_intervals,
                            const bool ignore_dn,
                            const s32 start_idx,
                            const s32 end_idx,
                            size_t& iv_ptr,
                            CoverageRecords& res) {
  const s32 read_tid = record->core.tid;
  const s64 read_start = record->core.pos;
  const s64 read_end = bam_endpos(record);
  const size_t num_intervals = chunk_intervals.size();

  // Advance iv_ptr past intervals that have already ended before this read.
  while (iv_ptr < num_intervals &&
         (chunk_intervals[iv_ptr].tid < read_tid ||
          (chunk_intervals[iv_ptr].tid == read_tid && chunk_intervals[iv_ptr].end <= read_start))) {
    ++iv_ptr;
  }

  for (size_t j = iv_ptr; j < num_intervals; ++j) {
    const auto& iv = chunk_intervals[j];
    if (iv.tid > read_tid || iv.start >= read_end) {
      break;
    }
    const s32 bases = CigarOverlapBases(record, iv.start, iv.end, ignore_dn);
    if (bases > 0) {
      // eidx is within this chunk's exclusive [start_idx, end_idx) range (see the top-of-function
      // note in CalculateCoverageRegion), so these accumulations are race-free across chunks.
      Expects(iv.eidx >= start_idx && iv.eidx < end_idx);
      res.total_coverage[iv.eidx] += static_cast<f64>(bases);
      res.count[iv.eidx] += 1.0;
      // Accumulate MAPQ sum; divided into a mean in the finalise pass below.
      res.mean_mapping_quality[iv.eidx] += static_cast<f64>(record->core.qual);
    }
  }
}

// Convert accumulated MAPQ sums into per-interval means; NaN where no read covered the interval.
void FinaliseMeanMapq(CoverageRecords& res, const s32 start_idx, const s32 end_idx) {
  for (s32 i = start_idx; i < end_idx; ++i) {
    const auto eidx = static_cast<Eigen::Index>(i);
    if (res.count[eidx] > 0.0) {
      res.mean_mapping_quality[eidx] /= res.count[eidx];
    } else {
      res.mean_mapping_quality[eidx] = std::numeric_limits<f64>::quiet_NaN();
    }
  }
}

}  // namespace

void CalculateCoverageRegion(const io::AlignmentReader& reader,
                             CoverageRecords& res,
                             const u32 exclude_flags,
                             const bool ignore_dn,
                             const s32 start_idx,
                             const s32 step) {
  // This chunk owns the half-open CoverageRecords index range [start_idx, end_idx): BuildRegList
  // queries only these intervals, so every eidx written below is unique to this task. Concurrent
  // chunks therefore write disjoint slices of the total_coverage/count/mean_mapping_quality
  // vectors without synchronisation; the only cross-chunk shared state is the two std::atomic<bool>
  // flags (has_duplicates, has_rescued_secondaries).
  const auto total_regions = static_cast<s32>(res.region.size());
  const s32 end_idx = std::min(start_idx + step, total_regions);

  if (end_idx <= start_idx) {
    return;
  }

  // Parse and validate all intervals up front so we fail fast before touching the BAM, matching
  // the old per-interval error behaviour. Returned sorted by BAM (tid, start).
  const vec<ChunkInterval> chunk_intervals = ParseChunkIntervals(res, reader, start_idx, end_idx);

  // One multi-region BAM query covers the entire chunk: reads are returned in
  // coordinate order, allowing a single sequential BGZF pass.
  auto [reglist, num_chroms] = BuildRegList(chunk_intervals, reader);
  const auto itr = io::SamItrRegion(reader.idx.get(), reader.hdr.get(), reglist.release(), num_chroms);

  // Streaming interval pointer — advances monotonically with reads since both
  // reads and intervals arrive in (tid, pos) order.
  size_t iv_ptr = 0;

  const io::Bam1Ptr record(bam_init1());
  while (io::SamItrMultiNext(reader.bam_file.get(), itr.get(), record.get()) >= 0) {
    if ((record->core.flag & static_cast<u32>(BAM_FDUP)) != 0U) {
      res.has_duplicates = true;
    }
    const auto [keep, is_rescued] = EvaluateReadFilter(record.get(), exclude_flags);
    if (!keep) {
      continue;
    }
    if (is_rescued) {
      res.has_rescued_secondaries = true;
    }
    AccumulateReadCoverage(record.get(), chunk_intervals, ignore_dn, start_idx, end_idx, iv_ptr, res);
  }

  FinaliseMeanMapq(res, start_idx, end_idx);
}

}  // namespace xoos::cnc
