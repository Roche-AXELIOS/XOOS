#include "vcf-to-bed.h"

#include <algorithm>
#include <mutex>
#include <string>

#include <taskflow/algorithm/for_each.hpp>

#include <xoos/error/error.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/log/logging.h>
#include <xoos/types/vec.h>

#include "util/file-util.h"
#include "util/log-util.h"
#include "util/vcf-util.h"

namespace xoos::svc {

/**
 * @brief Extract chromosome intervals based on variant positions in a VCF file using a single thread.
 * @param vcf_path VCF file path
 * @param target_intervals Map of chromosome to vector of target intervals
 * @param left_pad Padding to the left of variant start position
 * @param right_pad Padding to the right of tvariant end position
 * @param collapse_dist Distance to collapse neighboring intervals
 * @return Map of chromosome to vectors of variant intervals
 */
namespace {
// Forward declaration — defined below, after VcfToBedParallelRegion.
bool PassesTargetFilter(const vec<Interval>& target_intervals, vec<Interval>::iterator& itr, u64 pos, u64 ref_len);
}  // namespace

ChromIntervalsMap ExtractVariantIntervalsSingleThreaded(const fs::path& vcf_path,
                                                        const ChromIntervalsMap& target_intervals,
                                                        const u32 left_pad,
                                                        const u32 right_pad,
                                                        const u32 collapse_dist) {
  // Overview:
  // 1. Iterate through the VCF records in the file.
  // 2. For each variant, check if it overlaps with the target intervals (if provided).
  // 3. Create a new interval based on variant start position and reference allele length with the specified padding.
  // 4. Store the intervals in a map, grouped by chromosome.
  // 5. Merge overlapping or adjacent intervals within the specified collapse distance.
  // 6. Return the map of chromosome to vectors of variant intervals.

  const bool has_target_intervals{!target_intervals.empty()};
  ChromIntervalsMap chrom_to_var_intervals;

  const io::VcfReader vcf_reader(vcf_path);
  // temporary variables for the "previous" chromosome
  std::string prev_chrom;
  vec<Interval> intervals{};
  vec<Interval>::iterator intervals_itr{};
  vec<Interval> var_intervals{};
  while (const auto& record = vcf_reader.GetNextRecord(BCF_UN_STR)) {
    const std::string& chrom = record->Chromosome();
    if (record->Position() < 0) {
      WarnAsErrorIfSet("Found VCF record with variant position < 0");
      continue;
    }
    // VCF coordinates are 1-based, but this value is 0-based
    const auto pos = static_cast<u64>(record->Position());

    if (prev_chrom != chrom) {
      // this record has a different chromosome than the previous record
      // reset all temp variables for the new chromosome
      if (!var_intervals.empty()) {
        chrom_to_var_intervals[prev_chrom] = var_intervals;
        var_intervals = {};
      }
      prev_chrom = chrom;
      if (has_target_intervals) {
        if (target_intervals.contains(chrom)) {
          intervals = target_intervals.at(chrom);
          intervals_itr = intervals.begin();
        } else {
          intervals = {};
          intervals_itr = {};
        }
      }
    }

    const std::string& ref{record->Allele(0)};
    if (has_target_intervals && !PassesTargetFilter(intervals, intervals_itr, pos, ref.length())) {
      continue;
    }

    const u64 start = pos >= left_pad ? pos - left_pad : 0;
    const u64 end = pos + ref.length() + right_pad;
    var_intervals.emplace_back(start, end);
  }
  if (!var_intervals.empty()) {
    chrom_to_var_intervals[prev_chrom] = var_intervals;
  }

  ChromIntervalsMap output_chrom_intervals;
  for (auto& [_chrom, _intervals] : chrom_to_var_intervals) {
    // merge overlapping or adjacent intervals within the collapse distance threshold
    std::ranges::sort(_intervals);
    output_chrom_intervals[_chrom] = MergeIntervals(_intervals, collapse_dist);
  }

  return output_chrom_intervals;
}

struct VcfToBedParallelRegion {
  std::string chrom;
  s32 chrom_index;
  u64 start;
  u64 end;
};

/**
 * @brief Padding and collapse parameters for VCF-to-BED interval extraction.
 *
 * @var left_pad      Number of bases to extend each variant interval to the left.
 * @var right_pad     Number of bases to extend each variant interval to the right.
 * @var collapse_dist Maximum gap between adjacent intervals that triggers merging.
 */
struct VcfToBedParams {
  u64 left_pad;
  u64 right_pad;
  u64 collapse_dist;
};

namespace {
/**
 * @brief Advance the target interval iterator and check whether a variant overlaps.
 *
 * @details Advances @p itr past any intervals whose end is before @p pos, then checks
 * whether the variant at [@p pos, @p pos + @p ref_len) overlaps the current interval.
 *
 * @param target_intervals Sorted target intervals for the current chromosome.
 * @param itr              Iterator into @p target_intervals; advanced past consumed intervals.
 * @param pos              0-based variant start position.
 * @param ref_len          Length of the reference allele.
 * @return true if the variant overlaps a target interval and should be kept.
 */
bool PassesTargetFilter(const vec<Interval>& target_intervals,
                        vec<Interval>::iterator& itr,
                        const u64 pos,
                        const u64 ref_len) {
  if (target_intervals.empty() || itr == target_intervals.end()) {
    return false;
  }
  // advance past intervals that end before this variant
  while (itr != target_intervals.end() && itr->end < pos) {
    ++itr;
  }
  if (itr == target_intervals.end()) {
    return false;
  }
  return IntervalOverlap(*itr, pos, pos + ref_len);
}

bool RequiresSingleThreadedProcessing(const bool has_index, const fs::path& vcf_path, const u32 threads) {
  return threads <= 1 || !has_index || !VcfHasRecords(vcf_path);
}

void AppendWholeContigRegions(vec<VcfToBedParallelRegion>& regions,
                              const auto& contig_indexes,
                              const auto& contig_lengths) {
  for (const auto& [chrom, chrom_length] : contig_lengths) {
    if (!contig_indexes.contains(chrom) || chrom_length == 0) {
      continue;
    }
    regions.emplace_back(chrom, contig_indexes.at(chrom), 0, chrom_length);
  }
}

void AppendTargetRegions(vec<VcfToBedParallelRegion>& regions,
                         const auto& contig_indexes,
                         const ChromIntervalsMap& target_regions) {
  for (const auto& [chrom, intervals] : target_regions) {
    if (!contig_indexes.contains(chrom) || intervals.empty()) {
      continue;
    }
    regions.emplace_back(chrom, contig_indexes.at(chrom), intervals.begin()->start, (intervals.end() - 1)->end);
  }
}

vec<VcfToBedParallelRegion> BuildParallelRegions(const io::VcfReader& vcf_reader,
                                                 const ChromIntervalsMap& target_regions) {
  vec<VcfToBedParallelRegion> regions;
  const auto contig_indexes = vcf_reader.GetContigIndexes();
  const auto contig_lengths = vcf_reader.GetHeader()->GetContigLengths();
  if (target_regions.empty()) {
    AppendWholeContigRegions(regions, contig_indexes, contig_lengths);
    return regions;
  }
  AppendTargetRegions(regions, contig_indexes, target_regions);
  return regions;
}

/**
 * @brief Process a single VCF region, extracting padded variant intervals and merging them.
 *
 * @details Opens a VCF reader scoped to the genomic coordinates in @p region and iterates over
 * all records within that range. Each variant position is padded according to @p params to form
 * an interval. When @p target_regions is non-empty and contains the region's chromosome, only
 * variants that overlap a target interval are retained. The resulting intervals are sorted,
 * merged within @p params.collapse_dist, and appended to the shared @p chrom_to_var_intervals
 * map under a mutex lock. Empty results are not stored.
 *
 * This function is the per-region work unit called from @c ExtractVariantIntervalsParallelized
 * via a Taskflow parallel loop. Each invocation creates its own @c VcfReader, so concurrent
 * calls on different regions are safe as long as the VCF file is indexed.
 *
 * @param region        Genomic region (chromosome name, contig index, start, end) to query.
 * @param vcf_path      Path to an indexed VCF/BCF file.
 * @param target_regions Map of chromosome to sorted target intervals. If non-empty, only variants
 *                       overlapping a target interval are kept.
 * @param params        Padding and collapse distance parameters.
 * @param chrom_to_var_intervals_mutex Mutex guarding writes to @p chrom_to_var_intervals.
 * @param chrom_to_var_intervals       Shared output map; merged intervals for this region's
 *                                     chromosome are inserted under the mutex.
 *
 * @throws xoos::error::Error if VCF parsing fails for this region.
 */
void ProcessVcfRegion(const VcfToBedParallelRegion& region,
                      const fs::path& vcf_path,
                      const ChromIntervalsMap& target_regions,
                      const VcfToBedParams& params,
                      std::mutex& chrom_to_var_intervals_mutex,
                      ChromIntervalsMap& chrom_to_var_intervals) {
  vec<Interval> var_intervals{};
  const auto& chrom = region.chrom;
  try {
    io::VcfReader vcf_reader(vcf_path);
    if (vcf_reader.SetRegion(region.chrom_index, region.start, region.end)) {
      const bool has_target_regions = !target_regions.empty() && target_regions.contains(chrom);
      vec<Interval> target_intervals{};
      vec<Interval>::iterator target_intervals_itr{};
      if (has_target_regions) {
        target_intervals = target_regions.at(chrom);
        target_intervals_itr = target_intervals.begin();
      }

      while (const auto& record = vcf_reader.GetNextRecord(BCF_UN_STR)) {
        if (record->Position() < 0) {
          WarnAsErrorIfSet("Found VCF record with position < 0");
          continue;
        }
        // VCF coordinates are 1-based, but this value is 0-based
        const auto pos = static_cast<u64>(record->Position());
        const std::string& ref{record->Allele(0)};
        if (has_target_regions && !PassesTargetFilter(target_intervals, target_intervals_itr, pos, ref.length())) {
          continue;
        }

        const u64 start = pos >= params.left_pad ? pos - params.left_pad : 0;
        const u64 end = pos + ref.length() + params.right_pad;
        var_intervals.emplace_back(start, end);
      }
    }
  } catch (const std::exception& e) {
    throw error::Error("Error in vcf region '{} {} {}': '{}'", chrom, region.start, region.end, e.what());
  }

  if (!var_intervals.empty()) {
    // do not store the vector if it is empty
    std::ranges::sort(var_intervals);
    const auto out_var_intervals = MergeIntervals(var_intervals, params.collapse_dist);
    const std::scoped_lock lock{chrom_to_var_intervals_mutex};
    chrom_to_var_intervals[chrom] = out_var_intervals;
  }
}

}  // namespace

/**
 * @brief Extract chromosome intervals based on variant positions in a VCF file using one or more threads.
 * @param vcf_path VCF file path
 * @param target_regions Map of chromosome to vector of target intervals
 * @param threads Number of threads
 * @param left_pad Padding to the left of variant start position
 * @param right_pad Padding to the right of tvariant end position
 * @param collapse_dist Distance to collapse neighboring intervals
 * @return Map of chromosome to vectors of variant intervals
 */
ChromIntervalsMap ExtractVariantIntervalsParallelized(const fs::path& vcf_path,
                                                      const ChromIntervalsMap& target_regions,
                                                      const u32 threads,
                                                      const u32 left_pad,
                                                      const u32 right_pad,
                                                      const u32 collapse_dist) {
  // Overview:
  // 1. If 1 thread is requested or if the VCF file does not have an index, fall back to single-threaded processing.
  // 2. Read the VCF header and extract reference contig indexes and lengths.
  // 3. Set up 1 task per chromosome; if target intervals are available, limit to chromosomes in target intervals only.
  // 4. Set up Taskflow thread pool to process each chromosome in parallel.
  // 5. For each task, parse the VCF records at the corresponding chromosome.
  //   5.1. For each variant, check if it overlaps with the target intervals (if provided).
  //   5.2. Create a new interval based on the variant start position and reference allele length with the specified
  //        padding, and append it to a vector.
  //   5.3. Once all VCF records are parsed, store the vector of intervals in a map, grouped by chromosome.
  // 6. Merge overlapping or adjacent intervals within the specified collapse distance.
  // 7. Return the map of chromosome to vectors of variant intervals.

  // A VCF with no records has nothing to extract, and planning per-chromosome tasks from its index
  // would fail: an index over a record-less file covers no sequence, which `GetContigIndexes` rejects.
  // The single-threaded path reads the file directly and returns an empty map.
  io::VcfReader vcf_reader(vcf_path);
  const bool use_single_thread = RequiresSingleThreadedProcessing(vcf_reader.HasIndex(), vcf_path, threads);
  if (use_single_thread) {
    // fall back to using single thread
    return ExtractVariantIntervalsSingleThreaded(vcf_path, target_regions, left_pad, right_pad, collapse_dist);
  }
  const vec<VcfToBedParallelRegion> regions = BuildParallelRegions(vcf_reader, target_regions);

  // set up parallel tasks
  std::mutex chrom_to_var_intervals_mutex;
  ChromIntervalsMap chrom_to_var_intervals;
  const VcfToBedParams params{.left_pad = left_pad, .right_pad = right_pad, .collapse_dist = collapse_dist};
  auto task = [&vcf_path, &target_regions, &params, &chrom_to_var_intervals_mutex, &chrom_to_var_intervals](
                  const VcfToBedParallelRegion& region) {
    ProcessVcfRegion(region, vcf_path, target_regions, params, chrom_to_var_intervals_mutex, chrom_to_var_intervals);
  };

  // Execute parallel tasks
  tf::Executor executor(threads);
  tf::Taskflow flow;
  tf::Task map_task = flow.for_each(regions.begin(), regions.end(), task);
  map_task.name("Map VCF variant position extraction task to region");
  executor.run(flow).get();

  return chrom_to_var_intervals;
}

/**
 * @brief Convert a VCF file to a BED file.
 * @param param Parameters extracted from CLI
 */
void ConvertVcfToBed(const VcfToBedParam& param) {
  const ChromIntervalsMap& chrom_intervals =
      param.chrom_intervals_map.has_value() ? param.chrom_intervals_map.value() : ChromIntervalsMap{};
  Logging::Info("BED regions found: {}", !chrom_intervals.empty());

  Logging::Info("Parsing VCF file...");
  const ChromIntervalsMap out_map = ExtractVariantIntervalsParallelized(param.vcf_file,
                                                                        chrom_intervals,
                                                                        static_cast<u32>(param.threads),
                                                                        param.left_pad,
                                                                        param.right_pad,
                                                                        param.collapse_dist);

  if (out_map.empty()) {
    // Graceful empty output: an input VCF with no variants yields a BED with no regions. BED has no
    // header or container framing, so the valid empty representation of that file is zero bytes.
    WarnAsErrorIfSet("No variants found in VCF file '{}'. Writing empty outputs.", param.vcf_file);
  }

  CreateParentDirectoryIfNotExists(param.output_bed_file);
  Logging::Info("Writing BED file...");
  WriteBed(out_map, param.output_bed_file);
  Logging::Info("Converted VCF to BED file.");
}

}  // namespace xoos::svc
