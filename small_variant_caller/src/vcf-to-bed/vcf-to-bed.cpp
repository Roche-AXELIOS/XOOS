#include "vcf-to-bed.h"

#include <algorithm>
#include <string>

#include <taskflow/algorithm/for_each.hpp>

#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/log/logging.h>
#include <xoos/types/vec.h>

#include "util/log-util.h"

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
  std::string prev_chrom, chrom_seq;
  vec<Interval> intervals{};
  vec<Interval>::iterator intervals_itr{};
  vec<Interval> var_intervals{};
  while (const auto& record = vcf_reader.GetNextRecord(BCF_UN_STR)) {
    const std::string& chrom = record->Chromosome();
    if (record->Position() < 0) {
      WarnAsErrorIfSet("Found VCF record with variant position < 0");
      continue;
    }
    const auto pos = static_cast<u64>(record->Position());  // VCF coordinates are 1-based, but this value is 0-based

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
    if (has_target_intervals) {
      if (intervals.empty() || intervals_itr == intervals.end()) {
        continue;
      }
      // keep advancing the iterator until the current interval's end is either at or
      // greater than this variant's start position
      while (intervals_itr != intervals.end() && intervals_itr->end < pos) {
        ++intervals_itr;
      }
      if (!IntervalOverlap(*intervals_itr, pos, pos + static_cast<s64>(ref.length()))) {
        continue;  // variant and interval do not overlap
      }
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
  int chrom_index;
  u64 start;
  u64 end;
};

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

  vec<VcfToBedParallelRegion> regions;
  bool use_single_thread;
  {
    // Extract contig indexes and lengths from VCF
    io::VcfReader vcf_reader(vcf_path);
    use_single_thread = threads <= 1 || !vcf_reader.HasIndex();
    if (!use_single_thread) {
      const auto contig_indexes = vcf_reader.GetContigIndexes();
      const auto contig_lengths = vcf_reader.GetHeader()->GetContigLengths();
      if (target_regions.empty()) {
        // one parallel task for each chromosome
        for (const auto& [chrom, chrom_length] : contig_lengths) {
          if (contig_indexes.contains(chrom) && chrom_length > 0) {
            const int cid = contig_indexes.at(chrom);
            regions.emplace_back(chrom, cid, 0, chrom_length);
          }
        }
      } else {
        // one parallel task for each target chromosome
        for (const auto& [chrom, intervals] : target_regions) {
          if (contig_indexes.contains(chrom) && !intervals.empty()) {
            const int cid = contig_indexes.at(chrom);
            regions.emplace_back(chrom, cid, intervals.begin()->start, (intervals.end() - 1)->end);
          }
        }
      }
    }
  }

  if (use_single_thread) {
    // fall back to using single thread
    return ExtractVariantIntervalsSingleThreaded(vcf_path, target_regions, left_pad, right_pad, collapse_dist);
  }

  // set up parallel tasks
  std::mutex chrom_to_var_intervals_mutex;
  ChromIntervalsMap chrom_to_var_intervals;
  auto task = [&](const VcfToBedParallelRegion& region) {
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
          const auto pos =
              static_cast<u64>(record->Position());  // VCF coordinates are 1-based, but this value is 0-based
          const std::string& ref{record->Allele(0)};
          if (has_target_regions) {
            if (target_intervals.empty() || target_intervals_itr == target_intervals.end()) {
              continue;
            }
            // keep advancing the iterator until the current interval's end is either at or
            // greater than this variant's start position
            while (target_intervals_itr != target_intervals.end() && target_intervals_itr->end < pos) {
              ++target_intervals_itr;
            }
            if (!IntervalOverlap(*target_intervals_itr, pos, pos + static_cast<s64>(ref.length()))) {
              continue;  // variant and interval do not overlap
            }
          }

          const u64 start = pos >= left_pad ? pos - left_pad : 0;
          const u64 end = pos + ref.length() + right_pad;
          var_intervals.emplace_back(start, end);
        }
      }
    } catch (const std::exception& e) {
      throw error::Error("Error in vcf region '{} {} {}': '{}'", chrom, region.start, region.end, e.what());
    }

    if (!var_intervals.empty()) {
      // do not store the vector if it is empty
      std::ranges::sort(var_intervals);
      const auto out_var_intervals = MergeIntervals(var_intervals, collapse_dist);
      const std::scoped_lock lock{chrom_to_var_intervals_mutex};
      chrom_to_var_intervals[chrom] = out_var_intervals;
    }
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
  const ChromIntervalsMap out_map = ExtractVariantIntervalsParallelized(
      param.vcf_file, chrom_intervals, param.threads, param.left_pad, param.right_pad, param.collapse_dist);

  Logging::Info("Writing BED file...");
  WriteBed(out_map, param.output_bed_file);
  Logging::Info("Converted VCF to BED file.");
}

}  // namespace xoos::svc
