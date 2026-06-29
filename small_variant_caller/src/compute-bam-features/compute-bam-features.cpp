#include "compute-bam-features.h"

#include <optional>
#include <string>

#include <htslib/hts.h>

#include <taskflow/taskflow.hpp>

#include <xoos/io/fasta-reader.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/metadata-util.h>

#include "compute-bam-features/alignment-reader.h"
#include "progress-meter.h"
#include "util/file-util.h"
#include "util/locked-tsv-writer.h"
#include "util/log-util.h"
#include "util/region-util.h"
#include "xoos/io/htslib-util/kstring.h"

namespace xoos::svc {

u32 ParallelComputeBamFeatures(const ComputeBamFeaturesCliParams& cli_params) {
  return ComputeBamFeatures().ParallelComputeBamFeatures(cli_params);
}

void ComputeBamFeatures::ExecuteBamFeatureExtraction(const vec<Region>& regions_list,
                                                     const ComputeBamFeaturesCliParams& cli_params,
                                                     Progress& progress) {
  // For `taskflow` parallelization, A single ComputeFeatures call is used for all workflows, with specific enum flags
  // being set by the config to denote steps specific to each workflow. This reduces code duplication and streamlines
  // the codebase.

  tf::Executor executor(cli_params.threads);

  // Defines a function that each `taskflow` task runs. Each task will be passed a region, For each region a call to
  // compute features is made and the resulting features are appended to the output features file.
  tf::Taskflow flow;
  for (const auto& region : regions_list) {
    auto task = [this, &cli_params, &executor, &progress, &region]() {
      const s32 thread_id = executor.this_worker_id();
      ExtractBamFeatures(thread_id, region, cli_params, progress);
    };
    flow.emplace(task);
  }

  Logging::Info("Processing '{}' regions...", progress.total_region_count);
  executor.run(flow).get();
  Logging::Info("Completed processing '{}' regions...", progress.total_region_count);

  if (regions_list.size() != progress.total_region_count) {
    throw error::Error("Number of regions {} and tasks executed {} are not the same!",
                       regions_list.size(),
                       progress.total_region_count);
  }
}

StrMap<RefRegion> ComputeBamFeatures::GetAllRefRegions(const vec<fs::path>& bam_paths) {
  StrMap<RefRegion> result;
  for (const auto& path : bam_paths) {
    const auto reader = _alignment_reader_cache.Open(path, "r");
    for (const auto& [chrom, ref_region] : GetRefRegions(reader)) {
      const auto itr = result.find(chrom);
      if (itr == result.end()) {
        result.try_emplace(chrom, ref_region);
      } else {
        // If the contig is already present, update the start and end positions if necessary
        itr->second.start_position = std::min(itr->second.start_position, ref_region.start_position);
        itr->second.end_position = std::max(itr->second.end_position, ref_region.end_position);
        itr->second.length = std::max(itr->second.length, ref_region.length);
      }
    }
  }
  return result;
}

u32 ComputeBamFeatures::ParallelComputeBamFeatures(const ComputeBamFeaturesCliParams& cli_params) {
  // This is the top-level method for computing bam features. All calls to compute bam features, whether from the
  // compute-bam-features executable or from filter-variants go through here.

  // Use a vector of Intervals per chromosome vs BedRegion structs for storing BED regions. Makes lookups easier for
  // getting which intervals overlap a given position on a chromosome.
  const auto bed_regions = ToIntervals(cli_params.bed_regions);

  const auto contig_map = GetAllRefRegions(cli_params.bam_input);

  // For each contig/chrom with coverage we get the reference sequence from the reference genome
  LoadReferenceSequences(cli_params.genome, bed_regions, contig_map);

  // Break up the BAM into a set of regions for parallel processing. If we have BED regions these are used, otherwise
  // WGS is assumed and the regions are generated based on coverage and blocksize
  vec<Region> regions_list;
  PopulateRegionsQueue(contig_map, bed_regions, regions_list, cli_params.max_region_size_per_thread);

  // Setup logging for feature computation progress
  Progress progress(regions_list.size());

  // Adjust minimum family size if using duplex data
  u32 min_family_size = cli_params.min_family_size;
  if (IsDuplexProtocol(cli_params.sequencing_protocol) && cli_params.min_family_size > kDuplexReadFamilySize) {
    WarnAsErrorIfSet("Reduce min family size from {} to {} for duplex data.", min_family_size, kDuplexReadFamilySize);
    min_family_size = kDuplexReadFamilySize;
  }

  const ComputeBamFeaturesParams params{
      .feature_cols = cli_params.config.feature_cols,
      .min_bq = cli_params.min_bq,
      .min_mapq = cli_params.min_mapq,
      .min_allowed_distance_from_end = cli_params.min_allowed_distance_from_end,
      .min_family_size = min_family_size,
      .max_read_variant_count = cli_params.max_read_variant_count,
      .max_read_variant_count_normalized = cli_params.max_read_variant_count_normalized,
      .min_homopolymer_length = cli_params.min_homopolymer_length,
      .sequencing_protocol = cli_params.sequencing_protocol,
      .filter_homopolymer = cli_params.filter_homopolymer,
      .tumor_sample_name = cli_params.tumor_sample_name,
      .tumor_rg_ids = GetReadGroupIdsForSample(_alignment_reader_cache.Open(cli_params.bam_input, "r"),
                                               cli_params.tumor_sample_name),
      .decode_yc = cli_params.decode_yc,
      .min_base_type = cli_params.min_base_type,
  };

  SetupWorkersAndAlignmentReaders(cli_params, params);

  ExecuteBamFeatureExtraction(regions_list, cli_params, progress);

  return _bases_covered;
}

std::optional<StrUnorderedSet> GetReadGroupIdsForSample(const vec<AlignmentReader>& readers,
                                                        const std::optional<std::string>& sample_name) {
  if (!sample_name.has_value() || sample_name->empty()) {
    return std::nullopt;
  }
  if (readers.empty()) {
    throw error::Error("No alignment readers available to extract read groups");
  }
  StrUnorderedSet result;
  io::Kstring ks;
  for (const auto& reader : readers) {
    for (s32 i = 0; i < sam_hdr_count_lines(reader.hdr.get(), "RG"); ++i) {
      const char* const rg_id = sam_hdr_line_name(reader.hdr.get(), "RG", i);
      if (rg_id == nullptr) {
        continue;
      }
      if (sam_hdr_find_tag_pos(reader.hdr.get(), "RG", i, "SM", ks.Get()) != 0) {
        continue;
      }
      if (ks.Get() != nullptr && sample_name == ks.Str()) {
        result.emplace(rg_id);
      }
    }
  }
  return result;
}

void ComputeBamFeatures::SetupWorkersAndAlignmentReaders(const ComputeBamFeaturesCliParams& cli_params,
                                                         const ComputeBamFeaturesParams& params) {
  _alignment_readers.resize(cli_params.threads);
  for (u32 i = 0; i < cli_params.threads; ++i) {
    _alignment_readers[i] = _alignment_reader_cache.Open(cli_params.bam_input, "r");
  }

  if (cli_params.output_file.has_value()) {
    CreateParentDirectoryIfNotExists(cli_params.output_file.value());
    _writer = std::make_unique<LockedTsvWriter>(*cli_params.output_file);

    if (cli_params.command_line.has_value()) {
      // write version and command line as comment lines
      io::Comments cmt;
      const auto& cmd_info = cli_params.command_line.value();
      io::AddVersionAndCommandLineComment(cmt, cmd_info.version, cmd_info.command_line);
      _writer->AppendComments(cmt);
    }
    // write header row
    _writer->AppendRow(cli_params.config.bam_feature_names);
  }

  _workers.reserve(cli_params.threads);
  for (u32 i = 0; i < cli_params.threads; ++i) {
    _workers.emplace_back(_alignment_readers.at(i), params, _writer.get());
  }
}

void ComputeBamFeatures::ExtractBamFeatures(s32 thread_id,
                                            const Region& region,
                                            const ComputeBamFeaturesCliParams& cli_params,
                                            Progress& progress) {
  try {
    if (thread_id < 0 || std::cmp_greater_equal(thread_id, cli_params.threads)) {
      throw error::Error("Invalid thread id");
    }

    const auto skip_variants = cli_params.skip_variants_vcf
                                   ? std::make_optional(ExtractVariantKeySet(*cli_params.skip_variants_vcf, region))
                                   : std::nullopt;

    _workers.at(thread_id).ComputeBamFeatures(region, _ref_seqs[region.chrom], std::nullopt, skip_variants);

    progress.UpdateAndLog(log::LogLevel::kInfo);
    _bases_covered += region.end - region.start;
  } catch (const std::exception& e) {
    throw error::Error(
        "Error extracting BAM features in '{}:{}-{}': {}", region.chrom, region.start, region.end, e.what());
  }
}

void ComputeBamFeatures::PopulateRegionsQueue(const StrMap<RefRegion>& contig_map,
                                              const StrMap<vec<Interval>>& bed_regions,
                                              vec<Region>& regions_queue,
                                              const u64 max_region_size_per_thread) {
  if (!bed_regions.empty()) {
    GetGenomicRegionsForBedFile(contig_map, bed_regions, regions_queue, max_region_size_per_thread);
  } else {
    GetGenomicRegionsForChromosome(contig_map, regions_queue, max_region_size_per_thread);
  }
}

void ComputeBamFeatures::GetGenomicRegionsForBedFile(const StrMap<RefRegion>& contig_map,
                                                     const StrMap<vec<Interval>>& bed_regions,
                                                     vec<Region>& regions_queue,
                                                     const u64 max_region_size_per_thread) {
  // This function breaks up BED entries per contig into smaller regions of at most a set size that can be used for
  // parallel processing.
  std::optional<std::string> prev_chrom = std::nullopt;
  std::optional<Interval> prev_interval = std::nullopt;
  for (const auto& [chrom, intervals] : bed_regions) {
    // skip if the contig has no alignments
    if (!contig_map.contains(chrom)) {
      continue;
    }

    if (!prev_chrom.has_value() || prev_chrom.value() != chrom) {
      // this is the first interval in the chromosome; no previous interval
      prev_interval = std::nullopt;
    }
    prev_chrom = chrom;

    // Add intervals to the list of regions. If the interval is too big break up the interval into smaller chunks
    for (const auto& interval : intervals) {
      const u64 region_start = interval.start;
      const u64 region_end = interval.end;
      for (u64 start = region_start; start <= region_end + 1; start += max_region_size_per_thread) {
        const u64 end = std::min(start + max_region_size_per_thread, region_end + 1);
        regions_queue.emplace_back(chrom, start, end, prev_interval);
        prev_interval = Interval{start, end};
      }
    }
  }
}

void ComputeBamFeatures::GetGenomicRegionsForChromosome(const StrMap<RefRegion>& contig_map,
                                                        vec<Region>& regions_queue,
                                                        const u64 max_region_size_per_thread) {
  // This function breaks up contigs into smaller regions up to a certain size, which are used to compute features in
  // parallel
  std::optional<std::string> prev_chrom = std::nullopt;
  std::optional<Interval> prev_interval = std::nullopt;
  for (const auto& [chrom, region] : contig_map) {
    // leftmost aligned position of the contig
    const u64 region_start = region.start_position;
    // rightmost aligned position of the contig
    const u64 region_end = region.end_position;

    if (!prev_chrom.has_value() || prev_chrom.value() != chrom) {
      // this is the first interval in the chromosome; no previous interval
      prev_interval = std::nullopt;
    }
    prev_chrom = chrom;

    // Add regions to the list up to the max region size. Individual threads will run BAM feature computation per each
    // of these regions.
    for (u64 start = region_start; start <= region_end + 1; start += max_region_size_per_thread) {
      const u64 end = std::min(start + max_region_size_per_thread, region_end + 1);
      regions_queue.emplace_back(chrom, start, end, prev_interval);
      prev_interval = Interval{start, end};
    }
    // TODO: add dynamic partitioning that breaks up dense/high coverage regions
  }
}

void ComputeBamFeatures::LoadReferenceSequences(const fs::path& genome,
                                                const StrMap<vec<Interval>>& bed_regions,
                                                const StrMap<RefRegion>& contig_map) {
  io::FastaReader fasta_reader(genome);
  for (const auto& [chrom, ref_region] : contig_map) {
    if (ref_region.length == 0) {
      continue;
    }
    if (!bed_regions.empty() && !bed_regions.contains(chrom)) {
      // If a BED file is provided, do not store the sequences for chromosomes not listed in the BED file.
      continue;
    }

    _ref_seqs[chrom] = fasta_reader.GetSequence(chrom, 0, static_cast<s64>(ref_region.length));
  }
}

}  // namespace xoos::svc
