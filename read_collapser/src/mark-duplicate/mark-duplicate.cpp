#include "mark-duplicate/mark-duplicate.h"

#include <filesystem>
#include <ranges>

#include <htslib/sam.h>

#include <taskflow/taskflow.hpp>

#include <xoos/error/error.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/log/logging.h>
#include <xoos/util/math.h>

#include "core/read-collapser-options.h"
#include "core/region.h"
#include "io/alignment-io.h"
#include "mark-duplicate/mean-base-quality.h"
#include "metrics/metrics.h"

namespace xoos::read_collapser {

u32 MarkDuplicate(const ReadCollapserOptions& options) {
  // Reset Metrics instance so that test cases do not share the same data (due to the static Metrics instance member)
  ConcurrentMetrics::Reset();

  // Create the same number of super regions as the number of threads
  // so that the BAM output produced by each thread is already sorted
  // and we can just concatenate the output files to get a single sorted BAM file.
  auto super_regions = DetermineSuperRegions(options, options.threads);
  const auto super_region_count = super_regions.size();

  auto alignment_readers = OpenAlignmentReaders(options.bam_input, super_region_count);
  auto alignment_writers = vec<AlignmentWriter>();
  for (size_t i = 0; i < super_region_count; ++i) {
    const auto bam_filename = options.output_dir / fmt::format("output.{:04}.bam", i);
    const auto hdr = alignment_readers.at(i).hdr.get();
    io::SamHdrAddPgLine(
        hdr, io::PgHdrLine{options.program_name, options.program_name, options.version, options.command_line});
    // Only write index for the per-thread BAM files if we are not merging them later
    const bool write_index = !options.merge_output;
    alignment_writers.emplace_back(OpenAlignmentWriter(bam_filename, hdr, write_index));
  }

  tf::Taskflow taskflow;
  for (u32 super_region_id = 0; super_region_id < super_region_count; ++super_region_id) {
    auto& alignment_reader = alignment_readers.at(super_region_id);
    auto& alignment_writer = alignment_writers.at(super_region_id);

    auto& super_region = super_regions.at(super_region_id);
    std::function<void()> super_region_func =
        [&options, super_region_id, &super_region, &alignment_reader, &alignment_writer, &super_regions] {
          try {
            DuplicateMarkSuperRegion(
                options, super_region_id, super_region, alignment_reader, alignment_writer, super_regions);
          } catch (const std::exception& e) {
            Logging::Error("Error processing super region '{}': {}", super_region_id, e.what());
            throw;
          }
        };
    taskflow.emplace(super_region_func);
  }

  tf::Executor executor{options.threads};
  executor.run(taskflow).get();

  auto metrics_total = ConcurrentMetrics::GetTotal();
  metrics_total.WriteAllMetricsToTsv(options.output_dir,
                                     ReadCollapserMode::kMarkDuplicate,
                                     options.cluster_by_umi,
                                     options.version,
                                     options.command_line);
  return static_cast<u32>(super_region_count);
}

/**
 * Merge the per-thread BAM output files into a single BAM file and create its index.
 * Remove the temporary per-thread BAM files after merging.
 *
 * We assume that the per-thread BAM files are named as "output.XXXX.bam" under @p options.output_dir,
 * and are sorted within each file, and sorted across files based on the order of super regions specified.
 *
 * Note that file handles associated with the per-thread BAM files must be closed
 * before calling this function.
 *
 * @param options The ReadCollapserOptions containing output directory information.
 * @param super_region_count The number of super regions (i.e., the number of per-thread BAM files) including the super
 * region for unmapped reads.
 */
static void MergeOutputBamFiles(const ReadCollapserOptions& options, const u32 super_region_count) {
  Logging::Info("Merging output BAM files");
  vec<fs::path> bam_files;
  for (size_t i = 0; i < super_region_count; ++i) {
    if (fs::exists(options.output_dir / fmt::format("output.{:04}.bam", i))) {
      bam_files.emplace_back(options.output_dir / fmt::format("output.{:04}.bam", i));
    }
  }
  const auto merged_bam_filename = options.output_dir / "output.bam";
  ConcatenateBamFiles(bam_files, merged_bam_filename);
  // Index the merged BAM file
  const fs::path bai_filename = merged_bam_filename.string() + ".bai";
  io::SamIndexBuild3(merged_bam_filename, bai_filename, 0, static_cast<s32>(options.threads));
  Logging::Info("Removing temporary BAM files");
  for (const auto& bam_file : bam_files) {
    fs::remove(bam_file);
    // also remove the index file if it exists
    const fs::path bai_file = bam_file.string() + ".bai";
    if (fs::exists(bai_file)) {
      fs::remove(bai_file);
    }
  }
}

/**
 * Mark duplicates in each cluster in @p clusters. The alignment with the highest mean base quality in each cluster
 * is not marked as a duplicate, all other alignments in the cluster are marked as duplicates.
 */
static void DuplicateMarkCluster(const Clusters& clusters) {
  for (const auto& cluster : clusters | std::views::values) {
    // Update alignment duplicate status
    for (const auto& alignment : cluster->alignments) {
      alignment->is_duplicate = true;
    }
    if (auto alignment_not_duplicate = FindAlignmentWithMaxMeanBaseQ(cluster->alignments);
        alignment_not_duplicate != nullptr) {
      alignment_not_duplicate->is_duplicate = false;
    }
  }
}

/**
 * Determine if the alignment record should be skipped because it overlaps with a prior region.
 * This includes regions within the same super region or regions in a prior super region.
 *
 * @param record The alignment record to check.
 * @param super_regions The vector of all super regions.
 * @param super_region_id The ID of the current super region.
 * @param subregion_index The index of the current subregion within the super region.
 * @return True if the alignment overlaps with a prior region, false otherwise.
 */
static bool ShouldSkipAlignment(const bam1_t* record,
                                const vec<SuperRegion>& super_regions,
                                const size_t super_region_id,
                                const size_t subregion_index) {
  // Get the current region
  const auto& current_region = super_regions.at(super_region_id).at(subregion_index);

  // If the alignment starts within the current region, the current region should process it
  if (std::cmp_greater_equal(record->core.pos, current_region.start) &&
      std::cmp_less(record->core.pos, current_region.end)) {
    return false;
  }

  // Check for overlap with a PRIOR region
  const Region* prior_region_ptr = nullptr;
  if (subregion_index > 0) {
    // Check prior region within the SAME super region
    prior_region_ptr = &super_regions.at(super_region_id).at(subregion_index - 1);

  } else if (super_region_id > 0) {
    // Check the last region of the PRIOR super region
    prior_region_ptr = &super_regions.at(super_region_id - 1).back();
  }

  if (prior_region_ptr == nullptr) {
    return false;
  }

  // Skip the alignment if it overlaps with the prior region
  const auto& prior_region = *prior_region_ptr;
  if (record->core.tid == prior_region.tid && std::cmp_less(prior_region.start, bam_endpos(record)) &&
      std::cmp_less(record->core.pos, prior_region.end)) {
    return true;
  }

  return false;
}

// Perform clustering and consensus on a super region
void DuplicateMarkSuperRegion(const ReadCollapserOptions& options,
                              const u32 super_region_id,
                              const SuperRegion& super_region,
                              const AlignmentReader& alignment_reader,
                              const AlignmentWriter& alignment_writer,
                              const vec<SuperRegion>& super_regions) {
  auto& metrics = ConcurrentMetrics::Get();
  u32 clusters_id = 0;
  auto create_cluster_id = [&super_region_id, &clusters_id]() -> ClusterId {
    return {super_region_id, clusters_id++};  // NOLINT
  };

  for (size_t subregion_index = 0; subregion_index < super_region.size(); ++subregion_index) {
    const auto& region = super_region.at(subregion_index);
    // Treat each unmapped read as its own singleton cluster and directly write it out
    if (region.tid == HTS_IDX_NOCOOR) {
      const auto unmapped_alignment_writer = [&alignment_writer](const ClusterId& cluster_id,
                                                                 const bam1_t* const record) {
        WriteAlignment(record, false, cluster_id, 1, alignment_writer.bam);
      };
      ReadAndWriteUnmappedAlignments(options, create_cluster_id, alignment_reader, unmapped_alignment_writer);
      continue;
    }
    const auto itr = io::SamItrQueryI(alignment_reader.idx.get(), region.tid, region.start, region.end);
    vec<AlignmentPtr> alignments;
    // Read alignments for the region in batches
    // If the batch size is reached or there are no more reads to read, cluster and write the alignments
    while (true) {
      auto record = io::Bam1Ptr{bam_init1()};
      const bool has_more_reads = io::SamItrNext(alignment_reader.bam.get(), itr.get(), record.get());
      const bool reached_batch_size_limit = (alignments.size() >= options.batch_size);
      if ((!has_more_reads || reached_batch_size_limit) && !alignments.empty()) {
        Clusters clusters = ClusterAlignments(options, alignments, create_cluster_id);
        DuplicateMarkCluster(clusters);
        WriteAlignments(alignments, options.remove_duplicates, !options.exclude_cluster_tags, alignment_writer.bam);
        // Clear the alignments for the next batch
        alignments.clear();
      }
      if (!has_more_reads) {
        // No more records to read
        break;
      }
      // Skip alignments that are already covered by a prior region
      if (ShouldSkipAlignment(record.get(), super_regions, super_region_id, subregion_index)) {
        continue;
      }
      if (ShouldDiscardAlignment(options, record.get(), metrics)) {
        continue;
      }
      alignments.emplace_back(std::make_shared<Alignment>(std::move(record)));
    }
  }
}

void MarkDuplicateAndMergeOutput(const ReadCollapserOptions& options) {
  const auto super_region_count = MarkDuplicate(options);
  if (options.merge_output) {
    // Merging the per-thread BAM output files into a single BAM file
    // after finishing marking duplicates and closing the file handles
    MergeOutputBamFiles(options, super_region_count);
  }
  Logging::Info("Done");
}

}  // namespace xoos::read_collapser
