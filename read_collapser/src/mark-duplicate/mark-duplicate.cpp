#include "mark-duplicate/mark-duplicate.h"

#include <algorithm>
#include <ranges>

#include <htslib/sam.h>

#include <taskflow/taskflow.hpp>

#include <xoos/error/error.h>
#include <xoos/io/alignment-reader.h>
#include <xoos/io/bam/read-group.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/log/logging.h>
#include <xoos/types/fs.h>
#include <xoos/util/math.h>
#include <xoos/util/tmp-file.h>

#include "core/read-collapser-options.h"
#include "core/region-processor.h"
#include "core/region.h"
#include "core/thread-context.h"
#include "io/alignment-io.h"
#include "io/output-path.h"
#include "mark-duplicate/mean-base-quality.h"
#include "metrics/metrics.h"

namespace xoos::read_collapser {

/**
 * Phase 2: Re-read a temp BAM file and write to the final output, marking supplementary alignments
 * whose primary was marked as duplicate. The temp file and its index are deleted after processing.
 *
 * @param temp_bam_path Path to the temp BAM file from Phase 1.
 * @param final_bam_path Path to the final output BAM file.
 * @param duplicate_read_names Set of qnames of primary duplicates with SA tags.
 * @param remove_duplicates If true, supplementary duplicates are omitted from output.
 * @param write_index Whether to write an index for the final BAM file.
 * @return The number of supplementary alignments marked as duplicate.
 */
namespace {

u64 MarkSupplementaryAlignments(const fs::path& temp_bam_path,
                                const fs::path& final_bam_path,
                                const DuplicateReadNames& duplicate_read_names,
                                const bool remove_duplicates,
                                const bool write_index,
                                const bool remove_pg_tag) {
  // RAII guards ensure temp BAM and its index are removed after processing
  const TmpFile temp_bam(temp_bam_path);
  const TmpFile temp_bai(fs::path{temp_bam_path.string() + ".bai"});

  const auto reader = io::OpenAlignmentReader(temp_bam.Path());
  const auto writer = OpenAlignmentWriter(final_bam_path, reader.hdr.get(), write_index);

  u64 supplementary_duplicates_marked = 0;
  auto record = io::Bam1Ptr{bam_init1()};

  while (sam_read1(reader.bam_file.get(), reader.hdr.get(), record.get()) >= 0) {
    const bool is_supplementary = (record->core.flag & BAM_FSUPPLEMENTARY) != 0;

    if (is_supplementary && duplicate_read_names.contains(bam_get_qname(record.get()))) {
      record->core.flag |= BAM_FDUP;
      ++supplementary_duplicates_marked;
    }

    if (remove_pg_tag) {
      RemovePgTag(record.get());
    }

    if (!remove_duplicates || (record->core.flag & BAM_FDUP) == 0) {
      io::SamWrite1(writer.bam.get(), writer.bam->bam_header, record.get());
    }
  }

  return supplementary_duplicates_marked;
}

/**
 * Phase 2: Merge per-thread duplicate read name sets, then re-read each temp BAM file in parallel to propagate
 * the duplicate flag to supplementary alignments whose primary was marked as duplicate.
 *
 * @param options Configuration for duplicate marking.
 * @param super_region_count Number of super regions (temp BAM files).
 * @param per_thread_dup_names Per-thread sets of duplicate read names collected during Phase 1.
 * @param executor Taskflow executor for parallel processing.
 * @param metrics The metrics to record the supplementary duplicate count into.
 */
void PropagateSupplementaryDuplicates(const ReadCollapserOptions& options,
                                      const size_t super_region_count,
                                      vec<DuplicateReadNames>& per_thread_dup_names,
                                      tf::Executor& executor,
                                      Metrics& metrics) {
  DuplicateReadNames merged_dup_names;
  for (auto& thread_set : per_thread_dup_names) {
    merged_dup_names.merge(thread_set);
    // Free each per-thread set immediately to avoid holding both the merged and per-thread copies in memory
    thread_set = DuplicateReadNames{};
  }
  vec<DuplicateReadNames>{}.swap(per_thread_dup_names);

  Logging::Info("Marking supplementary duplicates ({} primary duplicates with SA tags)", merged_dup_names.size());

  tf::Taskflow supp_taskflow;
  auto supp_dup_counts = vec<u64>(super_region_count, 0);

  for (u32 i = 0; i < super_region_count; ++i) {
    const auto temp_bam_path = MarkdupRegionTempBamPath(options, i);
    if (!fs::exists(temp_bam_path)) {
      continue;
    }
    const auto final_bam_path = MarkdupRegionBamPath(options, i);
    const bool write_index = !options.merge_output;

    supp_taskflow.emplace([temp_bam_path,
                           final_bam_path,
                           &merged_dup_names,
                           &supp_dup_counts,
                           i,
                           remove_dups = options.markdup_options.remove_duplicates,
                           write_index,
                           remove_pg_tag = !options.markdup_options.keep_pg_tag] {
      supp_dup_counts.at(i) = MarkSupplementaryAlignments(
          temp_bam_path, final_bam_path, merged_dup_names, remove_dups, write_index, remove_pg_tag);
    });
  }

  executor.run(supp_taskflow).get();

  u64 total_supp_dups = 0;
  for (const auto count : supp_dup_counts) {
    total_supp_dups += count;
  }
  metrics.duplicate_supplementary_alignments = total_supp_dups;
}

}  // namespace

u32 MarkDuplicate(const ReadCollapserOptions& options) {
  // Create the same number of super regions as the number of threads
  // so that the BAM output produced by each thread is already sorted
  // and we can just concatenate the output files to get a single sorted BAM file.
  auto super_regions = DetermineSuperRegions(options, options.threads);
  const auto super_region_count = super_regions.size();

  // Size at least 1 so the Phase 2 supplementary-duplicate count always has a metrics to write into.
  vec<Metrics> super_region_metrics(std::max<size_t>(super_region_count, 1));

  // When mark_supplementary_alignments is enabled, Phase 1 writes to temp files and Phase 2 produces the final output.
  // Otherwise, Phase 1 writes directly to the final output.
  const bool mark_supplementary = options.markdup_options.mark_supplementary_alignments;
  // Only write index for the per-thread BAM files if we are not merging them later or we need to mark supplementary
  // alignments in Phase 2 because we need to read the per-thread BAM files again in Phase 2.
  const bool write_index = mark_supplementary || !options.merge_output;

  // Per-thread sets of duplicate read names with SA tags, used only when mark_supplementary_alignments is enabled.
  // Each thread gets its own set to avoid contention; sets are merged after Phase 1.
  auto per_thread_dup_names = vec<DuplicateReadNames>(super_region_count);

  // Executor is shared by Phase 1 and Phase 2.
  tf::Executor executor{options.threads};

  // Phase 1: mark duplicates in parallel across super regions. One context per super region (task i -> context i ->
  // output.i.bam) preserves the sorted-concatenation guarantee. Contexts are scoped so their BAM writers flush and
  // index by RAII before Phase 2 re-reads the temp BAMs.
  {
    vec<ThreadContext> contexts;
    contexts.reserve(super_region_count);
    for (size_t i = 0; i < super_region_count; ++i) {
      auto& context = contexts.emplace_back(options);
      const auto bam_filename =
          mark_supplementary ? MarkdupRegionTempBamPath(options, i) : MarkdupRegionBamPath(options, i);
      context.OpenBamWriter(bam_filename, write_index);
    }

    tf::Taskflow taskflow;
    for (u32 super_region_id = 0; super_region_id < super_region_count; ++super_region_id) {
      auto& context = contexts.at(super_region_id);
      auto& super_region = super_regions.at(super_region_id);
      auto& dup_names_set = per_thread_dup_names.at(super_region_id);

      taskflow.emplace([&context, &options, super_region_id, &super_region, &super_regions, &dup_names_set] {
        try {
          DuplicateMarkSuperRegion(context, options, super_region, super_region_id, super_regions, dup_names_set);
        } catch (const std::exception& e) {
          Logging::Error("Error processing super region '{}': {}", super_region_id, e.what());
          throw;
        }
      });
    }
    executor.run(taskflow).get();

    // Snapshot per-context metrics before the contexts (and their BAM writers) are destroyed at the end of this scope.
    for (size_t i = 0; i < super_region_count; ++i) {
      super_region_metrics.at(i) = contexts.at(i).GetMetrics();
    }
  }  // contexts destroyed here: all per-thread BAMs are flushed and indexed via AlignmentWriter RAII

  // Phase 2: Mark supplementary duplicates if enabled
  if (options.markdup_options.mark_supplementary_alignments) {
    // Phase 1 tasks have completed, so super_region_metrics 0 is no longer being written concurrently; record the
    // supplementary duplicate count there so it is included in the aggregated total.
    PropagateSupplementaryDuplicates(
        options, super_region_count, per_thread_dup_names, executor, super_region_metrics.at(0));
  }

  const auto metrics_total = SumMetrics(super_region_metrics);
  const io::TsvMetadata tsv_metadata{.command_line_info = options.command_line_info,
                                     .read_group = io::GetFirstReadGroupFromBam(options.bam_input)};
  metrics_total.WriteAllMetricsToTsv(ReadCollapserMode::kMarkDuplicate, options, tsv_metadata);
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
namespace {

void MergeOutputBamFiles(const ReadCollapserOptions& options, const u32 super_region_count) {
  Logging::Info("Merging output BAM files");
  vec<fs::path> bam_files;
  for (size_t i = 0; i < super_region_count; ++i) {
    const auto region_bam = MarkdupRegionBamPath(options, i);
    if (fs::exists(region_bam)) {
      bam_files.emplace_back(region_bam);
    }
  }
  const auto merged_bam_filename = MarkdupMergedBamPath(options);
  ConcatenateBamFiles(bam_files, merged_bam_filename);
  // Index the merged BAM file
  const fs::path bai_filename = merged_bam_filename.string() + ".bai";
  // NOLINTNEXTLINE(readability-suspicious-call-argument) - bam/index args are correctly ordered
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
void DuplicateMarkCluster(const Clusters& clusters) {
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
 * For each primary alignment that is marked as duplicate and has an SA tag, insert its qname into
 * @p duplicate_read_names. Only primary (non-supplementary, non-secondary) alignments are considered.
 */
void CollectDuplicateReadNamesWithSaTag(const vec<AlignmentPtr>& alignments, DuplicateReadNames& duplicate_read_names) {
  for (const auto& alignment : alignments) {
    if (!alignment->is_duplicate) {
      continue;
    }
    const u16 flag = alignment->record->core.flag;
    // Only collect primary alignments (not supplementary or secondary)
    static constexpr u16 kNonPrimaryFlagMask =
        static_cast<u16>(static_cast<u32>(BAM_FSUPPLEMENTARY) | static_cast<u32>(BAM_FSECONDARY));
    if ((flag & kNonPrimaryFlagMask) != 0) {
      continue;
    }
    if (alignment->HasSaTag()) {
      duplicate_read_names.emplace(bam_get_qname(alignment->record.get()));
    }
  }
}

}  // namespace

bool ShouldSkipAlignment(const bam1_t* record,
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
  return record->core.tid == prior_region.tid && std::cmp_less(prior_region.start, bam_endpos(record)) &&
         std::cmp_less(record->core.pos, prior_region.end);
}

// Perform clustering and duplicate marking on a super region
void DuplicateMarkSuperRegion(ThreadContext& context,
                              const ReadCollapserOptions& options,
                              const SuperRegion& super_region,
                              const u32 super_region_id,
                              const vec<SuperRegion>& super_regions,
                              DuplicateReadNames& duplicate_read_names) {
  const auto& alignment_reader = context.Reader();
  const auto& alignment_writer = context.BamWriter();
  auto& metrics = context.GetMetrics();
  u32 clusters_id = 0;
  auto create_cluster_id = [&super_region_id, &clusters_id]() -> ClusterId {
    return {super_region_id, clusters_id++};  // NOLINT
  };

  // When mark_supplementary_alignments is enabled, don't remove duplicates during Phase 1 —
  // Phase 2 handles removal after supplementary reads have been marked.
  const bool remove_duplicates_phase1 =
      !options.markdup_options.mark_supplementary_alignments && options.markdup_options.remove_duplicates;

  const bool remove_pg_tag = !options.markdup_options.keep_pg_tag;

  // Per-batch action: cluster the alignments, mark duplicates, optionally collect supplementary duplicate names,
  // then write the alignments to the per-thread BAM.
  const auto process_batch = [&options,
                              &create_cluster_id,
                              &metrics,
                              &duplicate_read_names,
                              &remove_duplicates_phase1,
                              &alignment_writer,
                              remove_pg_tag](vec<AlignmentPtr>& alignments) {
    const Clusters clusters = ClusterAlignments(options, alignments, create_cluster_id, metrics);
    DuplicateMarkCluster(clusters);
    if (options.markdup_options.mark_supplementary_alignments) {
      CollectDuplicateReadNamesWithSaTag(alignments, duplicate_read_names);
    }
    WriteAlignments(alignments,
                    remove_duplicates_phase1,
                    !options.markdup_options.exclude_cluster_tags,
                    alignment_writer.bam,
                    remove_pg_tag);
  };

  for (size_t subregion_index = 0; subregion_index < super_region.size(); ++subregion_index) {
    const auto& region = super_region.at(subregion_index);
    // Treat each unmapped read as its own singleton cluster and directly write it out
    if (region.tid == HTS_IDX_NOCOOR) {
      const auto unmapped_read_writer = [&alignment_writer, remove_pg_tag](const ClusterId& cluster_id,
                                                                           const bam1_t* const record) {
        WriteAlignment(record, false, cluster_id, 1, alignment_writer.bam, remove_pg_tag);
      };
      ReadAndWriteUnmappedReads(options, create_cluster_id, alignment_reader, unmapped_read_writer, metrics);
      continue;
    }
    // Ownership strategy: start-containment. An alignment is owned by the region containing its start position, so a
    // read overlapping a prior region is skipped here (that region already processed it). This lets per-thread BAMs
    // concatenate into a globally sorted file. See ShouldSkipAlignment; its precondition is sorted, non-overlapping
    // regions.
    const auto owns_read = [&super_regions, &super_region_id, &subregion_index](const bam1_t* const record) {
      return !ShouldSkipAlignment(record, super_regions, super_region_id, subregion_index);
    };
    ProcessRegionInBatches(context, options, region, owns_read, process_batch);
  }
}

void MarkDuplicateAndMergeOutput(const ReadCollapserOptions& options) {
  ValidateOutputFilesDoNotExist(options, ReadCollapserMode::kMarkDuplicate);
  const auto super_region_count = MarkDuplicate(options);
  if (options.merge_output) {
    // Merging the per-thread BAM output files into a single BAM file
    // after finishing marking duplicates and closing the file handles
    MergeOutputBamFiles(options, super_region_count);
  }
  Logging::Info("Done");
}

}  // namespace xoos::read_collapser
