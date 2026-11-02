#pragma once

#include <concepts>
#include <memory>
#include <utility>

#include <htslib/sam.h>

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/vec.h>

#include "core/read-collapser-options.h"
#include "core/region.h"
#include "core/thread-context.h"
#include "io/alignment-io.h"
#include "io/alignment.h"

namespace xoos::read_collapser {

/**
 * A predicate that decides whether the region currently being processed owns a given alignment record, and should
 * therefore keep it rather than skip it.
 *
 * The consensus and markdup paths share the same batch-read loop but diverge on ownership strategy:
 *  - consensus uses a closest-start rule (RegionLookupTable::FindClosestOverlappingRegion) to preserve cluster
 *    integrity at region boundaries;
 *  - markdup uses start-containment (ShouldSkipAlignment) so per-thread BAMs concatenate into a globally sorted file.
 *
 * The chosen strategy is injected as a callable so each path supplies its own rule.
 *
 * PRECONDITION: implementations may rely on assumptions about how regions are partitioned. In particular, the
 * markdup start-containment strategy (ShouldSkipAlignment) inspects only the immediately prior region. This is
 * correct as long as regions are sorted and non-overlapping: a read overlapping several regions overlaps a
 * contiguous run of them, so checking the immediate prior suffices to skip already-owned reads.
 */
template <typename F>
concept ReadOwnershipPredicate = std::predicate<F, const bam1_t*>;

/**
 * A per-batch action invoked with the accumulated alignments each time a batch is flushed (either the batch size is
 * reached or the region is exhausted). The alignments are passed by mutable reference because clustering assigns
 * cluster information to them in place.
 *
 * The two paths supply different actions:
 *  - consensus: cluster, optionally emit a cluster BAM, then call consensus and write the FASTQ results;
 *  - markdup: cluster, mark duplicates, optionally collect supplementary duplicate names, then write the alignments.
 */
template <typename F>
concept BatchProcessor = std::invocable<F, vec<AlignmentPtr>&>;

/**
 * Shared batch-read loop for a single region, used by both the consensus and markdup paths.
 *
 * Reads alignments overlapping @p region from @p context's reader in batches of at most `options.batch_size`. For each
 * record it consults @p owns_read to decide ownership, constructs an Alignment, applies the common ShouldDiscardRecord
 * filtering, and accumulates surviving alignments into `context.GetMetrics()`. Whenever the batch is full or the region
 * is exhausted, the accumulated alignments are handed to @p process_batch and the buffer is cleared.
 *
 * Only the mapped batch-read loop is unified here. Region partitioning, cross-region ownership strategy, unmapped-read
 * handling, and taskflow orchestration remain mode-specific and are supplied by the caller (the former two via
 * @p owns_read, the latter two around this call).
 *
 * @param context Per-thread resources; its reader is read from and its metrics accumulated into.
 * @param options Run options; supplies batch size and the alignment-construction / filtering parameters.
 * @param region The mapped region to read alignments from. Callers must handle unmapped (HTS_IDX_NOCOOR) regions
 *               separately; this loop assumes a coordinate-indexed query.
 * @param owns_read Ownership predicate; returns true if the current region owns (keeps) the record.
 * @param process_batch Per-batch action invoked with the accumulated alignments before the buffer is cleared.
 */
template <ReadOwnershipPredicate OwnsReadFn, BatchProcessor ProcessBatchFn>
void ProcessRegionInBatches(ThreadContext& context,
                            const ReadCollapserOptions& options,
                            const Region& region,
                            OwnsReadFn owns_read,
                            ProcessBatchFn process_batch) {
  const auto itr = io::SamItrQueryI(context.Reader().idx.get(), region.tid, region.start, region.end);
  vec<AlignmentPtr> alignments;
  // Read alignments for the region in batches.
  // When the batch size is reached or there are no more reads, hand the batch to process_batch before continuing.
  while (true) {
    auto record =
        io::Bam1Ptr{bam_init1()};  // NOSONAR(M23_286) — Bam1Ptr has a custom deleter; make_unique cannot apply
    const bool has_more_reads = io::SamItrNext(context.Reader().bam_file.get(), itr.get(), record.get());
    const bool reached_batch_size_limit = (alignments.size() >= options.batch_size);
    if ((!has_more_reads || reached_batch_size_limit) && !alignments.empty()) {
      // Flush the accumulated batch and reset the buffer for the next one.
      process_batch(alignments);
      alignments.clear();
    }
    if (!has_more_reads) {
      break;
    }
    // An alignment overlapping the region boundary must be processed by exactly one region. Skip records that the
    // current region does not own, per the injected ownership strategy, to avoid processing a read more than once.
    if (!owns_read(record.get())) {
      continue;
    }
    const auto alignment_ptr =
        std::make_shared<Alignment>(std::move(record), options.cluster_by_umi, options.ignore_read_name_parsing_errors);
    if (ShouldDiscardRecord(options, alignment_ptr.get(), context.GetMetrics())) {
      continue;
    }
    alignments.emplace_back(alignment_ptr);
  }
}

}  // namespace xoos::read_collapser
