#pragma once

#include <concepts>
#include <memory>
#include <optional>
#include <string>

#include <ankerl/unordered_dense.h>
#include <htslib/sam.h>

#include <xoos/io/alignment-reader.h>
#include <xoos/io/fasta-reader.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/int.h>

#include "alignment-metrics-options.h"
#include "core/alignment.h"
#include "core/pileup.h"
#include "core/region-lookup.h"
#include "core/super-region.h"
#include "metadata/dataset-metadata.h"
#include "metrics/metrics.h"
#include "metrics/stratified-metrics.h"

namespace xoos::alignment_metrics {

constexpr auto kUnknownChromosome = "*";

using RefPosQueryPosPair = std::pair<s64, u32>;

/**
 * A concept that defines a callable object that can generate alignment records when called.
 *
 * The callable should be invocable with a `bam1_t*` argument where the alignment record will be stored,
 * and it should return a `s32` indicating the success or failure of the operation.
 *
 * `AlignmentProvider` provides a consistent interface for retrieving alignment records from various sources
 * such as BAM files, in-memory data structures, or other custom sources, allowing the `MetricsWorker` to be
 * flexible without hard-coded dependency on a BAM file.
 */
template <typename F>
concept AlignmentProvider = std::invocable<F, bam1_t*> && std::same_as<std::invoke_result_t<F, bam1_t*>, s32>;

/**
 * The default alignment provider that retrieves alignments from a BAM file iterator.
 */
class BamItrAlignmentProvider {
 public:
  explicit BamItrAlignmentProvider(const SuperRegion& super_region, const io::AlignmentReader& alignment_reader);
  s32 operator()(bam1_t* bam1_ptr) const;

 private:
  // Non-owning reference to a reader owned by the caller for the provider's lifetime; stored by
  // reference by design, making the provider intentionally non-assignable.
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-const-or-ref-data-members)
  const io::AlignmentReader& _alignment_reader;
  io::HtsItrPtr _itr;
};

/**
 * A worker class that processes alignments and calculates metrics for a given super region.
 * The results are stored in a `ConsensusAccuracyMetrics` object that can be retrieved by calling `GetMetrics()`
 * after processing is complete.
 */
class MetricsWorker {
 public:
  /**
   * Initialize a new `MetricsWorker` instance with the given options, dataset metadata, and region lookup table.
   * Once initialzed, the worker can be reused to process multiple super regions by calling `ProcessRegion()` for each
   * region.
   */
  MetricsWorker(const AlignmentMetricsOptions& options,
                const DatasetMetadata& dataset_metadata,
                const RegionLookupTable& region_lookup_table);

  /**
   * Process a given super region and calculate metrics for all alignments that overlap with the subregions
   * within the super region. This function reads data from the BAM file and FASTA file provided through
   * `AlignmentMetricsOptions` when the worker is initialized.
   */
  void ProcessRegion(const SuperRegion& super_region, size_t super_region_id);

  /**
   * Generic processing function for a region using a custom alignment provider and reference sequence.
   *
   * @param super_region The super region to process.
   * @param super_region_id The ID of the super region.
   * @param alignment_provider The alignment provider to use for retrieving alignments.
   * @param ref_seq The reference sequence to use for the region. All bases should be uppercase.
   */
  void ProcessRegion(const SuperRegion& super_region,
                     const size_t super_region_id,
                     AlignmentProvider auto& alignment_provider,
                     std::optional<std::string>&& ref_seq) {
    Initialize(super_region, super_region_id, std::move(ref_seq));
    // Iterate through all alignments provided by the callback
    const io::Bam1Ptr bam1_ptr(bam_init1());
    while (true) {
      const auto sam_itr_res = alignment_provider(bam1_ptr.get());
      if (sam_itr_res == -1) {
        // No more alignments in the region
        break;
      }
      // Update the current subregion index based on the alignment position
      if (_current_subregion_index < super_region.subregions.size()) {
        const auto& current_subregion = super_region.subregions.at(_current_subregion_index);
        if (bam1_ptr->core.pos >= current_subregion.end) {
          while (_current_subregion_index + 1 < super_region.subregions.size() &&
                 bam1_ptr->core.pos >= super_region.subregions.at(_current_subregion_index + 1).end) {
            ++_current_subregion_index;
          }
        }
      }
      // Create an alignment object to hold the alignment pointer, metadata, and qualities
      Alignment alignment(bam1_ptr.get(), _options);
      // Process the alignment
      ProcessAlignment(alignment);
    }
    // After processing all alignments, finalize the metrics for the region. Aggregation routes each
    // valid pileup / subregion into the label containers indicated by its label bit vector.
    AggregateMetricsFromPileups(
        _pileups, super_region, *_labeled_metrics, _options.min_depth, _options.max_alt_allele_fraction);
    // Update mean coverage histogram for coverage uniformity metrics
    UpdateCoverageUniformityHistogram(_pileups, super_region, *_labeled_metrics);
  }

  /**
   * Process a region without coverage, fetching its reference sequence when configured.
   */
  void ProcessNoCoverageRegion(const std::string& chromosome, const Interval& region) const;

  /**
   * Process a region with no coverage and update coverage histogram
   */
  void ProcessNoCoverageRegion(const Interval& region, std::optional<std::string> ref_seq) const;

  /**
   * Process a single alignment and update the metrics based on this alignment.
   */
  void ProcessAlignment(Alignment& alignment);

  /**
   * Handle a match CIGAR operation for the given alignment, updating pileup entries and reference/query positions
   * accordingly.
   */
  RefPosQueryPosPair HandleMatchCigar(Alignment& alignment, s64 ref_pos, u32 query_pos, u32 cigar_len);

  /**
   * Handle an insertion CIGAR operation for the given alignment, updating pileup entries and reference/query positions
   * accordingly.
   */
  RefPosQueryPosPair HandleInsertionCigar(
      Alignment& alignment, s64 ref_pos, u32 query_pos, u32 cigar_len, bool is_last_cigar_op);

  /**
   * Handle a deletion CIGAR operation for the given alignment, updating pileup entries and reference/query positions
   * accordingly.
   */
  RefPosQueryPosPair HandleDeletionCigar(Alignment& alignment, s64 ref_pos, u32 query_pos, u32 cigar_len);

  /**
   * Handle homopolymer errors for the given alignment, updating the homopolymer error statistics accordingly.
   * This method is called after processing all CIGAR operations for the alignment to aggregate homopolymer error
   * information within the `Alignment` object.
   */
  void HandleHomopolymerErrors(const Alignment& alignment);

  /**
   * Retrieve the first (label index 0) metrics container calculated by this worker. In non-stratified
   * mode this is the sole container; in stratified mode it is the container for the first label. Most
   * callers should prefer GetStratifiedMetrics().
   */
  std::shared_ptr<Metrics> GetMetrics() const {
    return _metrics;
  }

  /**
   * Retrieve the per-label metrics containers calculated by this worker.
   */
  std::shared_ptr<StratifiedMetrics> GetStratifiedMetrics() const {
    return _labeled_metrics;
  }

  /**
   * Retrieve the pileups generated during processing of the current region.
   * Useful for inspecting the per-position metrics and coverage information
   * for debugging and testing purposes.
   */
  vec<Pileup> GetPileups() const {
    return _pileups;
  }

  /**
   * Initialize the worker to process a new super region, resetting internal state as needed.
   */
  void Initialize(const SuperRegion& super_region, size_t super_region_id, std::optional<std::string>&& ref_seq);

 private:
  /**
   * Increment total and spanning homopolymer coverage for the homopolymer positions covered by this read, and mark the
   * read's per-homopolymer error profiles that qualify for spanning coverage.
   */
  void UpdateHomopolymerSpanningCoverage(Alignment& alignment, s64 ref_pos);

  /**
   * Record the entire alignment to read-level metrics
   */
  void RecordRead(const Alignment& alignment, bool passes_read_filter) const;
  /**
   * Record the entire alignment to GC bias metrics.
   */
  void RecordReadForGcBias(const Alignment& alignment, bool passes_read_filter) const;

  /**
   * Process each CIGAR operation in the alignment, updating pileup entries and reference/query positions accordingly.
   * This method calls the appropriate handler for each CIGAR operation type (match, insertion, deletion, etc.).
   * and records the relevant information in the pileup.
   */
  void ProcessCigarOps(Alignment& alignment, s64 ref_pos, u32 query_pos);

  /**
   * Record a single matched base into the given pileup.
   * `query_offset` is the offset of the base within the current CIGAR operation (added to `query_pos`).
   */
  void RecordMatch(Pileup& pileup, const Alignment& alignment, s64 rpos, u32 query_pos, u32 query_offset);

  /**
   * Record homopolymer error information for a single matched base into the pileup's homopolymer error profile.
   * `query_offset` is the offset of the base within the current CIGAR operation (added to `query_pos`).
   */
  void RecordMatchHomopolymerError(
      const Pileup& pileup, Alignment& alignment, s64 rpos, u32 query_pos, u32 query_offset);

  /**
   * Record insertion into the given pileup.
   */
  void RecordInsertion(
      Pileup& pileup, const Alignment& alignment, const std::string& inserted_seq, u32 query_pos, u32 cigar_len);

  /**
   * Record homopolymer error information for an insertion into the pileup's homopolymer error profile.
   */
  void RecordInsertionHomopolymerError(const Pileup& pileup,
                                       Alignment& alignment,
                                       s64 ref_pos,
                                       const std::string& inserted_seq,
                                       u32 query_pos,
                                       u32 cigar_len) const;

  /**
   * Record deletion into the given pileup.
   */
  void RecordDeletion(Pileup& pileup, const Alignment& alignment, s64 ref_pos, u32 query_pos, u32 cigar_len);

  /**
   * Record homopolymer error information for a deletion into the pileup's homopolymer error profile.
   */
  void RecordDeletionHomopolymerError(const Pileup& pileup, Alignment& alignment, s64 ref_pos, u32 query_pos);

  /**
   * A pileup vector that holds depth, error, and homopolymer error profiles for each position in the current super
   * region. Only positions that are part of the subregions (regions of interest) within the super region have their
   * `valid` field set to true and are considered for metric calculations.
   */
  vec<Pileup> _pileups;
  // Indexes to track the current super region being processed
  size_t _current_super_region_index{0};
  // Indexes to track the current subregion within the super region being processed
  size_t _current_subregion_index{0};
  // The current super region being processed
  SuperRegion _current_super_region{};
  // The reference sequence for the current super region, used for error calculations
  std::optional<std::string> _ref_seq;
  // These are non-owning references/handles to state owned for the lifetime of the worker by the caller;
  // storing them by reference is a deliberate design choice, so the worker is intentionally non-assignable.
  // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
  // The options for calculating metrics
  const AlignmentMetricsOptions& _options;
  // Metadata about the dataset being processed
  const DatasetMetadata& _dataset_metadata;
  // A lookup table to quickly whether a read has already been processed when processing an earlier super region
  const RegionLookupTable& _region_lookup_table;
  // The alignment reader used to read alignments from the input BAM/CRAM file
  const std::shared_ptr<io::AlignmentReader> _alignment_reader;
  // One FASTA reader per worker, reused across all tasks assigned to this worker.
  std::unique_ptr<io::FastaReader> _fasta_reader;
  // The per-label metrics containers that hold the calculated metrics (one per label, or a single
  // container in non-stratified mode).
  const std::shared_ptr<StratifiedMetrics> _labeled_metrics;
  // Aliasing pointer to the label-0 container inside `_labeled_metrics`, kept for the GetMetrics()
  // convenience accessor and existing single-container tests. Valid because the container vector is
  // not resized after construction.
  const std::shared_ptr<Metrics> _metrics;
  // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)
};

}  // namespace xoos::alignment_metrics
