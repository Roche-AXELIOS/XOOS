#pragma once

#include <concepts>
#include <memory>
#include <string>

#include <ankerl/unordered_dense.h>
#include <htslib/sam.h>

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/int.h>

#include "alignment-metrics-options.h"
#include "core/alignment.h"
#include "core/pileup.h"
#include "core/region-lookup.h"
#include "core/super-region.h"
#include "io/alignment-reader.h"
#include "metadata/dataset-metadata.h"
#include "metrics/metrics.h"

namespace xoos::alignment_metrics {

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
  explicit BamItrAlignmentProvider(const SuperRegion& super_region, const AlignmentReader& alignment_reader);
  s32 operator()(bam1_t* bam1_ptr) const;

 private:
  const AlignmentReader& _alignment_reader;
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
   * @param ref_seq The reference sequence to use for the region.
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
    // After processing all alignments, finalize the metrics for the region
    AggregateMetricsFromPileups(_pileups, *_metrics, _options.min_depth, _options.max_alt_allele_fraction);
    // Update mean coverage histogram for coverage uniformity metrics
    UpdateCoverageUniformityHistogram(_pileups, super_region, *_metrics);
  }

  /**
   * Process a region with no coverage and update coverage histogram
   */
  void ProcessNoCoverageRegion(const Interval& region) const;

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
   * Retrieve the `ConsensusAccuracyMetrics` object that contains the metrics calculated by this worker.
   */
  std::shared_ptr<Metrics> GetMetrics() const {
    return _metrics;
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
  std::optional<std::string> _ref_seq{};
  // The options for calculating metrics
  const AlignmentMetricsOptions& _options;
  // Metadata about the dataset being processed
  const DatasetMetadata& _dataset_metadata;
  // A lookup table to quickly whether a read has already been processed when processing an earlier super region
  const RegionLookupTable& _region_lookup_table;
  // The alignment reader used to read alignments from the input BAM/CRAM file
  const std::shared_ptr<AlignmentReader> _alignment_reader;
  // The metrics object that holds the calculated metrics
  const std::shared_ptr<Metrics> _metrics;
};

}  // namespace xoos::alignment_metrics
