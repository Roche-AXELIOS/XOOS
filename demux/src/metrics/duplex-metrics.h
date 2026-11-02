#pragma once

#include <string_view>

#include <xoos/concurrent/enumerable-thread-local.h>
#include <xoos/histogram/histogram-summary.h>
#include <xoos/types/vec.h>

#include "adapter-design/adapter-design.h"
#include "core/demux-and-trim-pipeline.h"
#include "metrics-constraints.h"
#include "sequence/matcher/match-info.h"

namespace xoos::demux {

// +1 to include 0 length reads
using LengthHistogram = histogram::Histogram<u64>;
using PerSIDHistogram = vec<LengthHistogram>;
// counter for every SID
using PerSIDCount = vec<u64>;

// ============================================================================
// Duplex-specific Run Metrics
// ============================================================================
constexpr std::string_view kConcordantDuplexBases = "concordant_duplex_bases";
constexpr std::string_view kDiscordantDuplexBases = "discordant_duplex_bases";
constexpr std::string_view kAmbiguousDuplexBases = "ambiguous_duplex_bases";
constexpr std::string_view kDuplexBases = "duplex_bases";
constexpr std::string_view kFullDuplexReads = "full_duplex_reads";
constexpr std::string_view kPartialDuplexReads = "partial_duplex_reads";
constexpr std::string_view kNoHairpinReads = "no_hairpin_reads";
constexpr std::string_view kTooLongConsensusReads = "too_long_consensus_reads";
constexpr std::string_view kTooManyErrorsReads = "too_many_errors_reads";
constexpr std::string_view kNoEndadapterReads = "no_endadapter_reads";
constexpr std::string_view kFailedHairpinStemTrimReads = "failed_hairpin_stem_trim_reads";
constexpr std::string_view kFoundByStringCompare = "found_by_string_compare";
constexpr std::string_view kFoundByGlobalSymmetry = "found_by_global_symmetry";
constexpr std::string_view kFoundByLocalSymmetry = "found_by_local_symmetry";
constexpr std::string_view kLongerR2Reads = "longer_r2_reads";
constexpr std::string_view kLongerR2FullDuplexReads = "longer_r2_full_duplex_reads";
constexpr std::string_view kMultiCandidateAlignments = "multi_candidate_alignments";
constexpr std::string_view kBothUmiReads = "both_umi_reads";
constexpr std::string_view k5pUmiReads = "5p_umi_reads";
constexpr std::string_view k3pUmiReads = "3p_umi_reads";
constexpr std::string_view kStrandFwReads = "strand_fw_reads";
constexpr std::string_view kStrandRvReads = "strand_rv_reads";
constexpr std::string_view kStrandFwSigReads = "strand_fw_sig_reads";
constexpr std::string_view kStrandRvSigReads = "strand_rv_sig_reads";
constexpr std::string_view kNonDuplexBases = "non_duplex_bases";
constexpr std::string_view kFullDuplexReadsPercentage = "full_duplex_reads_percentage";
constexpr std::string_view kPartialDuplexReadsPercentage = "partial_duplex_reads_percentage";
constexpr std::string_view kFailedAssignedReadsPercentage = "failed_assigned_reads_percentage";
constexpr std::string_view kBothUmiReadsPercentage = "both_umi_reads_percentage";
constexpr std::string_view k5pUmiReadsPercentage = "5p_umi_reads_percentage";
constexpr std::string_view k3pUmiReadsPercentage = "3p_umi_reads_percentage";
constexpr std::string_view kNoEndadapterReadsPercentage = "no_endadapter_reads_percentage";
constexpr std::string_view kConcordantDuplexBasesPercentage = "concordant_duplex_bases_percentage";
constexpr std::string_view kDiscordantDuplexBasesPercentage = "discordant_duplex_bases_percentage";
constexpr std::string_view kAmbiguousDuplexBasesPercentage = "ambiguous_duplex_bases_percentage";
constexpr std::string_view kDuplexBasesPercentage = "duplex_bases_percentage";
constexpr std::string_view kNonDuplexBasesPercentage = "non_duplex_bases_percentage";

/**
 * Responsible for tracking metrics about the trimming of adapters and assignment
 * of reads to SIDs for Duplex HD data.
 *
 * The Metrics class uses a concurrent::EnumerableThreadLocal to keep a count of metrics per
 * thread to avoid communication between threads. Before reporting the metrics they should be
 * aggregated with SumTotal()
 *
 * To minimize overhead and complexity, I am using a struct rather than a class.
 */
struct DuplexMetrics {
  static thread_local concurrent::EnumerableThreadLocal<DuplexMetrics> instance;

  static DuplexMetrics& Instance();

  /**
   * Aggregate the thread specific metrics together into a total metrics count.
   */
  static DuplexMetrics SumTotal();

  DuplexMetrics();

  void WriteMetrics(const DemuxAndTrimParam& params, const SidPool& sid_pool, AdapterType adapter_type) const;

  void WriteRunMetrics(const DemuxAndTrimParam& params, const SidPool& sid_pool, AdapterType adapter_type) const;

  void WriteSampleMetrics(const DemuxAndTrimParam& params, const SidPool& sid_pool, AdapterType adapter_type) const;

  void WriteReadLengthDistributions(const DemuxAndTrimParam& params, const SidPool& sid_pool) const;

  /**
   * @brief Sum the passing consensus bases written to output for a single SID.
   *
   * concordant, ambiguous, discordant and simplex are disjoint partitions (bases downgraded to the
   * ambiguous quality score are removed from the concordant count), so all four are summed. This must
   * match the run-level total and the duplex_bases + non_duplex_bases partition.
   *
   * @param id SID index into the per-SID passing count vectors.
   * @return Total passing bases for the SID, or 0 if @p id is out of range.
   */
  u64 SumPassingTotalBases(u32 id) const;

  /**
   * @brief Count the number of SIDs that have at least one assigned read (passing or failed).
   *
   * @param sid_pool The pool of SIDs to check against the collected metrics.
   * @return Number of SIDs with a non-zero assigned read count.
   */
  u64 CountAssignedSids(const SidPool& sid_pool) const;

  void Add(const DuplexMetrics& other);

  static size_t MinConcordDupBases();

  /**
   * @brief Counts on why a read failed before it was assigned to a SID.
   * @note Should sum to the total number of unassigned reads.
   */
  struct UnassignedCounts {
    // Also called a non-hairpin, no hairpin was found and thus not assigned to a SID
    u64 no_hairpin_found{0};
    // Reads that were too long so we couldn't fit them in the buffer
    u64 read_too_long{0};
    // Reads that were too short according to the min read length filter
    u64 read_too_short{0};

    // Base count metrics
    // Bases counts of unassigned reads
    u64 raw_bases{0};
  };

  /**
   * @brief Counts counting properties of reads that passed the filter.
   * @note full_duplex and partial_duplex should sum to the total number of passing reads.
   */
  struct PassingCounts {
    // Full-length reads (if consensus alignment is equal to or less than the trimmed read length)
    PerSIDCount full_duplex{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
    // Partial duplex reads (if consensus alignment is less than the trimmed read length + 1)
    PerSIDCount partial_duplex{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
    // R1 is swapped with R2 before alignment because R2 was longer than R1
    PerSIDCount longer_r2{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
    // R1 is swapped with R2 before alignment because R2 was longer than R1 and the read was full length
    PerSIDCount longer_r2_full_duplex{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
    // Reads that have both UMIs
    PerSIDCount both_umi{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
    // Reads that have only 5p UMI
    PerSIDCount only_5p_umi{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
    // Reads that have only 3p UMI
    PerSIDCount only_3p_umi{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
    // Reads that didn't have an endadapter at all
    PerSIDCount no_endadapter{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
    // Reads where multiple tied alignment candidates required gap-based selection
    PerSIDCount multi_candidate_alignments{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};

    // Base count metrics
    // Base counts of passing reads before trimming
    PerSIDCount raw_bases{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
    // Concordant duplex bases
    PerSIDCount concordant_bases{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
    // Discordant bases
    PerSIDCount discordant_bases{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
    // Duplex bases downgraded to topologically ambiguous positions.
    PerSIDCount ambiguous_bases{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
    // i.e. Non-duplex bases
    PerSIDCount simplex_bases{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
  };

  /**
   * @brief Per sample counts on why a read failed even if it was assigned to a SID.
   * @note The sum of this plus passing reads should equal the total number of assigned reads.
   */
  struct FailedAssignedCounts {
    // trimmed read that were too short after processing (e.g. after endadapter trimming)
    PerSIDCount trimmed_read_too_short{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
    // consensus sequence didn't fit in consensus read buffer
    PerSIDCount consensus_too_long{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
    // consensus sequence had too many error according to the edit distance threshold
    PerSIDCount too_many_errors{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
    // Reads that had extra expected loop sequence but we failed to trim it
    PerSIDCount failed_hairpin_stem_trim_reads{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};

    // Base count metrics
    // Base counts for failed assigned reads
    PerSIDCount raw_bases{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
  };

  /**
   * @brief Counts on how we found the mid adapter.
   * @note Should sum to the total number of assigned reads.
   */
  struct MidAdapterCounts {
    PerSIDCount found_by_string_compare{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
    PerSIDCount found_by_global_symmetry{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
    PerSIDCount found_by_local_symmetry{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
  };

  /**
   * @brief Counts on strand detection outcomes. Only useful when strand detection is enabled.
   * @note Should only be used for reads that passed the filter but will not sum to total of passing reads because some
   * reads may be ambiguous.
   *
   */
  struct StrandCounts {
    // Number of forward reads (not significant)
    PerSIDCount fw{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
    // Number of reversed reads (not significant)
    PerSIDCount rv{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
    // Number of significant forward reads
    PerSIDCount fw_sig{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
    // Number of significant reversed reads
    PerSIDCount rv_sig{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
  };

  // Structs for count categories for organizing the metrics
  UnassignedCounts unassigned_counts;
  FailedAssignedCounts failed_assigned_counts;
  MidAdapterCounts midadapter_counts;
  StrandCounts strand_counts;
  PassingCounts passing_counts;

  // Distribution metrics (unknown SIDs)
  // TODO: Putting these in relevant structs above could remove redundancy (counts can be derived from distributions)
  // Raw read length distribution before any trimming
  LengthHistogram total_length_distr{metrics_constraints::max_logged_read_length + 1, 0};
  // Distribution of reads that were unassigned either due to no midadapter found, too short, or too long
  LengthHistogram unassigned_length_distr{metrics_constraints::max_logged_read_length + 1, 0};
  // Distribution of reads unassigned due to missing midadapter / hairpin
  LengthHistogram no_hairpin_length_distr{metrics_constraints::max_logged_read_length + 1, 0};
  LengthHistogram raw_failed_length_distr{metrics_constraints::max_logged_read_length + 1, 0};
  LengthHistogram untrimmed_passing_length_distr{metrics_constraints::max_logged_read_length + 1, 0};
  LengthHistogram untrimmed_full_passing_length_distr{metrics_constraints::max_logged_read_length + 1, 0};
  LengthHistogram untrimmed_partial_passing_length_distr{metrics_constraints::max_logged_read_length + 1, 0};
  LengthHistogram trimmed_passing_length_distr{metrics_constraints::max_logged_read_length + 1, 0};
  LengthHistogram trimmed_full_passing_length_distr{metrics_constraints::max_logged_read_length + 1, 0};
  LengthHistogram trimmed_partial_passing_length_distr{metrics_constraints::max_logged_read_length + 1, 0};
  LengthHistogram total_endadapter_position_distr{metrics_constraints::max_logged_read_length + 1, 0};

  void RecordRawFailedReadLength(u64 raw_length) {
    raw_failed_length_distr.AddCountToHistogram(raw_length, 1);
  }

  // Distribution metrics (known SIDs)
  PerSIDHistogram sample_trimmed_full_duplex_passing_read_length_distr{
      metrics_constraints::max_sid_id_index + 1, LengthHistogram(metrics_constraints::max_logged_read_length + 1, 0)};
  PerSIDHistogram sample_trimmed_partial_duplex_passing_read_length_distr{
      metrics_constraints::max_sid_id_index + 1, LengthHistogram(metrics_constraints::max_logged_read_length + 1, 0)};
  // Distribution of where the furthest right position of the endadapter was found in the read before trimming
  PerSIDHistogram sample_endadapter_position_distr{metrics_constraints::max_sid_id_index + 1,
                                                   LengthHistogram(metrics_constraints::max_logged_read_length + 1, 0)};
  // Distribution of the number of bases in the consensus sequence portion of the read
  PerSIDHistogram sample_trimmed_duplex_passing_read_length_distr{
      metrics_constraints::max_sid_id_index + 1, LengthHistogram(metrics_constraints::max_logged_read_length + 1, 0)};
  PerSIDHistogram sample_untrimmed_duplex_passing_read_length_distr{
      metrics_constraints::max_sid_id_index + 1, LengthHistogram(metrics_constraints::max_logged_read_length + 1, 0)};
  PerSIDHistogram sample_untrimmed_full_duplex_passing_read_length_distr{
      metrics_constraints::max_sid_id_index + 1, LengthHistogram(metrics_constraints::max_logged_read_length + 1, 0)};
  PerSIDHistogram sample_untrimmed_partial_duplex_passing_read_length_distr{
      metrics_constraints::max_sid_id_index + 1, LengthHistogram(metrics_constraints::max_logged_read_length + 1, 0)};
};

// Helper function to report strand metrics
void ReportStrandMetrics(const DuplexMetrics& global_results);

}  // namespace xoos::demux
