#pragma once

#include <xoos/concurrent/enumerable-thread-local.h>
#include <xoos/histogram/histogram-summary.h>
#include <xoos/io/metadata-util.h>

#include <unordered_map>
#include <vector>

#include "core/demux-and-trim-pipeline.h"
#include "metrics-constraints.h"
#include "sequence/matcher/match-info.h"
#include "xoos/types/str-container.h"

namespace fs = std::filesystem;

namespace xoos::demux {

// +1 to include 0 length reads
using LengthHistogram = histogram::Histogram<u64>;
using PerSIDHistogram = vec<LengthHistogram>;
// counter for every SID
using PerSIDCount = std::vector<u64>;

// ============================================================================
// Duplex-specific Run Metrics
// ============================================================================
constexpr auto kConcordantDuplexBases = "concordant_duplex_bases";
constexpr auto kDiscordantDuplexBases = "discordant_duplex_bases";
constexpr auto kDuplexBases = "duplex_bases";
constexpr auto kFullDuplexReads = "full_duplex_reads";
constexpr auto kPartialDuplexReads = "partial_duplex_reads";
constexpr auto kNoHairpinReads = "no_hairpin_reads";
constexpr auto kTooLongReads = "too_long_reads";
constexpr auto kTooLongConsensusReads = "too_long_consensus_reads";
constexpr auto kNoEndadapterReads = "no_endadapter_reads";
constexpr auto kFailedHairpinStemTrimReads = "failed_hairpin_stem_trim_reads";
constexpr auto kFoundByStringCompare = "found_by_string_compare";
constexpr auto kFoundByGlobalSymmetry = "found_by_global_symmetry";
constexpr auto kFoundByLocalSymmetry = "found_by_local_symmetry";
constexpr auto kLongerR2Reads = "longer_r2_reads";
constexpr auto kLongerR2FullDuplexReads = "longer_r2_full_duplex_reads";
constexpr auto kBothUmiReads = "both_umi_reads";
constexpr auto k5pUmiReads = "5p_umi_reads";
constexpr auto k3pUmiReads = "3p_umi_reads";
constexpr auto kStrandFwReads = "strand_fw_reads";
constexpr auto kStrandRvReads = "strand_rv_reads";
constexpr auto kStrandFwSigReads = "strand_fw_sig_reads";
constexpr auto kStrandRvSigReads = "strand_rv_sig_reads";
constexpr auto kTotalBases = "total_bases";
constexpr auto kNonDuplexBases = "non_duplex_bases";
constexpr auto kMeanPassingReadLength = "mean_passing_read_length";

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

  using SidPool = const std::unordered_map<uint, Barcode>;

  DuplexMetrics();

  void WriteMetrics(const DemuxAndTrimParam& params, const SidPool& sid_pool) const;

  void WriteRunMetrics(const DemuxAndTrimParam& params, const SidPool& sid_pool,
                       const StrMap<u32>& sid_name_to_id) const;

  void WriteSampleMetrics(const DemuxAndTrimParam& params, const SidPool& sid_pool,
                          const StrMap<u32>& sid_name_to_id) const;

  void WriteReadLengthDistributions(const DemuxAndTrimParam& params) const;

  void WriteSampleAssignmentMetrics(const DemuxAndTrimParam& params, const SidPool& sid_pool) const;

  std::unordered_set<u32> FindAssignedSids(const StrMap<u32>& sid_name_to_id) const;

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

  /**
   * @brief Metrics on the base counts in passing reads.
   * @note Every member should count per sample. The sum of all members should equal the total number of passing bases.
   */
  struct BaseCounts {
    // Concordant duplex bases
    PerSIDCount concordant{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
    // Discordant bases
    PerSIDCount discordant{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
    // i.e. Non-duplex bases
    PerSIDCount simplex{PerSIDCount(metrics_constraints::max_sid_id_index + 1, 0)};
  };

  // Structs for count categories for organizing the metrics
  UnassignedCounts unassigned_counts;
  FailedAssignedCounts failed_assigned_counts;
  MidAdapterCounts midadapter_counts;
  StrandCounts strand_counts;
  PassingCounts passing_counts;
  BaseCounts base_counts;

  // Distribution metrics (unknown SIDs)
  // TODO: Putting these in relevant structs above could remove redundancy (counts can be derived from distributions)
  // Raw read length distribution before any trimming
  LengthHistogram total_length_distr;
  // Distribution of reads that were unassigned either due to no midadapter found, too short, or too long
  LengthHistogram unassigned_length_distr;
  // Distribution of reads unassigned due to missing midadapter / hairpin
  LengthHistogram no_hairpin_length_distr;

  // Distribution metrics (known SIDs)
  PerSIDHistogram full_duplex_length_distr;
  PerSIDHistogram partial_duplex_length_distr;
  // Distribution of where the furthest right position of the endadapter was found in the read before trimming
  PerSIDHistogram endadapter_position_distr;
  // Distribution of the number of bases in the consensus sequence portion of the read
  PerSIDHistogram passing_length_distr;
};

// Helper function to report strand metrics
void ReportStrandMetrics(const DuplexMetrics& global_results);

// Helper function to convert from SID pool to a map of SID name to SID id for easier lookup when writing metrics
StrMap<u32> ConvertSidPoolToSidNameToIdMap(const detail::SidPool& sid_pool);

}  // namespace xoos::demux
