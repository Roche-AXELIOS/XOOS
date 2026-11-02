#pragma once

#include <string_view>
#include <unordered_map>

#include <xoos/concurrent/enumerable-thread-local.h>
#include <xoos/histogram/histogram-summary.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "metrics/metrics-constraints.h"  // IWYU pragma: keep
#include "metrics/read-completeness.h"
#include "sequence/matcher/match-info.h"
#include "sequence/sequence-enums.h"

namespace xoos::demux {

struct DemuxAndTrimParam;

// ============================================================================
// Simplex-specific Run Metrics
// ============================================================================
constexpr std::string_view kPreassignmentPassingReads = "preassignment_passing_reads";
constexpr std::string_view kFailedReads = "failed_reads";
constexpr std::string_view kFullSimplexReads = "full_simplex_reads";
constexpr std::string_view kPartialSimplexReads = "partial_simplex_reads";
constexpr std::string_view kBothSidDetectedReads = "both_sid_detected_reads";
constexpr std::string_view kSidDiscordantReads = "sid_discordant_reads";
constexpr std::string_view kSidConcordantReads = "sid_concordant_reads";
constexpr std::string_view kUmiFailedReads = "umi_failed_reads";
constexpr std::string_view kPerfectIndexReads = "perfect_index_reads";
constexpr std::string_view kFullSimplexReadsPercentage = "full_simplex_reads_percentage";
constexpr std::string_view kPartialSimplexReadsPercentage = "partial_simplex_reads_percentage";
constexpr std::string_view kBothSidDetectedReadsPercentage = "both_sid_detected_reads_percentage";
constexpr std::string_view kSidDiscordantReadsPercentage = "sid_discordant_reads_percentage";
constexpr std::string_view kSidConcordantReadsPercentage = "sid_concordant_reads_percentage";
constexpr std::string_view kUmiFailedReadsPercentage = "umi_failed_reads_percentage";
constexpr std::string_view kPerfectIndexReadsPercentage = "perfect_index_reads_percentage";
// Simplex-only base metric: bases of pre-assignment filtered reads. Duplex has no `filtered_reads`
// category (those reads fold into `unassigned`), so this name lives with the simplex constants.
constexpr std::string_view kFilteredBases = "filtered_bases";

enum class SequenceFound {
  kYes = 0,
  kNo = 1,
};

using enum SequenceFound;

// Helper functions to handle optional members with compile-time detection
template <class T>
SequenceFound GetUmi5p(const T& record) {
  if constexpr (requires { record.umi_5p; }) {
    return record.umi_5p ? kYes : kNo;
  } else {
    return kNo;
  }
}

template <class T>
SequenceFound GetUmi3p(const T& record) {
  if constexpr (requires { record.umi_3p; }) {
    return record.umi_3p ? kYes : kNo;
  } else {
    return kNo;
  }
}

struct SequencesFound {
  SequenceFound sid_5p{kNo};
  SequenceFound umi_5p{kNo};
  SequenceFound umi_3p{kNo};
  SequenceFound sid_3p{kNo};

  bool operator==(const SequencesFound& other) const = default;
};

struct SequencesFoundHash {
  size_t operator()(const SequencesFound& sf) const;
};

struct AdapterMetric {
  std::optional<u32> sid;
  SequencesFound found;

  bool operator==(const AdapterMetric& other) const = default;
};

struct AdapterMetricHash {
  size_t operator()(const AdapterMetric& am) const;
};

struct SidCollision {
  ReadEnd end{};
  u32 expected_sid{};
  u32 collided_sid{};

  bool operator==(const SidCollision& other) const = default;
};

struct SidCollisionHash {
  size_t operator()(const SidCollision& sc) const;
};

/**
 * Responsible for tracking metrics about the trimming of adapters and assignment
 * of reads to SIDs.
 *
 * The Metrics class uses a concurrent::EnumerableThreadLocal to keep a count of metrics per
 * thread to avoid communication between threads. Before reporting the metrics they should be
 * aggregated with SumTotal()
 */
class SimplexMetrics {
 public:
  static SimplexMetrics& Instance();

  /**
   * Aggregate the thread specific metrics together into a total metrics count.
   */
  static SimplexMetrics SumTotal();

  SimplexMetrics();

  void IncrementInputReadCount();
  void RecordInputReadLength(u64 raw_length);
  void RecordRawFailedReadLength(u64 raw_length);
  void WriteMetrics(const DemuxAndTrimParam& param, const SidPool& sid_pool) const;

  /**
   * @brief Records per-read detection metrics for every assigned read. The passing-only read-length
   * histograms and scalar full/partial counters are recorded separately after post-assignment filters.
   *
   * @tparam T Any trim-info type with sid_5p/sid_3p fields and optionally umi_5p/umi_3p.
   * @param[in] record The trim-info record to extract metrics from.
   * @param[in] untrimmed_read_length Retained for compatibility with existing callers; no distribution is
   * recorded until RecordPassingReadLength().
   */
  template <class T>
  void RecordReadMetrics(const T& record, const std::optional<u64> /*untrimmed_read_length*/) {
    auto sf = SequencesFound{
        .sid_5p = record.sid_5p ? kYes : kNo,
        .umi_5p = GetUmi5p(record),
        .umi_3p = GetUmi3p(record),
        .sid_3p = record.sid_3p ? kYes : kNo,
    };

    IncrementAdapterCount(AdapterMetric{record.sid, sf}, u64{1});

    // Adapter counts above are recorded for every read; the remaining metrics are per-SID, so
    // no-SID reads have nothing further to record.
    if (!record.sid) {
      return;
    }

    const u32 sid = *record.sid;

    IncrementCount(&_sid_read_counts, sid, u64{1});
    if (IsPerfectMatch(record)) {
      IncrementCount(&_sid_perfect_index_read_counts, sid, u64{1});
    }

    HandleFullPartialCount(record, sid);
  }

  /**
   * @brief Records passing-only run-level and sample-level length distributions.
   *
   * The caller must invoke this only after all post-assignment filters have accepted the read.
   * `insert.epos` is the exclusive zero-based boundary before a 3' adapter. The end-adapter
   * distribution uses its one-based position, or zero when the read has no 3' adapter.
   */
  template <class T>
  void RecordPassingReadLength(const T& record, const u64 raw_length) {
    if (!record.sid.has_value()) {
      return;
    }
    const u32 sid = *record.sid;
    const auto trimmed_length = static_cast<u64>(record.insert.Length());
    const auto endadapter_position = record.insert.epos < raw_length ? static_cast<u64>(record.insert.epos) + 1 : 0;

    _untrimmed_passing_read_len_dist.AddCountToHistogram(raw_length, u64{1});
    _trimmed_passing_read_len_dist.AddCountToHistogram(trimmed_length, u64{1});
    _sample_untrimmed_simplex_passing_read_length_distr.at(sid).AddCountToHistogram(raw_length, u64{1});
    _sample_trimmed_simplex_passing_read_length_distr.at(sid).AddCountToHistogram(trimmed_length, u64{1});
    _total_endadapter_position_distr.AddCountToHistogram(endadapter_position, u64{1});
    _sample_endadapter_position_distr.at(sid).AddCountToHistogram(endadapter_position, u64{1});

    if (IsFullRead(record)) {
      _untrimmed_full_passing_read_len_dist.AddCountToHistogram(raw_length, u64{1});
      _trimmed_full_passing_read_len_dist.AddCountToHistogram(trimmed_length, u64{1});
      _sample_untrimmed_full_simplex_passing_read_length_distr.at(sid).AddCountToHistogram(raw_length, u64{1});
      _sample_trimmed_full_simplex_passing_read_length_distr.at(sid).AddCountToHistogram(trimmed_length, u64{1});
    } else {
      _untrimmed_partial_passing_read_len_dist.AddCountToHistogram(raw_length, u64{1});
      _trimmed_partial_passing_read_len_dist.AddCountToHistogram(trimmed_length, u64{1});
      _sample_untrimmed_partial_simplex_passing_read_length_distr.at(sid).AddCountToHistogram(raw_length, u64{1});
      _sample_trimmed_partial_simplex_passing_read_length_distr.at(sid).AddCountToHistogram(trimmed_length, u64{1});
    }
  }

  /**
   * @brief Records the full/partial read counters for a passing, SID-assigned read; the caller must
   * invoke this only for reads that pass all post-assignment gates.
   *
   * @tparam T Any trim-info type with sid_5p/sid_3p fields and optionally umi_5p/umi_3p.
   * @param[in] record The trim-info record of a passing, SID-assigned read.
   */
  template <class T>
  void RecordCompleteness(const T& record) {
    // Contract: only passing, SID-assigned reads reach here. Guard against misuse so an
    // unassigned record fails safely instead of dereferencing an empty optional.
    if (!record.sid.has_value()) {
      return;
    }
    const u32 sid = *record.sid;
    IncrementCount(IsFullRead(record) ? &_sid_full_read_counts : &_sid_partial_read_counts, sid, u64{1});
  }

  void WriteRunMetrics(const SidPool& sid_pool, const DemuxAndTrimParam& param) const;

  void WriteSampleMetrics(const SidPool& sid_pool, const DemuxAndTrimParam& param) const;

  void WriteReadLengthDistributions(const SidPool& sid_pool, const DemuxAndTrimParam& param) const;

  void WriteSampleAssignmentMetrics(const SidPool& sid_pool, const DemuxAndTrimParam& param) const;

  u64 TotalInputReadCount() const;

  u64 GetSidReadCount(u32 id) const;

  /// Per-SID `total_bases`: insert (post-trim) bases of that SID's passing reads.
  u64 GetSidTotalBaseCount(u32 id) const;
  /// Per-SID `trimmed_bases`: bases removed by trimming from that SID's passing reads.
  u64 GetSidTrimmedBaseCount(u32 id) const;

  u64 TotalAssignedReadCount(const SidPool& sid_pool) const;
  u64 TotalFoundSids(const SidPool& sid_pool) const;

 private:
  /**
   * @brief Records the both-SID-detected, concordant/discordant, and SID-collision counters for an
   * assigned read.
   *
   * @tparam T Any trim-info type with sid_5p/sid_3p fields.
   * @param[in] record The trim-info record of an assigned read.
   * @param[in] sid The resolved SID index the read is assigned to.
   */
  template <class T>
  void HandleFullPartialCount(const T& record, const u32 sid) {
    if (IsBothSidDetected(record)) {
      IncrementCount(&_sid_both_sid_detected_read_counts, sid, u64{1});
      // Concordant is the complement of discordant within both-SID-detected reads.
      if (IsIndexDiscordant(record)) {
        IncrementCount(&_sid_sid_discordant_read_counts, sid, u64{1});
      } else {
        IncrementCount(&_sid_sid_concordant_read_counts, sid, u64{1});
      }
    }

    if (record.sid_5p && sid != record.sid_5p) {
      ++_sid_collisions[SidCollision{ReadEnd::k5p, sid, *record.sid_5p}];
    }

    if (record.sid_3p && sid != record.sid_3p) {
      ++_sid_collisions[SidCollision{ReadEnd::k3p, sid, *record.sid_3p}];
    }
  }

  template <typename T>
  void IncrementCount(vec<T>* elements, size_t position, T amount) const {
    if (position >= elements->size()) {
      elements->resize(position + 1, 0);
    }
    elements->at(position) += amount;
  }

  template <typename T>
  T GetCount(const vec<T>& elements, size_t position) const {
    return position < elements.size() ? elements.at(position) : 0;
  }

  static u64 TotalCountForPool(const SidPool& sid_pool, const std::function<u64(u32)>& func);

  u64 TotalFilteredReadCount() const;

  u64 TotalUnassignedReadCount() const;

  u64 GetSidFullReadCount(u32 id) const;

  static u64 TotalSidCount(const SidPool& sid_pool);
  u64 TotalPerfectIndexReadCount(const SidPool& sid_pool) const;
  u64 GetSidPerfectIndexReadCount(u32 id) const;

  u64 GetSidPartialReadCount(u32 id) const;

  u64 GetSidPassingReadCount(u32 id) const;

  u64 TotalFullReadCount(const SidPool& sid_pool) const;

  u64 TotalPartialReadCount(const SidPool& sid_pool) const;

  u64 TotalBothSidDetectedReadCount(const SidPool& sid_pool) const;
  u64 GetSidBothSidDetectedReadCount(u32 id) const;

  u64 TotalSidDiscordantReadCount(const SidPool& sid_pool) const;
  u64 GetSidSidDiscordantReadCount(u32 id) const;

  u64 TotalSidConcordantReadCount(const SidPool& sid_pool) const;
  u64 GetSidSidConcordantReadCount(u32 id) const;

  u64 GetSidUmiFailedReadCount(u32 id) const;

  void IncrementAdapterCount(const AdapterMetric& am, u64 count);

  u64 GetAdapterCount(const AdapterMetric& am) const;

  u64 GetSidCollisionCount(ReadEnd end, u32 expected_sid, u32 collided_sid) const;

  void Add(const SimplexMetrics& other);

  // variables
  static constexpr auto kSidPoolSizeHint = 1000;
  static thread_local concurrent::EnumerableThreadLocal<SimplexMetrics> instance;
  static const vec<SequencesFound> kAllSequencesFound;

 public:
  /// Passing reads before assignment
  u64 preassign_passing_read_count{};

  /// The reads being written to the final FASTQ.
  u64 passing_read_count{};

  /// Increment the per-SID passing-read counter at the success point.
  void IncrementSidPassingReadCount(u32 id);

  /// Increment the per-SID umi_failed counter for a read assigned to `id` that failed the UMI gate.
  void IncrementSidUmiFailedReadCount(u32 id);

  /// Accumulate a passing read's insert (total) and trimmed base counts for a SID, driving both the
  /// run-level and per-SID base totals from a single point so their sums cannot diverge.
  void AddSidPassingBases(u32 id, u64 total_bases, u64 trimmed_bases);

  /// The total number of short reads from input file
  u64 too_short_read_count{};

  /// The total number of reads too long to fit in the fixed record buffer.
  u64 too_long_read_count{};

  /// The total number of reads filtered by the pre-assignment input-length checks (too short or too long).
  u64 filtered_read_count{};

  /// The total number of reads that passed the pre-filtering step but were not assigned a SID.
  u64 unassigned_read_count{};

  /// The total number of filtered trimmed reads from input file
  u64 too_short_trimmed_read_count{};

  /// The total number of reads discarded due to discordant SID mode
  u64 sid_discordant_discarded_read_count{};

  /// The total number of SID-assigned partial reads discarded in full-only mode.
  u64 full_only_discarded_read_count{};

  /// The total number of SID-assigned UMI-simplex (YSU) reads discarded for lacking a SID+UMI combination.
  u64 umi_failed_read_count{};

  /// Raw (pre-trim) bases summed over every input read (all `total_reads`). Percentage denominator for all base
  /// metrics.
  u64 raw_base_count{};

  /// Raw (pre-trim) bases of reads counted as `filtered_reads` (pre-assignment too-short / too-long).
  u64 filtered_base_count{};

  /// Raw (pre-trim) bases of reads counted as `unassigned_reads`.
  u64 unassigned_base_count{};

  /// Raw (pre-trim) bases of reads counted as `failed_assigned_reads`.
  u64 failed_assigned_base_count{};

  /// Trimmed (post-trim) bases of `passing_reads`; reported as `total_bases`.
  u64 total_base_count{};

  /// Bases removed by trimming from `passing_reads` (raw minus trimmed); reported as `trimmed_bases`.
  u64 trimmed_base_count{};

 private:
  using LengthHistogram = histogram::Histogram<u64>;
  using PerSampleHistogram = vec<LengthHistogram>;

  // Run-level histograms.
  LengthHistogram _total_read_length_distr{metrics_constraints::max_logged_read_length + 1, 0};
  LengthHistogram _raw_failed_read_length_distr{metrics_constraints::max_logged_read_length + 1, 0};
  LengthHistogram _untrimmed_passing_read_len_dist{metrics_constraints::max_logged_read_length + 1, 0};
  LengthHistogram _untrimmed_full_passing_read_len_dist{metrics_constraints::max_logged_read_length + 1, 0};
  LengthHistogram _untrimmed_partial_passing_read_len_dist{metrics_constraints::max_logged_read_length + 1, 0};
  LengthHistogram _trimmed_passing_read_len_dist{metrics_constraints::max_logged_read_length + 1, 0};
  LengthHistogram _trimmed_full_passing_read_len_dist{metrics_constraints::max_logged_read_length + 1, 0};
  LengthHistogram _trimmed_partial_passing_read_len_dist{metrics_constraints::max_logged_read_length + 1, 0};
  LengthHistogram _total_endadapter_position_distr{metrics_constraints::max_logged_read_length + 1, 0};

  // Per-sample histograms; index is the SID id in [0, max_sid_id_index].
  PerSampleHistogram _sample_untrimmed_simplex_passing_read_length_distr;
  PerSampleHistogram _sample_untrimmed_full_simplex_passing_read_length_distr;
  PerSampleHistogram _sample_untrimmed_partial_simplex_passing_read_length_distr;
  PerSampleHistogram _sample_trimmed_simplex_passing_read_length_distr;
  PerSampleHistogram _sample_trimmed_full_simplex_passing_read_length_distr;
  PerSampleHistogram _sample_trimmed_partial_simplex_passing_read_length_distr;
  PerSampleHistogram _sample_endadapter_position_distr;

  /// The total number of reads from input file
  u64 _total_input_read_count{};

  /// The number of reads assigned to each sid
  vec<u64> _sid_read_counts = vec<u64>(kSidPoolSizeHint);

  /// The number of reads with both 5p sid and 3p sid detected
  vec<u64> _sid_both_sid_detected_read_counts = vec<u64>(kSidPoolSizeHint);

  /// The number of full reads assigned to each sid, full means both umi detected and one sid detected
  vec<u64> _sid_full_read_counts = vec<u64>(kSidPoolSizeHint);

  /// The number of partial reads assigned to each sid, partial means one umi detected and one sid detected
  vec<u64> _sid_partial_read_counts = vec<u64>(kSidPoolSizeHint);

  /// The number of reads assigned to each sid that pass all filters.
  vec<u64> _sid_passing_read_counts = vec<u64>(kSidPoolSizeHint);

  /// Per-SID insert (post-trim) bases of passing reads; sums to the run-level `total_base_count`.
  vec<u64> _sid_total_base_counts = vec<u64>(kSidPoolSizeHint);

  /// Per-SID bases removed by trimming from passing reads; sums to the run-level `trimmed_base_count`.
  vec<u64> _sid_trimmed_base_counts = vec<u64>(kSidPoolSizeHint);

  /// The number of reads with both sid detected, zero mismatch in both sid, and the same sid
  vec<u64> _sid_perfect_index_read_counts = vec<u64>(kSidPoolSizeHint);

  /// The number of reads with both sid detected, sid are not the same
  vec<u64> _sid_sid_discordant_read_counts = vec<u64>(kSidPoolSizeHint);

  /// The number of reads with both sid detected and the same sid (concordant)
  vec<u64> _sid_sid_concordant_read_counts = vec<u64>(kSidPoolSizeHint);

  /// The number of reads assigned to each sid that failed the UMI gate (umi_failed)
  vec<u64> _sid_umi_failed_read_counts = vec<u64>(kSidPoolSizeHint);

  /**
   * Track the number of adapters of each sample which followed a certain pattern of:
   * (5' sid found, 5' umi found, 3' umi found, 3' sid found)
   */
  std::unordered_map<AdapterMetric, u64, AdapterMetricHash> _adapter_counts;

  std::unordered_map<SidCollision, u64, SidCollisionHash> _sid_collisions;
};

}  // namespace xoos::demux
