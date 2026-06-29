#pragma once

#include <xoos/concurrent/enumerable-thread-local.h>
#include <xoos/histogram/histogram-summary.h>
#include <xoos/io/metadata-util.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/util/hash.h>

#include <string>
#include <unordered_map>

#include "metrics/metrics-constraints.h"  // IWYU pragma: keep
#include "sequence/matcher/match-info.h"
#include "sequence/matcher/seq-matcher.h"

namespace fs = std::filesystem;

namespace xoos::demux {

struct DemuxAndTrimParam;

// ============================================================================
// Simplex-specific Run Metrics
// ============================================================================
constexpr auto kPreassignmentPassingReads = "preassignment_passing_reads";
constexpr auto kFailedReads = "failed_reads";
constexpr auto kFullReads = "full_reads";
constexpr auto kPartialReads = "partial_reads";
constexpr auto kBothSidDetectedReads = "both_sid_detected_reads";
constexpr auto kSidDiscordantReads = "sid_discordant_reads";
constexpr auto kIndexHoppingReads = "index_hopping_reads";
constexpr auto kPerfectIndexReads = "perfect_index_reads";

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

  using SidPool = const std::unordered_map<u32, Barcode>;

  SimplexMetrics();

  void IncrementInputReadCount();
  void IncrementMinTrimmedReadLenFilteredCount();
  void WriteMetrics(const DemuxAndTrimParam& param, const SidPool& sid_pool) const;

  template <class T>
  void AddTrimmedRead(const T& record, const std::optional<u64> untrimmed_read_length) {
    _total_read_count += 1;
    auto sf = SequencesFound{
        .sid_5p = record.sid_5p ? kYes : kNo,
        .umi_5p = GetUmi5p(record),
        .umi_3p = GetUmi3p(record),
        .sid_3p = record.sid_3p ? kYes : kNo,
    };

    IncrementAdapterCount(AdapterMetric{record.sid, sf}, u64{1});

    if (record.sid) {
      u32 sid = *record.sid;

      IncrementCount(&_sid_read_counts, sid, u64{1});

      if (IsFullRead(record)) {
        IncrementCount(&_sid_full_read_counts, sid, u64{1});

        if (untrimmed_read_length) {
          _untrimmed_full_read_len_dist.at(sid).AddCountToHistogram(*untrimmed_read_length, u64{1});
        }
        _trimmed_full_read_len_dist.at(sid).AddCountToHistogram(record.insert.Length(), u64{1});
      }

      if (IsPartialRead(record)) {
        IncrementCount(&_sid_partial_read_counts, sid, u64{1});

        if (untrimmed_read_length) {
          _untrimmed_partial_read_len_dist.at(sid).AddCountToHistogram(*untrimmed_read_length, u64{1});
        }
        _trimmed_partial_read_len_dist.at(sid).AddCountToHistogram(record.insert.Length(), u64{1});
      }

      if (IsPerfectMatch(record)) {
        IncrementCount(&_sid_perfect_index_read_counts, sid, u64{1});
      }

      if (IsBothSidDetected(record)) {
        IncrementCount(&_sid_both_sid_detected_read_counts, sid, u64{1});
      }

      if (IsIndexDiscordant(record)) {
        IncrementCount(&_sid_sid_discordant_read_counts, sid, u64{1});
      }

      if (record.sid_5p && sid != record.sid_5p) {
        _sid_collisions[SidCollision{ReadEnd::k5p, sid, *record.sid_5p}] += 1;
      }

      if (record.sid_3p && sid != record.sid_3p) {
        _sid_collisions[SidCollision{ReadEnd::k3p, sid, *record.sid_3p}] += 1;
      }
    }
  }

  void WriteRunMetrics(const SidPool& sid_pool, const DemuxAndTrimParam& param) const;

  void WriteSampleMetrics(const SidPool& sid_pool, const DemuxAndTrimParam& param) const;

  void WriteReadLengthDistributions(const SidPool& sid_pool, const DemuxAndTrimParam& param) const;

  void WriteSampleAssignmentMetrics(const SidPool& sid_pool, const DemuxAndTrimParam& param) const;

  u64 TotalReadCount() const;

  u64 TotalInputReadCount() const;

  u64 GetSidReadCount(u32 id) const;

  u64 TotalAssignedReadCount(const SidPool& sid_pool) const;
  u64 TotalFoundSids(const SidPool& sid_pool) const;

 private:
  template <class T>
  bool IsFullRead(const T& record) const {
    if constexpr (requires { record.umi_5p && record.umi_3p; }) {
      return (record.sid_5p || record.sid_3p) && (record.umi_5p && record.umi_3p);
    } else {
      // The adapter doesn't have UMI information
      // We define a full read as having both SIDs here
      return (record.sid_5p && record.sid_3p);
    }
  }

  template <class T>
  bool IsPartialRead(const T& record) const {
    if constexpr (requires { record.umi_5p && record.umi_3p; }) {
      return !IsFullRead(record) && (record.sid_5p || record.sid_3p) && (record.umi_5p || record.umi_3p);
    } else {
      // The adapter doesn't have UMI information
      // We define a partial read as having at least one SID here
      return !IsFullRead(record) && (record.sid_5p || record.sid_3p);
    }
  }

  template <class T>
  bool IsPerfectMatch(const T& record) const {
    return record.sid_5p_edist && record.sid_5p_edist == 0 && record.sid_3p_edist && record.sid_3p_edist == 0 &&
           record.sid_5p && record.sid_3p && record.sid_5p == record.sid_3p;
  }

  template <class T>
  bool IsBothSidDetected(const T& record) const {
    return record.sid_5p && record.sid_3p;
  }

  template <class T>
  bool IsIndexDiscordant(const T& record) const {
    return record.sid_5p && record.sid_3p && record.sid_5p != record.sid_3p;
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

  static f32 CalculatePercentage(u64 numerator, u64 denominator);
  static std::string FormatPercent(f32 data, s32 precision = 2);

  f32 PercentAssignedReads(const SidPool& sid_pool) const;

  u64 TotalFilteredReadCount() const;
  f32 PercentFilteredReads() const;

  u64 TotalUnassignedReadCount(const SidPool& sid_pool) const;

  u64 GetSidFullReadCount(u32 id) const;
  f32 PercentFullReads(const SidPool& sid_pool) const;

  u64 TotalIndexHoppingReadCount(const SidPool& sid_pool) const;
  u64 GetSidIndexHoppingReadCount(const SidPool& sid_pool, u32 id) const;
  f32 PercentIndexHoppingReads(const SidPool& sid_pool) const;

  static u64 TotalSidCount(const SidPool& sid_pool);
  u64 TotalPerfectIndexReadCount(const SidPool& sid_pool) const;
  u64 GetSidPerfectIndexReadCount(u32 id) const;
  f32 PercentPerfectIndexReads(const SidPool& sid_pool) const;

  u64 GetSidPartialReadCount(u32 id) const;
  f32 PercentPartialReads(const SidPool& sid_pool) const;

  u64 TotalFullReadCount(const SidPool& sid_pool) const;

  u64 TotalPartialReadCount(const SidPool& sid_pool) const;

  u64 TotalBothSidDetectedReadCount(const SidPool& sid_pool) const;
  u64 GetSidBothSidDetectedReadCount(u32 id) const;
  f32 PercentBothSidDetectedReads(const SidPool& sid_pool) const;

  u64 TotalSidDiscordantReadCount(const SidPool& sid_pool) const;
  u64 GetSidSidDiscordantReadCount(u32 id) const;
  f32 PercentSidDiscordantReads(const SidPool& sid_pool) const;

  void IncrementAdapterCount(const AdapterMetric& am, u64 count);

  u64 GetAdapterCount(const AdapterMetric& am) const;

  u64 GetSidCollisionCount(ReadEnd end, u32 expected_sid, u32 collided_sid) const;

  void Add(const SimplexMetrics& other);

  // variables
  static constexpr auto kSidPoolSizeHint = 1000;
  static thread_local concurrent::EnumerableThreadLocal<SimplexMetrics> instance;
  static const vec<SequencesFound> kAllSequencesFound;

  /// The total number of reads processed
  u64 _total_read_count{};

  // Per-SID read length histograms; index is SID id in [0, max_sid_id_index]
  // Histograms are sized to max_logged_read_length + 1 bins; outliers store lengths > max_logged_read_length
  vec<histogram::Histogram<u64>> _untrimmed_full_read_len_dist{};
  vec<histogram::Histogram<u64>> _untrimmed_partial_read_len_dist{};
  vec<histogram::Histogram<u64>> _trimmed_full_read_len_dist{};
  vec<histogram::Histogram<u64>> _trimmed_partial_read_len_dist{};

  /// The total number of reads from input file
  u64 _total_input_read_count{};

  /// The total number of filtered trimmed reads from input file
  u64 _min_trimmed_read_len_count{};

  /// The number of reads assigned to each sid
  vec<u64> _sid_read_counts = vec<u64>(kSidPoolSizeHint);

  /// The number of reads with both 5p sid and 3p sid detected
  vec<u64> _sid_both_sid_detected_read_counts = vec<u64>(kSidPoolSizeHint);

  /// The number of full reads assigned to each sid, full means both umi detected and one sid detected
  vec<u64> _sid_full_read_counts = vec<u64>(kSidPoolSizeHint);

  /// The number of partial reads assigned to each sid, partial means one umi detected and one sid detected
  vec<u64> _sid_partial_read_counts = vec<u64>(kSidPoolSizeHint);

  /// The number of reads with both sid detected, zero mismatch in both sid, and the same sid
  vec<u64> _sid_perfect_index_read_counts = vec<u64>(kSidPoolSizeHint);

  /// The number of reads with both sid detected, sid are not the same
  vec<u64> _sid_sid_discordant_read_counts = vec<u64>(kSidPoolSizeHint);

  /**
   * Track the number of adapters of each sample which followed a certain pattern of:
   * (5' sid found, 5' umi found, 3' umi found, 3' sid found)
   */
  std::unordered_map<AdapterMetric, u64, AdapterMetricHash> _adapter_counts{};

  std::unordered_map<SidCollision, u64, SidCollisionHash> _sid_collisions{};
};

}  // namespace xoos::demux
