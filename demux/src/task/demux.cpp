#include "task/demux.h"

#include <utility>

#include <xoos/log/logging.h>

#include "adapters/duplex/trim-info-duplex.h"
#include "metrics/duplex-metrics.h"
#include "metrics/read-completeness.h"
#include "metrics/simplex-metrics.h"
#include "task/flow-context.h"
#include "utility/stop-watch.h"

namespace xoos::demux {

namespace {

/// Records a post-assignment simplex failure and prevents the read from being written.
void HandleSimplexPostAssignmentFailure(FixedReadRecord& read,
                                        SimplexMetrics& metrics,
                                        const FixedReadRecord::Status status,
                                        const u32 sid) {
  using enum FixedReadRecord::Status;

  switch (status) {
    case kUmiFail:
      metrics.IncrementSidUmiFailedReadCount(sid);
      break;
    case kFullOnlyFail:
      ++metrics.full_only_discarded_read_count;
      break;
    case kDiscordantSidFail:
      ++metrics.sid_discordant_discarded_read_count;
      break;
    case kTrimmedTooShortFail:
      ++metrics.too_short_trimmed_read_count;
      break;
    default:
      throw error::Error("Invalid post-assignment simplex failure status {}", std::to_underlying(status));
  }

  metrics.failed_assigned_base_count += read.SeqLen();
  metrics.RecordRawFailedReadLength(read.SeqLen());
  read.SetStatus(status);
}

}  // namespace

Demux::Demux(FlowContext& exec, size_t batch_nr) : Task(exec, batch_nr, fmt::format("Demux{}", batch_nr)) {
}

/**
 * @brief Apply the post-assignment filters to a SID-assigned read, updating status and metrics.
 * @return true if the read was filtered out (caller should stop processing it), false otherwise.
 **/
template <typename TrimInfo>
bool Demux::SimplexPostAssignmentFilters(FixedReadRecord& read,
                                         SimplexMetrics& metrics,
                                         const FlowManager& mgr,
                                         const TrimInfo& trim_info) {
  // Failed-assigned reads contribute their raw (pre-trim) bases to failed_assigned_bases. The read
  // is not mutated before this point, so read.SeqLen() is still the raw length.

  // UMI gate (UMI adapters only): a SID-assigned read lacking a SID+UMI combination on any end is a
  // umi_failed read. Runs first so it takes precedence over later filters (e.g. discordant).
  if constexpr (requires { trim_info.umi_5p; }) {
    if (!HasSidUmiCombination(trim_info)) {
      // trim_info.sid is guaranteed set here: the caller only invokes the filters for assigned reads.
      HandleSimplexPostAssignmentFailure(read, metrics, FixedReadRecord::Status::kUmiFail, trim_info.sid.value());
      return true;
    }
  }
  // In full-only mode a partial read is still assigned a SID (so it is recorded in assigned/partial
  // metrics) but must not be written. A partial has a single SID and so can never be discordant,
  // so this check precedes the discordant gate.
  if (mgr.Param().read_length_mode == ReadLengthMode::kFullOnly && IsPartialRead(trim_info)) {
    HandleSimplexPostAssignmentFailure(read, metrics, FixedReadRecord::Status::kFullOnlyFail, trim_info.sid.value());
    return true;
  }
  // Discordance is checked before the too-short gate so a read that is both is counted as discordant.
  if (ShouldDiscardDiscordant(trim_info, mgr.Param().discordant_sid_mode)) {
    HandleSimplexPostAssignmentFailure(
        read, metrics, FixedReadRecord::Status::kDiscordantSidFail, trim_info.sid.value());
    return true;
  }
  // Assigned but too short after trimming: a failed-assigned read.
  if (trim_info.insert.Length() < mgr.Param().min_trimmed_read_len) {
    HandleSimplexPostAssignmentFailure(
        read, metrics, FixedReadRecord::Status::kTrimmedTooShortFail, trim_info.sid.value());
    return true;
  }
  return false;
}

/**
 * @brief RunSimplexAdapter is a helper function to run the demuxing and trimming for simplex adapters
 * @param read the read to be demultiplexed and trimmed
 * @param metrics the Metrics instance to update the metrics
 * @param raw_length
 * @param mgr the FlowManager instance to access the adapter and SID pool information
 * @param trim_info the TrimInfo object containing the results of the demuxing and trimming for the read
 **/
template <typename TrimInfo>
void Demux::RunSimplexAdapter(
    FixedReadRecord& read, SimplexMetrics& metrics, u32 raw_length, const FlowManager& mgr, const TrimInfo& trim_info) {
  // Sanity check: a SID from the LUT must be in range for the SID pool.
  if (trim_info.sid.has_value() && trim_info.sid.value() >= mgr.SidPoolValue().size()) [[unlikely]] {
    throw error::Error("SID {} from LUT exceeds SidPool size {}", trim_info.sid.value(), mgr.SidPoolValue().size());
  }

  ++metrics.preassign_passing_read_count;

  // Record per-read SID detection metrics before classifying the outcome; only mutated when
  // trim_info.sid is set. Passing-only length distributions and scalar full/partial counters are
  // recorded separately below.
  metrics.RecordReadMetrics(trim_info, raw_length);

  if (!trim_info.sid.has_value()) {
    // No SID assigned: the read is unassigned.
    ++metrics.unassigned_read_count;
    metrics.unassigned_base_count += raw_length;
    metrics.RecordRawFailedReadLength(raw_length);
    return;
  }
  const auto sid{trim_info.sid.value()};

  // Apply the post-assignment filters. A filtered read is a
  // failed-assigned read and is not written.
  if (SimplexPostAssignmentFilters(read, metrics, mgr, trim_info)) {
    return;
  }

  // Record full/partial counters here: the read has passed all gates, and these are passing-only.
  metrics.RecordCompleteness(trim_info);
  metrics.RecordPassingReadLength(trim_info, raw_length);

  // sid_index is the stable SID identity written to the FASTQ header. In kAllSplit mode,
  // partial reads are routed to a separate set of sinks at offset sid_pool.size() from
  // the normal (full) sinks, but retain the same header SID as full reads.
  read.sid_index = 1 + sid;
  const auto partial_offset =
      static_cast<u32>((mgr.Param().read_length_mode == ReadLengthMode::kAllSplit && IsPartialRead(trim_info))
                           ? mgr.SidPoolValue().size()
                           : 0);
  read.sink_id = read.sid_index + partial_offset;
  read.SetStatus(FixedReadRecord::Status::kDemultiplexed);

  // Passing read: count it only after routing is fully resolved above.
  metrics.IncrementSidPassingReadCount(sid);
  // total_bases are the trimmed (insert) bases of passing reads; trimmed_bases are what trimming removed.
  // Route both through AddSidPassingBases so the run-level and per-SID base totals stay in lockstep.
  const auto trimmed_length = static_cast<u64>(trim_info.insert.Length());
  metrics.AddSidPassingBases(sid, trimmed_length, raw_length - trimmed_length);
}

/**
 * @brief RunDuplexAdapter is a helper function to run the demuxing for duplex adapters.
 * @param read the read to be demultiplexed
 * @param duplex_metrics the DuplexMetrics instance to update the metrics
 * @param mgr the FlowManager instance to access the adapter and SID pool information
 **/
void Demux::RunDuplexAdapter(FixedReadRecord& read, DuplexMetrics& duplex_metrics, const FlowManager& mgr) {
  // actual demuxing done here
  mgr.HairpinFinderObj().FindHairpin(read, duplex_metrics);

  if (read.trim_info_duplex.duplex_status == TrimInfoDuplex::DuplexStatus::kMidAdapterFound) {
    read.sid_index = 1 + read.trim_info_duplex.matches[DuplexMatch::kSID3p].match.barcode_id;
    read.sink_id = read.sid_index;
    read.SetStatus(FixedReadRecord::Status::kDemultiplexed);
  } else {
    // all other states fail
    read.SetStatus(FixedReadRecord::Status::kDuplexMidAdapterFail);
    duplex_metrics.RecordRawFailedReadLength(read.SeqLen());
  }
}

bool Demux::PrefilterRawRead(FixedReadRecord& read,
                             SimplexMetrics& metrics,
                             DuplexMetrics& duplex_metrics,
                             const FlowManager& mgr) {
  if (mgr.IsDuplex()) {
    duplex_metrics.total_length_distr.AddCountToHistogram(read.SeqLen(), 1);
    if (read.SeqLen() < mgr.Param().min_read_len) {
      read.SetStatus(FixedReadRecord::Status::kTooShortFail);
      duplex_metrics.unassigned_counts.read_too_short += 1;
      duplex_metrics.unassigned_counts.raw_bases += read.SeqLen();
      duplex_metrics.unassigned_length_distr.AddCountToHistogram(read.SeqLen(), 1);
      duplex_metrics.RecordRawFailedReadLength(read.SeqLen());
      return true;
    }
    if (read.GetStatus() == FixedReadRecord::Status::kTooLongFail) {
      duplex_metrics.unassigned_counts.read_too_long += 1;
      duplex_metrics.unassigned_counts.raw_bases += read.SeqLen();
      duplex_metrics.unassigned_length_distr.AddCountToHistogram(read.SeqLen(), 1);
      return true;
    }
    return false;
  }

  metrics.IncrementInputReadCount();
  metrics.RecordInputReadLength(read.SeqLen());
  // Raw bases are accumulated over every input read (all total_reads), including reads
  // filtered below, so raw_bases is the percentage denominator for all base metrics.
  metrics.raw_base_count += read.SeqLen();
  if (read.SeqLen() < mgr.Param().min_read_len) {
    ++metrics.too_short_read_count;
    ++metrics.filtered_read_count;
    metrics.filtered_base_count += read.SeqLen();
    metrics.RecordRawFailedReadLength(read.SeqLen());
    read.SetStatus(FixedReadRecord::Status::kTooShortFail);
    // read is too short, skip read
    return true;
  }
  if (read.GetStatus() == FixedReadRecord::Status::kTooLongFail) {
    ++metrics.too_long_read_count;
    ++metrics.filtered_read_count;
    metrics.filtered_base_count += read.SeqLen();
    // read is too long, skip read
    return true;
  }
  return false;
}

void Demux::operator()() {
  try {
    StopWatch sw;
    const auto& mgr{context.GetManager()};
    auto& batch{context.GetBatchData(batch_nr)};
    auto& metrics{SimplexMetrics::Instance()};
    auto& duplex_metrics{DuplexMetrics::Instance()};
    if (batch.records) {
      for (size_t i = 0; i < batch.Size(); ++i) {
        auto& read = (*(batch.records))[i];
        // filter and record raw reads here
        if (PrefilterRawRead(read, metrics, duplex_metrics, mgr)) {
          continue;
        }

        // Demux and trim the read according to the adapter type
        switch (mgr.AdapterTypeValue()) {
          using enum AdapterType;
          case kYsu: {
            const auto raw_length = read.SeqLen();
            // actual demuxing done here
            read.trim_info_ysu = mgr.DemuxObjectYsu()(read);
            RunSimplexAdapter(read, metrics, raw_length, mgr, read.trim_info_ysu);
            break;
          }
          case kYs: {
            const auto raw_length = read.SeqLen();
            // actual demuxing done here
            read.trim_info_simplex = mgr.DemuxObjectYs()(read);
            RunSimplexAdapter(read, metrics, raw_length, mgr, read.trim_info_simplex);
            break;
          }
          case kSimplex: {
            const auto raw_length = read.SeqLen();
            RunSimplexAdapter(read, metrics, raw_length, mgr, mgr.DemuxObjectSimplex()(read));
            break;
          }
          case kSimplex10x: {
            const auto raw_length = read.SeqLen();
            read.trim_info_simplex = mgr.DemuxObjectSimplex10x()(read);
            RunSimplexAdapter(read, metrics, raw_length, mgr, read.trim_info_simplex);
            break;
          }
          case kDuplex:
          case kDuplexUMI:
          case kDuplexStem:
            RunDuplexAdapter(read, duplex_metrics, mgr);
            break;
          default:
            throw error::Error("Unknown adapter type found. Saw adapter type: {}",
                               std::to_underlying(mgr.AdapterTypeValue()));
        }
      }
    }

    context.AddToDemuxTime(sw.ElapsedTime());
  } catch (const std::exception& e) {
    Logging::Error("Demux::operator() failed: {}", e.what());
    SetTaskException(std::current_exception());
  }
}
}  // namespace xoos::demux
