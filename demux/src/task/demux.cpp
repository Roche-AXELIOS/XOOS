#include "task/demux.h"

#include <xoos/log/logging.h>

#include "adapters/duplex/trim-info-duplex.h"
#include "metrics/duplex-metrics.h"
#include "metrics/simplex-metrics.h"
#include "task/flow-context.h"
#include "utility/stop-watch.h"

namespace xoos::demux {

Demux::Demux(FlowContext& exec, size_t batch_nr) : Task(exec, batch_nr, fmt::format("Demux{}", batch_nr)) {}

/**
 * @brief RunSimplexAdapter is a helper function to run the demuxing and trimming for simplex adapters
 * @param read the read to be demultiplexed and trimmed
 * @param metrics the Metrics instance to update the metrics
 * @param raw_length
 * @param mgr the FlowManager instance to access the adapter and SID pool information
 * @param trim_info the TrimInfo object containing the results of the demuxing and trimming for the read
 * @return true if the read was successfully demultiplexed and trimmed, false if the read was filtered out due to being
 * too short after trimming
 **/
template <typename TrimInfo>
void Demux::RunSimplexAdapter(FixedReadRecord& read, SimplexMetrics& metrics, u32 raw_length, const FlowManager& mgr,
                              const TrimInfo& trim_info) {
  if (trim_info.insert.Length() < mgr.param.min_trimmed_read_len) {
    // read is too short after trimming, skip writing
    metrics.IncrementMinTrimmedReadLenFilteredCount();
    read.SetStatus(FixedReadRecord::Status::kTrimmedTooShortFail);
    return;
  }

  // Note: we already initialized the file_sid to 0, so by default we'll write to the 'none' file
  metrics.AddTrimmedRead(trim_info, raw_length);
  if (trim_info.sid) {
    // SID found, but verify whether it is valid
    if (const auto sid{trim_info.sid.value()}; mgr.sid_pool.contains(sid)) {
      // This is the file sid that we'll use to write the sample out
      read.file_sid = 1 + sid;
      read.SetStatus(FixedReadRecord::Status::kDemultiplexed);
    }
  }
}

/**
 * @brief RunDuplexAdapter is a helper function to run the demuxing for duplex adapters.
 * @param read the read to be demultiplexed
 * @param duplex_metrics the DuplexMetrics instance to update the metrics
 * @param mgr the FlowManager instance to access the adapter and SID pool information
 **/
void Demux::RunDuplexAdapter(FixedReadRecord& read, DuplexMetrics& duplex_metrics, const FlowManager& mgr) {
  // actual demuxing done here
  mgr.DemuxObjectDuplex().Demux(read, duplex_metrics);

  if (read.trim_info_duplex.duplex_status == TrimInfoDuplex::DuplexStatus::kMidAdapterFound) {
    read.file_sid = 1 + read.trim_info_duplex.matches[DuplexMatch::kSID3p].match.barcode_id;
    read.SetStatus(FixedReadRecord::Status::kDemultiplexed);
  } else {
    // all other states fail
    read.SetStatus(FixedReadRecord::Status::kDuplexMidAdapterFail);
  }
}

void Demux::operator()() {
  try {
    StopWatch sw;
    auto& mgr{context.GetManager()};
    auto& batch{context.GetBatchData(batch_nr)};
    auto& metrics{SimplexMetrics::Instance()};
    auto& duplex_metrics{DuplexMetrics::Instance()};
    if (batch.records) {
      for (size_t i = 0; i < batch.Size(); ++i) {
        auto& read = (*(batch.records))[i];
        // filter and record raw reads here
        if (mgr.is_duplex) {
          duplex_metrics.total_length_distr.AddCountToHistogram(read.SeqLen(), 1);
          if (read.SeqLen() < mgr.param.min_read_len) {
            read.SetStatus(FixedReadRecord::Status::kTooShortFail);
            duplex_metrics.unassigned_counts.read_too_short += 1;
            duplex_metrics.unassigned_length_distr.AddCountToHistogram(read.SeqLen(), 1);
            continue;
          }
          if (read.GetStatus() == FixedReadRecord::Status::kTooLongFail) {
            duplex_metrics.unassigned_counts.read_too_long += 1;
            duplex_metrics.unassigned_length_distr.AddCountToHistogram(read.SeqLen(), 1);
            continue;
          }
        } else {
          metrics.IncrementInputReadCount();
          if (read.SeqLen() < mgr.param.min_read_len) {
            read.SetStatus(FixedReadRecord::Status::kTooShortFail);
            // read is too short, skip read
            continue;
          }
          if (read.GetStatus() == FixedReadRecord::Status::kTooLongFail) {
            // read is too long, skip read
            continue;
          }
        }

        // Demux and trim the read according to the adapter type
        switch (mgr.adapter_type) {
          using enum AdapterType;
          case kYsuSl: {
            const auto raw_length = read.SeqLen();
            // actual demuxing done here
            read.trim_info_ysu_sl = mgr.DemuxObjectYsuSl()(read);
            RunSimplexAdapter(read, metrics, raw_length, mgr, read.trim_info_ysu_sl);
            break;
          }
          case kYsuTe: {
            const auto raw_length = read.SeqLen();
            // actual demuxing done here
            read.trim_info_ysu_te = mgr.DemuxObjectYsuTe()(read);
            RunSimplexAdapter(read, metrics, raw_length, mgr, read.trim_info_ysu_te);
            break;
          }
          case kYs: {
            const auto raw_length = read.SeqLen();
            // actual demuxing done here
            read.trim_info_ys = mgr.DemuxObjectYs()(read);
            RunSimplexAdapter(read, metrics, raw_length, mgr, read.trim_info_ys);
            break;
          }
          case kDuplex:
          case kDuplexUMI:
          case kDuplexStem:
            RunDuplexAdapter(read, duplex_metrics, mgr);
            break;
          default:
            throw error::Error("Unknown adapter type found. Saw adapter type: {}", static_cast<u32>(mgr.adapter_type));
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
