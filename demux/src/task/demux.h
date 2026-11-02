#pragma once

#include "task/flow-manager.h"
#include "task/task.h"

namespace xoos::demux {
/// @brief Demux demuxes the data present in the input buffer and adds information about the barcodes.
class Demux : public Task {
 public:
  Demux(FlowContext& exec, size_t batch_nr);
  ~Demux() override = default;

  /// @brief This function is called by the task scheduler to execute the demux task.
  void operator()();

 private:
  /// @brief Apply the raw-read length pre-filters (too-short, too-long) before demultiplexing.
  /// Records the read against the relevant simplex or duplex metrics and sets its status.
  /// @return true if the read was filtered out (caller should skip it), false otherwise.
  static bool PrefilterRawRead(FixedReadRecord& read,
                               SimplexMetrics& metrics,
                               DuplexMetrics& duplex_metrics,
                               const FlowManager& mgr);

  template <typename TrimInfo>
  static void RunSimplexAdapter(FixedReadRecord& read,
                                SimplexMetrics& metrics,
                                u32 raw_length,
                                const FlowManager& mgr,
                                const TrimInfo& trim_info);

  /// @brief Apply the post-assignment filters (discordant SID, too-short-after-trim) to a
  /// SID-assigned read. Updates the read status and the corresponding failed-assigned metric.
  /// @return true if the read was filtered out (caller should stop processing it), false otherwise.
  template <typename TrimInfo>
  static bool SimplexPostAssignmentFilters(FixedReadRecord& read,
                                           SimplexMetrics& metrics,
                                           const FlowManager& mgr,
                                           const TrimInfo& trim_info);

  static void RunDuplexAdapter(FixedReadRecord& read, DuplexMetrics& duplex_metrics, const FlowManager& mgr);
};
}  // namespace xoos::demux
