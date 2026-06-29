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
  template <typename TrimInfo>
  static void RunSimplexAdapter(FixedReadRecord& read, SimplexMetrics& metrics, u32 raw_length, const FlowManager& mgr,
                                const TrimInfo& trim_info);

  static void RunDuplexAdapter(FixedReadRecord& read, DuplexMetrics& duplex_metrics, const FlowManager& mgr);
};
}  // namespace xoos::demux
