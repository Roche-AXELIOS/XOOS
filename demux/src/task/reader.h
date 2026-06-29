#pragma once
#include "task/task.h"

namespace xoos::demux {
/// @brief Reader task decodes data from the input file and writes it into a shared memory area
class Reader : public Task {
 public:
  Reader(FlowContext& exec, size_t batch_nr);
  ~Reader() override = default;

  /// @brief This function is called by the task scheduler to execute the reader task, which will read a batch
  /// of records from the input file and write them into the shared memory area.
  void operator()();
};
}  // namespace xoos::demux
