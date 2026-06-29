#include "task/reader.h"

#include <xoos/log/logging.h>

#include "fmt/format.h"
#include "task/flow-context.h"

namespace xoos::demux {

Reader::Reader(FlowContext& exec, size_t batch_nr) : Task(exec, batch_nr, fmt::format("reader{}", batch_nr)) {}

void Reader::operator()() {
  try {
    // Read input data - all the heavy lifting is actually done in the context class
    context.GetBatch(batch_nr);
  } catch (const std::exception& e) {
    Logging::Error("Error in Reader::operator(): {}", e.what());
    SetTaskException(std::current_exception());
  }
}
}  // namespace xoos::demux
