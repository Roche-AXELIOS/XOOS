#include "progress-meter.h"

namespace xoos::svc {

Progress::Progress(u64 total_region_count, u64 progress_update_percent)
    : total_region_count{total_region_count},
      region_count_per_progress_update{
          total_region_count > progress_update_percent ? total_region_count / progress_update_percent : 1} {
  // Constructor calculates the interval at which to log progress updates based on the total region count
  // and desired update percentage. Ensures at least 1 region per update to avoid division by zero.
}

void Progress::UpdateAndLog(const log::LogLevel level) {
  // Logging function for taskflow progress. Outputs a logging message only when a certain fraction of regions have been
  // processed.
  // Nothing to do if the total region count is zero, there is no progress to report and this
  // avoids a later division by zero.
  if (total_region_count == 0) {
    return;
  }

  // Atomically fetch and increment the count of regions processed, if we are a multiple of the progress update
  // then we log the progress.
  auto regions_processed_count = total_regions_processed_count.fetch_add(1);
  if (regions_processed_count == 0 || regions_processed_count % region_count_per_progress_update != 0) {
    return;
  }

  auto regions_processed_percent =
      static_cast<double>(regions_processed_count) / static_cast<double>(total_region_count) * 100;
  Logging::Log(level, "Processed '{:.2f}%'...", regions_processed_percent);
}

}  // namespace xoos::svc
