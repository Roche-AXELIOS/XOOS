#pragma once

#include <atomic>

#include <xoos/log/logging.h>
#include <xoos/types/int.h>

namespace xoos::svc {

/**
 * @brief Progress tracking utility for monitoring completion of parallel region processing tasks.
 *        Provides thread-safe progress updates and periodic logging.
 */
struct Progress {
  const u64 total_region_count;
  const u64 region_count_per_progress_update;  // update progress log every 10%
  std::atomic_uint32_t total_regions_processed_count{0};

  /**
   * @brief Constructs a Progress tracker for monitoring region processing.
   * @param total_region_count Total number of regions to be processed
   * @param progress_update_percent Percentage interval for progress updates (default: 10%)
   */
  explicit Progress(u64 total_region_count, u64 progress_update_percent = 10);

  /**
   * @brief Thread-safe method to increment progress counter and log periodic updates.
   *        Logs progress only at specified percentage intervals to avoid spam.
   * @param level Log level for progress messages
   */
  void UpdateAndLog(log::LogLevel level);
};

}  // namespace xoos::svc
