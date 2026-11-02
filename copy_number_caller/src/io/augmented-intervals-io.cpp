#include "io/augmented-intervals-io.h"

#include <fstream>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/types/fs.h>

#include "cli/copy-number-caller-options.h"
#include "intervals/augment-intervals.h"
#include "intervals/augmented-intervals.h"

namespace xoos::cnc {

AugmentedIntervals LoadOrGenerateWholeGenomeIntervals(const CopyNumberCallerOptions& options,
                                                      const fs::path& augmented_intervals_out) {
  if (options.augmented_intervals_fname.has_value()) {
    std::ifstream ifs(options.augmented_intervals_fname.value());
    return AugmentedIntervals(ifs);
  }
  auto intervals = GenerateAndAugmentWholeGenomeIntervals(options);
  Logging::Info("Writing augmented intervals to: {}", augmented_intervals_out.string());
  std::ofstream ofs(augmented_intervals_out);
  if (!ofs.is_open()) {
    throw error::Error("Failed to open augmented intervals output file: {}", augmented_intervals_out.string());
  }
  intervals.Write(ofs, options.command_line_info);
  return intervals;
}

}  // namespace xoos::cnc
