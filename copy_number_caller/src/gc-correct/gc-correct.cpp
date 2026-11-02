#include "gc-correct/gc-correct.h"

#include <fstream>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/util/file-functions.h>

#include "cli/copy-number-caller-options.h"
#include "coverage/coverage.h"
#include "gc-correct/gc-correct-algorithm.h"
#include "intervals/augmented-intervals.h"
#include "io/copy-number-caller-default-filenames.h"

namespace xoos::cnc {
void GCCorrectMain(const CopyNumberCallerOptions& options) {
  const auto corrected_coverage_out = options.output_dir / kDefaultCorrectedCoverageOutput;
  file::CheckOutputFilePathsForWriting({corrected_coverage_out}, options.overwrite);
  Logging::Info("Loading and GC correcting coverage from {}", options.normal_coverage_fname.value().string());
  std::ifstream cov_istream(options.normal_coverage_fname.value());
  std::ifstream augmented_intervals_istream(options.augmented_intervals_fname.value());
  const CoverageRecords covs(cov_istream);
  const AugmentedIntervals augmented_intervals = AugmentedIntervals(augmented_intervals_istream);
  if (!VerifyRegions(covs.region, augmented_intervals.GetRegions())) {
    throw error::Error("Regions in coverage and intervals files do not match");
  }
  const GCCorrectResult res = GCCorrect(covs, augmented_intervals, options.gc_correct_options.first_span);
  CoverageRecords corrected_covs;
  corrected_covs.region = covs.region;
  corrected_covs.total_coverage = res.total_coverage;
  corrected_covs.average_coverage = res.average_coverage;
  corrected_covs.count = res.count;
  corrected_covs.on_target = covs.on_target;
  Logging::Info("Printing corrected coverage to {}", corrected_coverage_out.string());
  std::ofstream ofs(corrected_coverage_out);
  corrected_covs.Write(ofs, options.command_line_info);
}
}  // namespace xoos::cnc
