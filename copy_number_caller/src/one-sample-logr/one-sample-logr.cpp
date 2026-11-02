#include "one-sample-logr/one-sample-logr.h"

#include <fstream>

#include <xoos/error/error.h>
#include <xoos/types/vec.h>
#include <xoos/util/file-functions.h>

#include "cli/copy-number-caller-options.h"
#include "coverage/coverage.h"
#include "io/column-names.h"
#include "io/copy-number-caller-default-filenames.h"
#include "observations/observations.h"
#include "one-sample-logr/self-normalize.h"

namespace xoos::cnc {
void OneSampleLogRMain(const CopyNumberCallerOptions& options) {
  const auto logrs_out = options.output_dir / kDefaultLogRsOutput;
  const auto logrs_bw_out = options.output_dir / kDefaultLogRsBwOutput;
  vec<fs::path> output_paths = {logrs_out};
  if (!options.reference_genome_fai_fname.empty()) {
    output_paths.push_back(logrs_bw_out);
  }
  file::CheckOutputFilePathsForWriting(output_paths, options.overwrite);
  std::ifstream cov_ifs(options.normal_coverage_fname.value());
  const CoverageRecords cov(cov_ifs);
  const Observations obvs(SelfNormalizeCounts(cov));

  std::ofstream ofs(logrs_out);
  if (!ofs.is_open()) {
    throw error::Error("Failed to open logR output file: {}", logrs_out.string());
  }
  obvs.Write(ofs, kColumnLogRatio, false, options.command_line_info);

  if (!options.reference_genome_fai_fname.empty()) {
    obvs.WriteBigWig(logrs_bw_out, options.reference_genome_fai_fname);
  }
}
}  // namespace xoos::cnc
