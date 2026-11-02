#include "two-sample-logr/two-sample-logr.h"

#include <fstream>
#include <stdexcept>

#include <fmt/format.h>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/types/vec.h>
#include <xoos/util/file-functions.h>

#include "cli/copy-number-caller-cli-option-names.h"
#include "cli/copy-number-caller-options.h"
#include "coverage/coverage.h"
#include "io/column-names.h"
#include "io/copy-number-caller-default-filenames.h"
#include "observations/observations.h"
#include "two-sample-logr/tumor-normal.h"

namespace xoos::cnc {

namespace opt = cli_opt_name;

void TwoSampleLogRMain(const CopyNumberCallerOptions& options) {
  const auto logrs_out = options.output_dir / kDefaultLogRsOutput;
  const auto logrs_bw_out = options.output_dir / kDefaultLogRsBwOutput;
  vec<fs::path> output_paths = {logrs_out};
  if (!options.reference_genome_fai_fname.empty()) {
    output_paths.push_back(logrs_bw_out);
  }
  file::CheckOutputFilePathsForWriting(output_paths, options.overwrite);
  if (!options.tumor_coverage_fname.has_value() || !options.normal_coverage_fname.has_value()) {
    Logging::Error(fmt::format("Must provide {} and {} options!", opt::kTumorCoverageFile, opt::kNormalCoverageFile));
    throw std::runtime_error("missing option");
  }
  std::ifstream tumor_cov_ifs(options.tumor_coverage_fname.value());
  const CoverageRecords tumor_cov(tumor_cov_ifs);
  std::ifstream normal_cov_ifs(options.normal_coverage_fname.value());
  const CoverageRecords normal_cov(normal_cov_ifs);
  const Observations obvs(
      ProcessTumorNormal(tumor_cov, normal_cov, options.normal_min_coverage, options.tumor_min_coverage));

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
