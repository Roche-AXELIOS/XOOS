#include "predict-sex/predict-sex.h"

#include <fstream>
#include <ostream>
#include <string>

#include <xoos/enum/enum-util.h>
#include <xoos/sex_predict/sex.h>
#include <xoos/util/file-functions.h>

#include "cli/copy-number-caller-options.h"
#include "coverage/coverage.h"

namespace xoos::cnc {

constexpr auto* kPredictSexOutputFilename = "predict-sex.tsv";

void PredictSexMain(const CopyNumberCallerOptions& options) {
  const auto output_fname = options.output_dir / kPredictSexOutputFilename;
  file::CheckOutputFilePathsForWriting({output_fname}, options.overwrite);
  std::ifstream cov_istream(options.predict_sex_coverage_fname);
  const CoverageRecords covs(cov_istream);

  const XYRatioResult xy_ratio_result = covs.GetXYRatio();
  const sex_predict::Sex sex = covs.PredictSex();
  const std::string sex_str = enum_util::FormatEnumName(sex);

  std::ofstream ofs(output_fname);
  ofs << "X_average_coverage" << "\t" << "Y_average_coverage" << "\t" << "XY_ratio" << "\t" << "Sex" << '\n';
  ofs << xy_ratio_result.chr_x_avg_cov << "\t" << xy_ratio_result.chr_y_avg_cov << "\t" << xy_ratio_result.xy_ratio
      << "\t" << sex_str << '\n';
}

}  // namespace xoos::cnc
