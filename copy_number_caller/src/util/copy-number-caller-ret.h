#pragma once
#include <optional>

#include <xoos/sex_predict/sex.h>
#include <xoos/types/vec.h>

#include "coverage/coverage-check.h"
#include "likelihood/total-copy-number-prior.h"
#include "observations/observations.h"
#include "seg-to-vcf/vcf-purity-source.h"
#include "segmentation/genomic-segments.h"

namespace xoos::cnc {
struct CopyNumberCallerRet {
  Observations ref_depths;
  Observations alt_depths;
  vec<GenomicSegment> segments;
  sex_predict::Sex sex{sex_predict::Sex::kUnknown};
  std::optional<f64> tumor_purity;
  std::optional<f64> tumor_ploidy;
  std::optional<VcfPuritySource> vcf_purity_source;
  std::optional<TotalCopyNumberPrior> total_copy_number_prior;
  std::optional<CoverageCheckResult> coverage_check_result;
};

}  // namespace xoos::cnc
