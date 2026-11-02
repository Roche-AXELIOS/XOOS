#include "two-sample-logr/calculate-logrs.h"

#include <Eigen/Dense>
#include <cmath>
#include <cstddef>

#include <xoos/error/error.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/util/math.h>

#include "coverage/coverage.h"
#include "observations/observations.h"
#include "util/eigen-math.h"
#include "util/parse-util.h"

namespace xoos::cnc {
void FilterAndNormalizeTumorNormal(CoverageRecords& tumor,
                                   CoverageRecords& normal,
                                   size_t normal_min_cov,
                                   size_t tumor_min_cov) {
  vec<u64> keep_idxs;
  const auto n = static_cast<size_t>(normal.average_coverage.size());
  for (size_t i = 0; i < n; ++i) {
    const auto idx = static_cast<Eigen::Index>(i);
    if (tumor.average_coverage[idx] >= static_cast<f64>(tumor_min_cov) &&
        normal.average_coverage[idx] >= static_cast<f64>(normal_min_cov)) {
      keep_idxs.push_back(i);
    }
  }
  if (keep_idxs.empty()) {
    throw error::Error("No regions passed coverage thresholds (normal >= {}, tumor >= {}); cannot normalize",
                       normal_min_cov,
                       tumor_min_cov);
  }
  normal = normal.FilterRow(keep_idxs);
  tumor = tumor.FilterRow(keep_idxs);
  const f64 normal_sum = eigen::Sum(normal.count);
  const f64 tumor_sum = eigen::Sum(tumor.count);
  if (math::IsCloseToZero(normal_sum) || math::IsCloseToZero(tumor_sum)) {
    throw error::Error(
        "Total count is zero after filtering (normal sum: {}, tumor sum: {}); cannot normalize", normal_sum, tumor_sum);
  }
  normal.count = normal.count / normal_sum;
  tumor.count = tumor.count / tumor_sum;
}

Observations CalculateLogrs(const CoverageRecords& tumor, const CoverageRecords& normal_normalized) {
  CoverageRecords tumor_normalized = tumor.FilterRegion(normal_normalized.region);
  Observations res;
  const auto n = static_cast<Eigen::Index>(normal_normalized.region.size());
  res.regions.resize(normal_normalized.region.size());
  res.contigs.resize(normal_normalized.region.size());
  res.starts.resize(n);
  res.ends.resize(n);
  res.obvs.resize(n);
  for (size_t i = 0; i < normal_normalized.region.size(); ++i) {
    if (normal_normalized.region[i] != tumor_normalized.region[i]) {
      throw error::Error("Tumor and normal regions do not match at index {}: '{}' vs '{}'",
                         i,
                         tumor_normalized.region[i],
                         normal_normalized.region[i]);
    }
    const auto& [contig, start, end] = ParseRegionString(normal_normalized.region[i]);
    const auto ei = static_cast<Eigen::Index>(i);
    res.contigs[i] = contig;
    res.starts[ei] = start;
    res.ends[ei] = end;
    res.regions[i] = normal_normalized.region[i];
  }
  // calculate log2-ratio; per-element std::log2.
  {
    const Eigen::VectorXd& tumor_count = tumor_normalized.count;
    const Eigen::VectorXd& normal_count = normal_normalized.count;
    res.obvs.resize(tumor_count.size());
    // NOSONAR(cpp:S5566) index loop over parallel vectors; per-element std::log2
    for (Eigen::Index i = 0; i < tumor_count.size(); ++i) {
      res.obvs[i] = std::log2(tumor_count[i]) - std::log2(normal_count[i]);
    }
  }

  // replace log2-ratio values that are -inf with the minimum log2-ratio value
  // -inf can happen if the tumour coverage is 0X
  const eigen::VectorXu64 finite_indices = eigen::FindFinite(res.obvs);
  const Eigen::VectorXd finite_values = eigen::Elem(res.obvs, finite_indices);

  const f64 min_finite = finite_values.minCoeff();
  for (f64& val : res.obvs) {
    if (std::isinf(val) && val < 0) {
      val = min_finite;
    }
  }
  // median-center the log2-ratios
  res.obvs = res.obvs.array() - eigen::Median(res.obvs);
  return res;
}
}  // namespace xoos::cnc
