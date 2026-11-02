#include "one-sample-logr/self-normalize.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include <xoos/log/logging.h>
#include <xoos/sex_predict/sex.h>
#include <xoos/types/float.h>

#include "coverage/coverage.h"
#include "observations/observations.h"
#include "util/eigen-math.h"
#include "util/parse-util.h"

namespace xoos::cnc {

// Expected logR at total-copy-number 0 for a near-pure germline sample (purity 0.99, ploidy 2):
// log2 of the ratio of expected coverage at TCN 0 to expected coverage at the ploidy baseline.
const f64 kExpectedLogRForGermlineAtTCN0 = -6.643856;

/**
 * @brief self-normalize coverage counts by dividing by the median count on autosomes, and log2-transforming
 * @param coverage coverage records with count field populated
 * @return Observations object with logR values
 */
Observations SelfNormalizeCounts(const CoverageRecords& coverage) {
  // Collect the autosome counts by iterating over the records rather than
  // materializing a full autosome copy of CoverageRecords (region strings, all
  // Eigen columns, and the region->row map) just to read the count column.
  // Count autosome records first so the buffer is sized exactly, avoiding an
  // over-allocation to the full record count and a shrinking resize.
  const auto is_autosome = [](const std::string& region) { return !sex_predict::IsInAllosome(region); };
  const auto n_autosome = std::ranges::count_if(coverage.region, is_autosome);
  Eigen::VectorXd autosome_counts(n_autosome);
  Eigen::Index next = 0;
  for (size_t i = 0; i < coverage.region.size(); ++i) {
    if (is_autosome(coverage.region[i])) {
      autosome_counts[next] = coverage.count[static_cast<Eigen::Index>(i)];
      ++next;
    }
  }
  const f64 median_count = eigen::Median(autosome_counts);
  if (median_count <= 0) {
    Logging::Error("median count for corrected coverage is 0!");
    throw std::runtime_error("no coverage!!");
  }
  const auto n = static_cast<size_t>(coverage.count.size());
  Observations ret(n);
  for (size_t i = 0; i < n; ++i) {
    const auto eidx = static_cast<Eigen::Index>(i);
    const auto [contig, start, end] = ParseRegionString(coverage.region[i]);
    ret.regions[i] = coverage.region[i];
    ret.contigs[i] = contig;
    ret.starts[eidx] = start;
    ret.ends[eidx] = end;
    if (coverage.count[eidx] > 0) {
      ret.obvs[eidx] = std::log2(coverage.count[eidx]) - std::log2(median_count);
      //  targets with 1 count can have logR value lower than the minimum we should expect; so just truncate these
      ret.obvs[eidx] = std::max(ret.obvs[eidx], kExpectedLogRForGermlineAtTCN0);
    } else {
      ret.obvs[eidx] = kExpectedLogRForGermlineAtTCN0;
    }
  }
  return ret;
}

}  // namespace xoos::cnc
