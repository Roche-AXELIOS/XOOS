#include "coverage/coverage.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <exception>
#include <fstream>
#include <numeric>
#include <sstream>
#include <tuple>

#include <csv.hpp>

#include <xoos/error/error.h>
#include <xoos/genome/region.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/util/missing-value.h>
#include <xoos/util/parse-float.h>
#include <xoos/util/parse-int.h>
#include <xoos/util/string-functions.h>

#include "io/column-names.h"
#include "util/eigen-math.h"
#include "util/parse-util.h"

namespace xoos::cnc {

/**
 * @brief load a coverage BED file.
 * Expected BED format: Contig\tStart\tEnd\tCounts\tTotalCoverage\tOnTarget (0-based half-open).
 * `#`-prefixed metadata/comments and a `#`-prefixed column header are supported and skipped.
 * @param ifs input stream for coverage BED file
 */
void CoverageRecords::LoadCoverageFile(std::istream& ifs) {
  std::string line;
  vec<f64> total_coverages;
  vec<f64> counts;
  while (std::getline(ifs, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::string contig;
    std::string start_str;
    std::string end_str;
    std::string total_cov_str;
    std::string count_str;
    std::string on_target_str;
    std::istringstream ss(line);
    if (!std::getline(ss, contig, '\t') || !std::getline(ss, start_str, '\t') || !std::getline(ss, end_str, '\t') ||
        !std::getline(ss, count_str, '\t') || !std::getline(ss, total_cov_str, '\t') ||
        !std::getline(ss, on_target_str, '\t')) {
      throw error::Error(
          "Coverage BED file must have 6 columns (Contig, Start, End, Counts, TotalCoverage, OnTarget): \"{}\"", line);
    }
    // reconstruct 1-based region string for internal use
    const auto [start, end, t, c] = [&]() -> std::tuple<u64, u64, f64, f64> {
      try {
        return {util::ParseU64(start_str),
                util::ParseU64(end_str),
                IsMissingValue(total_cov_str) ? 0 : util::ParseF64(total_cov_str),
                IsMissingValue(count_str) ? 0 : util::ParseF64(count_str)};
      } catch (const std::exception& e) {
        throw error::Error("Coverage BED file parse error: {} (line: \"{}\")", e.what(), line);
      }
    }();
    std::string r = contig + ":" + std::to_string(start + 1) + "-" + std::to_string(end);
    if (t < 0) {
      Logging::Error("Invalid total_coverage value {} for region {}", total_cov_str, r);
      throw error::Error("Invalid total_coverage value in coverage file");
    }
    if (c < 0) {
      Logging::Error("Invalid counts value {} for region {}", count_str, r);
      throw error::Error("Invalid counts value in coverage file");
    }
    region.push_back(r);
    total_coverages.push_back(t);
    counts.push_back(c);
    on_target.push_back(string::StringToBool(on_target_str));
  }
  total_coverage = eigen::FromStdVec(total_coverages);
  count = eigen::FromStdVec(counts);
  CalculateAverageCoverages();
  UpdateRegionMap();
}

/**
 * @brief repopulates `average_coverage`
 * @details repopulates the `average_coverage` member std::vector by dividing each
 * member in `total_coverage` by the genomic span derived from the corresponding
 * `region` element
 */
void CoverageRecords::CalculateAverageCoverages() {
  average_coverage = Eigen::VectorXd::Zero(total_coverage.size());
  for (size_t i = 0; i < region.size(); ++i) {
    const auto idx = static_cast<Eigen::Index>(i);
    const auto reg = region[i];
    const auto [contig, start, end] = ParseRegionString(reg);
    average_coverage[idx] = total_coverage[idx] / (static_cast<f64>(end) - static_cast<f64>(start));
  }
}

/**
 * @brief returns new CoverageRecords from a subset of row indexes
 * @param idxs - list of row indexes to subset
 * @return subsetted CoverageRecords
 */
CoverageRecords CoverageRecords::FilterRow(const vec<u64>& idxs) const {
  CoverageRecords new_covs;
  for (const auto i : idxs) {
    new_covs.region.push_back(region[i]);
    new_covs.on_target.push_back(on_target[i]);
  }
  new_covs.total_coverage = eigen::Elem(total_coverage, idxs);
  new_covs.average_coverage = eigen::Elem(average_coverage, idxs);
  new_covs.count = eigen::Elem(count, idxs);
  if (mean_mapping_quality.size() != 0) {
    new_covs.mean_mapping_quality = eigen::Elem(mean_mapping_quality, idxs);
  }
  new_covs.UpdateRegionMap();
  return new_covs;
}

/**
 * @brief returns new CoverageRecords from a subset containing given chromosome
 * @param chrom - chromsome to subset on
 * @return subsetted CoverageRecords
 */
CoverageRecords CoverageRecords::FilterChrom(const std::string& chrom) const {
  vec<u64> idxs;
  for (size_t i = 0; i < region.size(); ++i) {
    if (genome::IsEqualContig(region[i], chrom)) {
      idxs.push_back(i);
    }
  }
  return FilterRow(idxs);
}

/**
 * @brief returns new CoverageRecords from a subset containing the exact regions specified in regions_to_keep (this
 * function performs an lookup of the region string to retrieve the filtered subset; it does NOT try to look up based on
 * coordinate!
 * @param regions_to_keep - regions to subset on
 * @return subsetted CoverageRecords
 */
CoverageRecords CoverageRecords::FilterRegion(const vec<std::string>& regions_to_keep) const {
  vec<u64> idxs;
  for (const std::string& region : regions_to_keep) {
    const auto it = region_to_row.find(region);
    if (it == region_to_row.end()) {
      throw error::Error("Could not find {} in CoverageRecords", region);
    }
    idxs.push_back(it->second);
  }
  return FilterRow(idxs);
}

/**
 * @brief returns new CoverageRecords from on- or off-target regions
 * @param o - if true, then on-targets only. else off-targets.
 * @return subsetted CoverageRecords
 */
CoverageRecords CoverageRecords::GetOnOrOffTargetRegions(bool o) const {
  vec<u64> idxs;
  for (size_t i = 0; i < on_target.size(); ++i) {
    if ((o && on_target[i]) || (!o && !on_target[i])) {
      idxs.push_back(i);
    }
  }
  return FilterRow(idxs);
}

/**
 * @brief returns new CoverageRecords from regions containing allosomes
 * @return subsetted CoverageRecords
 */
CoverageRecords CoverageRecords::GetAllosomes() const {
  vec<u64> idxs;
  for (size_t i = 0; i < region.size(); ++i) {
    if (sex_predict::IsInAllosome(region[i])) {
      idxs.push_back(i);
    }
  }
  return FilterRow(idxs);
}

/**
 * @brief returns new CoverageRecords from regions containing autosomes
 * @return subsetted CoverageRecords
 */
CoverageRecords CoverageRecords::GetAutosomes() const {
  vec<u64> idxs;
  for (size_t i = 0; i < region.size(); ++i) {
    if (!sex_predict::IsInAllosome(region[i])) {
      idxs.push_back(i);
    }
  }
  return FilterRow(idxs);
}

/**
 * @brief print coverarge to file
 */
void CoverageRecords::Write(std::ofstream& ofs, const io::CommandLineInfo& command_line_info) const {
  auto writer = csv::make_tsv_writer(ofs);
  if (!command_line_info.version.empty()) {
    io::WriteTsvMetadata(ofs, command_line_info);
  }
  const vec<std::string> column_names{
      kColumnContigAsComment, kColumnStart, kColumnEnd, kColumnCounts, kColumnTotalCoverage, kColumnOnTarget};
  writer << column_names;
  const std::string na{util::kNA};
  for (size_t i = 0; i < region.size(); ++i) {
    const auto idx = static_cast<Eigen::Index>(i);
    const auto& [contig, start, end] = ParseRegionString(region[i]);
    const std::string total_cov_str = (total_coverage.size() != 0) ? std::to_string(total_coverage[idx]) : na;
    const std::string count_str = (count.size() != 0) ? std::to_string(count[idx]) : na;
    std::string on_target_str = na;
    if (!on_target.empty()) {
      on_target_str = on_target[i] ? "1" : "0";
    }
    writer << std::make_tuple(contig, start, end, count_str, total_cov_str, on_target_str);
  }
}

// mininum number of elements needed in std::vector to remove outliers
// these are implementation details which do not need to be exposed in the header file
const s32 kMinOutlierNelem = 5;
const f64 kIqrMinQuantile = 0.25;
const f64 kIqrMaxQuantile = 0.75;
const f64 kIqrCoeff = 1.5;

/**
 * @brief uses the inter-quantile range (IQR) method to remove outliers. anything outside the IQR (i.e. below the 25%
 * percentile or above the 75% percentile) will be deemed an outlier. If the span of the IQR times IQR_COEFF is less
 * than MIN_H, then the IQR range is too small and outliers will not be removed.
 * @param v std::vector of doubles from which to remove outliers
 * @return Eigen::VectorXd v with outliers removed
 */
namespace {
Eigen::VectorXd RemoveOutliers(const Eigen::VectorXd& v) {
  if (v.size() < kMinOutlierNelem) {
    return v;
  }
  Eigen::VectorXd probs(2);
  probs(0) = kIqrMinQuantile;
  probs(1) = kIqrMaxQuantile;
  const Eigen::VectorXd quants = eigen::Quantile(v, probs);
  const f64 iqr = quants[1] - quants[0];
  const f64 h = kIqrCoeff * iqr;
  const f64 lo = quants[0] - h;
  const f64 hi = quants[1] + h;
  Eigen::Matrix<bool, Eigen::Dynamic, 1> mask(v.size());
  for (Eigen::Index i = 0; i < v.size(); ++i) {
    mask(i) = (v[i] >= lo && v[i] <= hi);
  }
  return eigen::Elem(v, eigen::Find(mask));
}
}  // namespace

/**
 * @brief get chrX / chrY ratio. helpful for sex determination.
 * Supports both prefixed (chrX/chrY) and unprefixed (X/Y) contig names.
 */
XYRatioResult CoverageRecords::GetXYRatio() const {
  const CoverageRecords on_target_covs = GetOnOrOffTargetRegions(true);
  Eigen::VectorXd full_chr_x_cov = on_target_covs.FilterChrom("chrX").average_coverage;
  Eigen::VectorXd full_chr_y_cov = on_target_covs.FilterChrom("chrY").average_coverage;
  // Fall back to unprefixed contig names.
  if (full_chr_x_cov.size() == 0) {
    full_chr_x_cov = on_target_covs.FilterChrom("X").average_coverage;
  }
  if (full_chr_y_cov.size() == 0) {
    full_chr_y_cov = on_target_covs.FilterChrom("Y").average_coverage;
  }
  const Eigen::VectorXd chr_x_cov = RemoveOutliers(full_chr_x_cov);
  const Eigen::VectorXd chr_y_cov = RemoveOutliers(full_chr_y_cov);
  f64 chr_x_avg_cov = NAN;
  f64 chr_y_avg_cov = NAN;
  if (chr_x_cov.size() == 0 && chr_y_cov.size() == 0) {
    return {.xy_ratio = NAN, .chr_x_avg_cov = NAN, .chr_y_avg_cov = NAN};
  }
  // With few chrY observations the mean is unstable, so fall back to the median.
  constexpr u64 kMinChrYElemsForRatio = 10;
  if (static_cast<u64>(chr_y_cov.size()) < kMinChrYElemsForRatio) {
    chr_x_avg_cov = eigen::Median(chr_x_cov);
    chr_y_avg_cov = (chr_y_cov.size() != 0) ? eigen::Median(chr_y_cov) : 0;
  } else {
    chr_x_avg_cov = eigen::Mean(chr_x_cov);
    chr_y_avg_cov = (chr_y_cov.size() != 0) ? eigen::Mean(chr_y_cov) : 0;
  }
  const f64 xy_ratio = chr_x_avg_cov / (chr_y_avg_cov + 0.0001);
  return {.xy_ratio = xy_ratio, .chr_x_avg_cov = chr_x_avg_cov, .chr_y_avg_cov = chr_y_avg_cov};
}

/**
 * @brief Predict sex of a sample through ratio of chrX coverage over chrY coverage.
 * If the ratio > MIN_RATIO_FOR_SEX_DET, then return 'F'(Female)
 * if the ratio is between MIN_RATIO_FOR_SEX_DETR and MIN_RATIO_FOR_SEC_NA, then return 'N' (NA)
 * else return 'M' (Male)
 */
sex_predict::Sex CoverageRecords::PredictSex() const {
  using enum sex_predict::Sex;
  const XYRatioResult xy_ratio_result = GetXYRatio();
  const f64 xy_ratio = xy_ratio_result.xy_ratio;
  if (std::isnan(xy_ratio)) {
    return kUnknown;
  }
  if (xy_ratio > kMinRatioForSexDet) {
    return kFemale;
  }
  if (xy_ratio > kMinRatioForSexNA) {
    return kUnknown;
  }
  return kMale;
}

/**
 * @brief update region-to-row map
 */
void CoverageRecords::UpdateRegionMap() {
  region_to_row.clear();
  for (u64 i = 0; i < region.size(); ++i) {
    region_to_row[region[i]] = i;
  }
}

/**
 * @brief sort records by given indices
 */
void CoverageRecords::Sort(const vec<u64>& argsort) {
  vec<std::string> region_resort(region.size());
  vec<bool> on_target_resort(on_target.size());
  for (size_t i = 0; i < region.size(); ++i) {
    region_resort[i] = region[argsort[i]];
    on_target_resort[i] = on_target[argsort[i]];
  }
  region = region_resort;
  on_target = on_target_resort;
  total_coverage = eigen::Elem(total_coverage, argsort);
  if (average_coverage.size() != 0) {
    average_coverage = eigen::Elem(average_coverage, argsort);
  }
  count = eigen::Elem(count, argsort);
  UpdateRegionMap();
}

/**
 * @brief sort records by region
 */
void CoverageRecords::Sort() {
  vec<u64> argsort(region.size());
  // Seed with sequential indices [0, size) so regions can be sorted indirectly.
  std::ranges::iota(argsort, 0);
  std::ranges::stable_sort(argsort, [this](size_t x, size_t y) { return region[x] < region[y]; });
  Sort(argsort);
}

/**
 * @brief sort records by given regions
 */
void CoverageRecords::Sort(const vec<std::string>& ordered_regions) {
  vec<u64> argsort(region.size());
  size_t j = 0;
  for (const auto& r : ordered_regions) {
    auto it = region_to_row.find(r);
    if (it == region_to_row.end()) {
      Logging::Debug("{} does not appear in coverage object. skipping...", r);
      continue;
    }
    argsort[j++] = it->second;
  }
  if (j != region.size()) {
    throw error::Error("Intervals file is missing a region that was encountered in coverage file");
  }
  Sort(argsort);
}

bool VerifyRegions(const vec<std::string>& lhs, const vec<std::string>& rhs) {
  // we assume that the intervals and the coverage records have completely the same regions (no subsets)
  if (lhs.size() != rhs.size()) {
    return false;
  }
  for (size_t i = 0; i < lhs.size(); ++i) {
    if (lhs[i] != rhs[i]) {
      Logging::Error("region mismatch: interval {} vs coverage {} records", lhs[i], rhs[i]);
      return false;
    }
  }
  return true;
}

}  // namespace xoos::cnc
