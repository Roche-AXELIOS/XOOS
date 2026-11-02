#pragma once
#include <Eigen/Dense>
#include <atomic>
#include <string>
#include <utility>

#include <xoos/io/metadata-util.h>
#include <xoos/sex_predict/sex.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>

namespace xoos::cnc {

constexpr f32 kMinRatioForSexDet = 25;
constexpr f32 kMinRatioForSexNA = 20;

/**
 * @struct XYRatioResult
 * @brief Represents the result of calculating the ratio of coverage between
 *        chromosome X and chromosome Y, along with their average coverages.
 *
 * @var XYRatioResult::ratio
 * The calculated ratio of average coverage between chromosome X and chromosome Y.
 *
 * @var XYRatioResult::chr_x_avg_cov
 * The average coverage of chromosome X.
 *
 * @var XYRatioResult::chr_y_avg_cov
 * The average coverage of chromosome Y.
 */
struct XYRatioResult {
  f64 xy_ratio;
  f64 chr_x_avg_cov;
  f64 chr_y_avg_cov;
};

/**
 * @class CoverageRecords
 * @brief container for records in a coverage file
 *
 * Input file should have the following format:
 * line 1: header
 * line 2-n: <region>\t<contig>\t<total_coverage>\t<count>...\n
 */
struct CoverageRecords {
  CoverageRecords() = default;

  explicit CoverageRecords(std::istream& ifs) {
    LoadCoverageFile(ifs);
  }

  // NOTE: we need to define custom move/copy constructors because the `_region_to_row` must be updated properly. The
  // only way we can do this without a segfault is to update this map manually (i.e. looping over each region and adding
  // an entry of the region, row_idx pair to the map).
  CoverageRecords(const CoverageRecords& rhs)
      : region(rhs.region),
        total_coverage(rhs.total_coverage),
        average_coverage(rhs.average_coverage),
        count(rhs.count),
        on_target(rhs.on_target),
        mean_mapping_quality(rhs.mean_mapping_quality),
        has_duplicates(rhs.has_duplicates.load()),
        has_rescued_secondaries(rhs.has_rescued_secondaries.load()) {
    UpdateRegionMap();
  }

  CoverageRecords(CoverageRecords&& rhs) noexcept
      : region(std::move(rhs.region)),
        total_coverage(std::move(rhs.total_coverage)),
        average_coverage(std::move(rhs.average_coverage)),
        count(std::move(rhs.count)),
        on_target(std::move(rhs.on_target)),
        mean_mapping_quality(std::move(rhs.mean_mapping_quality)),
        has_duplicates(rhs.has_duplicates.load()),
        has_rescued_secondaries(rhs.has_rescued_secondaries.load()) {
    UpdateRegionMap();
  }

  CoverageRecords& operator=(const CoverageRecords& rhs) {
    if (this == &rhs) {
      return *this;
    }
    region = rhs.region;
    total_coverage = rhs.total_coverage;
    average_coverage = rhs.average_coverage;
    count = rhs.count;
    on_target = rhs.on_target;
    mean_mapping_quality = rhs.mean_mapping_quality;
    has_duplicates = rhs.has_duplicates.load();
    has_rescued_secondaries = rhs.has_rescued_secondaries.load();
    UpdateRegionMap();
    return *this;
  }

  CoverageRecords& operator=(CoverageRecords&& rhs) noexcept {
    region = std::move(rhs.region);
    total_coverage = std::move(rhs.total_coverage);
    average_coverage = std::move(rhs.average_coverage);
    count = std::move(rhs.count);
    on_target = std::move(rhs.on_target);
    mean_mapping_quality = std::move(rhs.mean_mapping_quality);
    has_duplicates = rhs.has_duplicates.load();
    has_rescued_secondaries = rhs.has_rescued_secondaries.load();
    UpdateRegionMap();
    return *this;
  }

  ~CoverageRecords() noexcept = default;

  void LoadCoverageFile(std::istream& ifs);
  void CalculateAverageCoverages();

  void Write(std::ofstream& ofs, const io::CommandLineInfo& command_line_info) const;
  CoverageRecords FilterRegion(const vec<std::string>& /*regions_to_keep*/) const;
  CoverageRecords FilterRow(const vec<u64>& idxs) const;
  CoverageRecords FilterChrom(const std::string& chrom) const;
  CoverageRecords GetOnOrOffTargetRegions(bool /*o*/) const;
  CoverageRecords GetAutosomes() const;
  CoverageRecords GetAllosomes() const;
  XYRatioResult GetXYRatio() const;
  sex_predict::Sex PredictSex() const;
  void UpdateRegionMap();

  /// @brief whether there are no coverage records
  bool Empty() const {
    return region.empty();
  }

  void Sort(const vec<u64>& argsort);
  void Sort();
  void Sort(const vec<std::string>& ordered_regions);
  std::string sample_name;
  vec<std::string> region;
  Eigen::VectorXd total_coverage;
  Eigen::VectorXd average_coverage;
  Eigen::VectorXd count;
  vec<bool> on_target;
  Eigen::VectorXd mean_mapping_quality;
  std::atomic<bool> has_duplicates = false;
  std::atomic<bool> has_rescued_secondaries = false;
  StrUnorderedMap<u64> region_to_row;
};

bool VerifyRegions(const vec<std::string>& lhs, const vec<std::string>& rhs);
}  // namespace xoos::cnc
