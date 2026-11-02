#pragma once
#include <Eigen/Dense>
#include <memory>
#include <string>
#include <utility>

#include <xoos/io/metadata-util.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "cli/sample-metadata-options.h"
#include "cli/vcf-parsing-options.h"
#include "observations/vcf-check.h"
#include "util/eigen-math.h"

namespace xoos::cnc {

// A biallelic SNV has a REF and an ALT allele, so each sample occupies two consecutive AD entries.
constexpr size_t kBiallelicAlleleCount = 2;

// Data-carrier struct: the coordinate/value vectors are intentionally public for direct access by the
// segmentation and likelihood code, while `_is_baf_seg_obs` is kept private behind accessors because it is
// an invariant flag rather than raw data. The resulting public/private data-member mix is deliberate.
struct Observations {  // NOSONAR(cpp:S5414,cpp:M23_170) — intentional data-carrier design (see note above)
  Observations() = default;

  // Convenience constructor from std::vector coordinates/values. Enables brace-init call sites
  // (e.g. Observations({"chr1"}, {0}, {100}, {0.5})). Conversion is element-order-preserving and exact.
  // NOLINTNEXTLINE — used by tests only; kept as the ergonomic brace-init entry point.  NOSONAR
  Observations(vec<std::string> c, const vec<u64>& s, const vec<u64>& e, const vec<f64>& o)
      : contigs(std::move(c)),
        regions({}),
        starts(eigen::FromStdVecU64(s)),
        ends(eigen::FromStdVecU64(e)),
        obvs(eigen::FromStdVec(o)),
        dps({}) {
  }

  explicit Observations(size_t size) {
    const auto n = static_cast<Eigen::Index>(size);
    contigs.resize(size);
    regions.resize(size);
    // Eigen resize does not zero-fill; use setZero to zero-initialize.
    starts = eigen::VectorXu64::Zero(n);
    ends = eigen::VectorXu64::Zero(n);
    obvs = Eigen::VectorXd::Zero(n);
    dps.resize(size);
  }

  Observations(const Observations& rhs) = default;
  Observations(Observations&& rhs) noexcept = default;
  Observations& operator=(const Observations& rhs) = default;
  Observations& operator=(Observations&& rhs) noexcept = default;
  ~Observations() = default;

  void SetBAFSegObvsStatus(const bool status) {
    _is_baf_seg_obs = status;
  }

  void PopulateFieldsFromRegions();

  bool IsBAFSegObvs() const {
    return _is_baf_seg_obs;
  }

  bool IsSorted() const;

  // NOLINTBEGIN(misc-non-private-member-variables-in-classes,cppcoreguidelines-non-private-member-variables-in-classes)
  vec<std::string> contigs;
  vec<std::string> regions;
  eigen::VectorXu64 starts;
  eigen::VectorXu64 ends;
  Eigen::VectorXd obvs;
  vec<s32> dps;
  // NOLINTEND(misc-non-private-member-variables-in-classes,cppcoreguidelines-non-private-member-variables-in-classes)
  Observations FilterByRange(size_t i, size_t j) const;
  Observations FilterByIdxs(const vec<size_t>& idxs) const;

  void Write(std::ofstream& ofs,
             const std::string& value_column_name,
             bool value_as_int,
             const io::CommandLineInfo& command_line_info) const;
  void WriteBigWig(const fs::path& filename, const fs::path& fai_filename) const;

 private:
  bool _is_baf_seg_obs{false};
};

struct RefAltObservations {
  Observations ref_obvs;
  Observations alt_obvs;

  // optional vector of somatic VAFs, only populated if the VCF is from a tumor-normal pair and contains somatic
  // variants.
  std::optional<vec<f64>> somatic_vafs;

  // VCF data-quality check result, populated for tumor-normal VCF parsing.
  std::optional<VcfCheckResult> vcf_check;
};

struct VcfHeaderInfo {
  std::string normal_sample;
  std::string tumour_sample;
  bool has_tumor_normal{false};
  s32 normal_index{-1};
  s32 tumor_index{-1};
};

Observations ReadObservations(std::istream& in_stream);
Observations GetDhFromDepths(const Observations& ref_depths, const Observations& alt_depths);
Observations GetDhFromVcf(io::VcfReader& vcf_rdr,
                          bool keep_all_variants,
                          bool has_tumor_normal,
                          const BafFilterOptions& baf_filter_options,
                          const SampleMetadataOptions& sample_metadata_options);
RefAltObservations GetDepthsFromVcf(io::VcfReader& vcf_rdr,
                                    bool keep_all_variants,
                                    bool has_tumor_normal,
                                    const BafFilterOptions& baf_filter_options,
                                    const SampleMetadataOptions& sample_metadata_options);

VcfHeaderInfo GetTumorNormalSampleInfoFromHeader(const std::shared_ptr<io::VcfHeader>& header,
                                                 const std::string& normal_sample_name,
                                                 const std::string& tumor_sample_name);
VcfHeaderInfo GetSingleSampleInfoFromHeader(const std::shared_ptr<io::VcfHeader>& header,
                                            const std::string& normal_sample_name);
RefAltObservations GetDepthsFromFile(const fs::path& fname);

constexpr f64 kExtremeBafLowThreshold = 0.03;
constexpr f64 kExtremeBafHighThreshold = 0.97;
constexpr f64 kExtremeBafProportionThreshold = 0.10;
constexpr f64 kExtremeBafPurityThreshold = 0.80;

/*
 * @brief Compute the proportion of heterozygous SNPs (based on normal) with extreme tumor BAF (> 0.97 or < 0.03).
 * @param ref_depths Reference allele depths from the tumor
 * @param alt_depths Alternate allele depths from the tumor
 * @return Proportion of SNPs with extreme BAF, or 0.0 if there are no valid depth loci
 * @throws error::Error if ref_depths and alt_depths have mismatched sizes
 */
f64 ComputeExtremeBafProportion(const Observations& ref_depths, const Observations& alt_depths);

/*
 * @brief Check if the extreme BAF proportion is high and log a warning if so.
 * @param ref_depths Reference allele depths from the tumor
 * @param alt_depths Alternate allele depths from the tumor
 * @param purity Estimated tumor purity
 * @return True if the extreme BAF proportion is high, false otherwise
 */
bool CheckExtremeBafProportionAndWarn(const Observations& ref_depths, const Observations& alt_depths, f64 purity);
}  // namespace xoos::cnc
