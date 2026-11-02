#include "copy-number-stats.h"

#include <cmath>

#include <xoos/types/float.h>

namespace xoos::stats {

namespace {
// Default copy number of a diploid (non-allosomal) locus in a normal cell.
constexpr f64 kDiploidCopyNumber = 2;
}  // namespace

/**
 * @brief expected log ratio between tumor and normal coverage given purity tumor ploidy and total copy number
 * @param purity
 * @param tumor_ploidy
 * @param total_copy_number
 * @return
 */
f64 ExpectedLogR(f64 purity, f64 tumor_ploidy, f64 total_copy_number) {
  const f64 expected_read_depth_ratio = ((kDiploidCopyNumber * (1 - purity)) + (purity * total_copy_number)) /
                                        ((kDiploidCopyNumber * (1 - purity)) + (purity * tumor_ploidy));
  return log2(expected_read_depth_ratio);
}

namespace {

/**
 * @brief expected log ratio between tumor and normal coverage given purity tumor ploidy and total copy number, assuming
 * that the default ploidy of the region is haploid
 * @param purity
 * @param tumor_ploidy
 * @param total_copy_number
 * @return
 */
f64 ExpectedLogRHaploid(const f64 purity, const f64 tumor_ploidy, const f64 total_copy_number) {
  const f64 expected_read_depth_ratio =
      ((1 - purity) + (purity * total_copy_number)) / ((1 - purity) + (purity * tumor_ploidy));
  return log2(expected_read_depth_ratio);
}

}  // namespace

/**
 * @brief expected log ratio between tumor and normal coverage given purity tumor ploidy and total copy number
 * @param purity
 * @param tumor_ploidy
 * @param total_copy_number
 * @param expect_haploid  - whether to treat the logrs as representing counts from a region with a default haploid state
 * (ex. a male allosome)
 * @return
 */
f64 ExpectedLogR(f64 purity, f64 tumor_ploidy, f64 total_copy_number, bool expect_haploid) {
  return expect_haploid ? ExpectedLogRHaploid(purity, tumor_ploidy, total_copy_number)
                        : ExpectedLogR(purity, tumor_ploidy, total_copy_number);
}

/**
 * @brief expected total copy number given purity, tumor ploidy and mean logr, and assuming that the default ploidy
 * of region is haploid
 * @param purity
 * @param tumor_ploidy
 * @param mean_logr
 * @return
 */
namespace {

f64 ExpectedTotalCopyNumberHaploid(const f64 purity, const f64 tumor_ploidy, const f64 mean_logr) {
  const f64 a = 1 - purity;
  const f64 numerator = (std::pow(2, mean_logr) * (a + (purity * tumor_ploidy))) - a;
  const f64 denominator = purity;
  return numerator / denominator;
}

/**
 * @brief expected total copy number given purity, tumor ploidy and mean logr
 * @param purity
 * @param tumor_ploidy
 * @param mean_logr
 * @return
 */
f64 ExpectedTotalCopyNumber(const f64 purity, const f64 tumor_ploidy, const f64 mean_logr) {
  const f64 a = 1 - purity;
  const f64 numerator =
      (std::pow(2, mean_logr) * ((kDiploidCopyNumber * a) + (purity * tumor_ploidy))) - (kDiploidCopyNumber * a);
  const f64 denominator = purity;
  return numerator / denominator;
}

}  // namespace

/**
 * @brief expected total copy number given purity, tumor ploidy and mean logr
 * @param purity
 * @param tumor_ploidy
 * @param mean_logr
 * @param expect_haploid  - whether to treat the logrs as representing counts from a region with a default haploid state
 * (ex. a male allosome)
 * @return
 */
f64 ExpectedTotalCopyNumber(f64 purity, f64 tumor_ploidy, f64 mean_logr, bool expect_haploid) {
  return expect_haploid ? ExpectedTotalCopyNumberHaploid(purity, tumor_ploidy, mean_logr)
                        : ExpectedTotalCopyNumber(purity, tumor_ploidy, mean_logr);
}

/**
 * @brief expected minor allele frequency of a germline heterozygous SNP given tumor purity, copy number and
 * minor-allele multiplicity
 * @param purity
 * @param copy_number
 * @param multiplicity
 * @return
 */
f64 ExpectedAF(f64 purity, f64 copy_number, f64 multiplicity) {
  return (purity * multiplicity + (1 - purity)) / (purity * copy_number + (kDiploidCopyNumber * (1 - purity)));
  // m/c if purity is 1
}

f64 ExpectedAF(f64 purity, f64 copy_number, f64 multiplicity, bool expect_haploid) {
  return expect_haploid ? NAN : ExpectedAF(purity, copy_number, multiplicity);
}

}  // namespace xoos::stats
