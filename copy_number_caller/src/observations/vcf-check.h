#pragma once
#include <optional>

#include <xoos/types/float.h>
#include <xoos/types/int.h>

namespace xoos::cnc {

/// Minimum fraction of heterozygous germline SNPs required in a somatic-flagged VCF.
constexpr f64 kDefaultMinHetSnpFraction = 0.5;

struct VcfCheckResult {
  u32 num_het_snps;
  u32 num_somatic_variants;
  // std::nullopt when denominator is zero
  std::optional<f64> het_snp_fraction;
  bool is_sufficient;
};

/// Returns true when a VCF check was performed and the result was insufficient.
inline bool IsVcfInsufficient(const std::optional<VcfCheckResult>& result) {
  return result.has_value() && !result->is_sufficient;
}

}  // namespace xoos::cnc
